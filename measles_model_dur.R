## =====================================================================
## measles_model_dur.R
## Same as measles_model.R, but with the NEW stage durations
## (log-normal latent period, gamma prodromal period; patched package).
##
## measles_model.R itself is NOT changed: this file loads the patched
## package, then sources measles_model.R, then
##   - makes every model it builds use the new durations, and
##   - writes results/figures to results_dur/ and figures_dur/.
##
## Use:  source("measles_model_dur.R")   instead of   source("measles_model.R")
## Needs the patched package in MEASLES_LIB (default ~/R_libs_dur).
## =====================================================================

## 1. Load the patched package FIRST, so measles_model.R cannot load the original
local({
  lib  <- path.expand(Sys.getenv("MEASLES_LIB", unset = "~/R_libs_dur"))
  libs <- c(lib, .libPaths(), path.expand("~/R_libs_crn"), path.expand("~/R_libs"))
  .libPaths(unique(libs[dir.exists(libs)]))
  for (p in c("epiworldR", "measles")) suppressPackageStartupMessages(
    if (dir.exists(file.path(lib, p))) library(p, lib.loc = lib, character.only = TRUE)
    else library(p, character.only = TRUE))
})

## Stop if the COMPILED package is not the patched one
local({
  so <- file.path(find.package("measles"), "libs", "measles.so")
  if (!length(grepRaw("Latent dispersion", readBin(so, "raw", file.info(so)$size), fixed = TRUE)))
    stop("The loaded measles package is not the patched one: ", so,
         "\nInstall it into ~/R_libs_dur (run_dur.sh does this) or set MEASLES_LIB.")
})

## 2. Everything else exactly as in measles_model.R
##    (data.table / ggplot2 are loaded first and the library folders put back
##    afterwards, because measles_model.R switches to ~/R_libs, whose compiled
##    packages need a newer compiler than the job has)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
.keep_libs <- .libPaths()
source("measles_model.R")
.libPaths(.keep_libs)

## 2a. Optional: every school with the same number of students, e.g.
##     SCHOOL_SIZE=400. Grade sizes are scaled in proportion to ENROLLMENT
##     (rounded so they add up exactly), so the share of students in each
##     grade - and the per-capita contact matrix - stay as measured.
SCHOOL_SIZE <- as.integer(Sys.getenv("SCHOOL_SIZE", unset = 0))
if (SCHOOL_SIZE > 0) {
  scale_to <- function(n, total) {
    x <- n / sum(n) * total; f <- floor(x)
    extra <- total - sum(f)
    if (extra > 0) { i <- order(x - f, decreasing = TRUE)[seq_len(extra)]; f[i] <- f[i] + 1 }
    setNames(as.integer(f), names(n))
  }
  ENROLLMENT <- lapply(ENROLLMENT, scale_to, total = SCHOOL_SIZE)
  if (exists("CONTACT_RAW") && exists("symmetrise_recip"))
    SCHOOL_CMAT <- setNames(lapply(SCHOOLS, function(s)
      symmetrise_recip(CONTACT_RAW[[s]], ENROLLMENT[[s]])), SCHOOLS)
  for (s in SCHOOLS)
    if (any(diag(SCHOOL_CMAT[[s]]) > ENROLLMENT[[s]]))
      stop("SCHOOL_SIZE=", SCHOOL_SIZE, " is too small for ", s,
           ": a grade would have fewer students than its same-grade contacts per day")
  TRANS_RATE <- vapply(SCHOOL_CMAT, function(S) R0 / (INF_DAYS * max(abs(eigen(S)$values))), numeric(1))
}

## 2b. Optional: a different R0 (default 12, as in measles_model.R).
##     R0=7 or R0=9 recomputes the transmission rate of every school the same
##     way measles_model.R does: beta = R0 / (INF_DAYS x largest eigenvalue).
R0_ENV <- Sys.getenv("R0", unset = "")
if (nzchar(R0_ENV)) {
  R0 <- as.numeric(R0_ENV)
  if (is.na(R0) || R0 <= 0) stop("R0 must be a positive number, e.g. R0=9")
  TRANS_RATE <- vapply(SCHOOL_CMAT,
                       function(S) R0 / (INF_DAYS * max(abs(eigen(S)$values))),
                       numeric(1))
}

## 3. Optional: contagious rash. RASH_DAYS=4 makes a student with an unnoticed
##    rash who is still at school contagious for the first 4 rash days
##    (needs the rash-patched package). 0 (default) = rash not contagious.
RASH_DAYS <- as.numeric(Sys.getenv("RASH_DAYS", unset = 0))
if (RASH_DAYS > 0) {
  so <- file.path(find.package("measles"), "libs", "measles.so")
  if (!length(grepRaw("Rash contagious days", readBin(so, "raw", file.info(so)$size), fixed = TRUE)))
    stop("RASH_DAYS > 0 but the compiled measles package has no contagious-rash option: ", so)
  ## Recalibrate beta: with RASH_DAYS = 4 each case is contagious at school for
  ## about 5.25 days instead of 4 (checked: secondary cases from the first case
  ## then match the rash-not-contagious model). Other RASH_DAYS: re-check.
  INF_DAYS_RASH <- if (RASH_DAYS == 4) 5.25 else
    stop("INF_DAYS_RASH is only calibrated for RASH_DAYS = 4")
  TRANS_RATE <- vapply(SCHOOL_CMAT,
                       function(S) R0 / (INF_DAYS_RASH * max(abs(eigen(S)$values))),
                       numeric(1))
}

## 4. Every model built from here on uses the new durations (and the rash setting)
ModelMeaslesMixingRiskQuarantine <- function(...) {
  m <- measles::ModelMeaslesMixingRiskQuarantine(...)
  set_param(m, "Duration model", 1)
  if (RASH_DAYS > 0) set_param(m, "Rash contagious days", RASH_DAYS)
  m
}

## 5. Keep these results apart from the original ones
##    (R0 other than 12 gets its own folders too, e.g. results_dur_R0_9)
TAG <- paste0(if (RASH_DAYS > 0) "_rash" else "", if (R0 != 12) sprintf("_R0_%g", R0) else "",
              if (SCHOOL_SIZE > 0) sprintf("_size%d", SCHOOL_SIZE) else "")
RESULTS_DIR <- paste0("results_dur", TAG)
FIG_DIR     <- paste0("figures_dur", TAG)

cat("measles from:", find.package("measles"), "| Duration model: 1 (new durations)",
    "| Rash contagious days:", RASH_DAYS, "| R0:", R0, "| results in:", RESULTS_DIR, "\n")
cat("Students per school:", paste(names(ENROLLMENT), vapply(ENROLLMENT, sum, 0), collapse = ", "), "\n")
for (s in names(ENROLLMENT)) cat("  ", s, ":", paste(names(ENROLLMENT[[s]]), ENROLLMENT[[s]], sep = "=", collapse = " "), "\n")
