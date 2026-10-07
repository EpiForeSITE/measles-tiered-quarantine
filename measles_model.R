## =====================================================================
## measles_model.R
##
## Shared code for the measles tiered-quarantine analyses.
## Every .qmd in this folder starts with   source("measles_model.R")
##
##   1. Packages and global settings
##   2. School contact matrices, enrollment, transmission rates
##   3. simulator()       one scenario -> vector of final outbreak sizes
##   4. Design helpers    strategies, vaccination grids, run_or_load()
##   5. Small utilities   Wilson CI, plot theme, save_fig()
## =====================================================================


## ---------------------------------------------------------------------
## 1. Packages and global settings
## ---------------------------------------------------------------------

## Personal library used on the cluster (skipped if it doesn't exist)
local({
  lib <- file.path(Sys.getenv("HOME"), "R_libs")
  if (dir.exists(lib)) .libPaths(lib)
})

suppressPackageStartupMessages({
  library(epiworldR)
  library(measles)
  library(data.table)
  library(ggplot2)
})

R0       <- 12     # target basic reproduction number
INF_DAYS <- 4      # infectious (prodromal) window, days
N_DAYS   <- 100    # simulated days per run
SEED     <- 221

## Simulations per scenario. Override for a quick test, e.g.
##   N_SIMS=200 quarto render 01_tiered_quarantine_tables.qmd
N_SIMS <- as.integer(Sys.getenv("N_SIMS", unset = 4000))

## Threads: use the SLURM allocation if there is one
N_THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                   unset = parallel::detectCores()))

SCHOOLS <- c("Elementary", "Middle", "High")

RESULTS_DIR <- "results"   # cached simulations (CSV) and text output
FIG_DIR     <- "figures"   # PNG copies of every figure


## ---------------------------------------------------------------------
## 2. Schools: real per-capita SMART contact matrices
##    (rows/cols = grades; symmetrised for reciprocity, real magnitude)
## ---------------------------------------------------------------------

grade_matrix <- function(values, grades) {
  matrix(values, nrow = length(grades), byrow = TRUE,
         dimnames = list(grades, grades))
}

CONTACT_RAW <- list(
  Elementary = grade_matrix(c(
    18.74, 1.32, 0.25, 0.09, 3.77, 0.14,
     1.14,24.52, 3.06, 0.02, 0.17, 0.00,
     0.29, 4.21,13.96, 6.60, 1.60, 0.00,
     0.08, 0.02, 5.03,17.49, 9.57, 0.02,
     3.21, 0.16, 1.15, 9.00,25.37, 0.79,
     0.21, 0.00, 0.00, 0.03, 1.39,32.42),
    grades = c("K", "1", "2", "3", "4", "5")),

  Middle = grade_matrix(c(
    35.15, 0.67, 0.00, 0.00,
     0.61,27.91, 0.00, 0.00,
     0.00, 0.00,65.94, 0.11,
     0.00, 0.00, 0.17,31.32),
    grades = c("5", "6", "7", "8")),

  High = grade_matrix(c(
    25.98,14.45, 3.94, 1.38,
    12.12,20.71, 6.04, 2.42,
     4.56, 8.33,14.26, 6.31,
     1.89, 3.97, 7.48,14.97),
    grades = c("9", "10", "11", "12"))
)

symmetrise  <- function(C) 0.5 * (C + t(C))
SCHOOL_CMAT <- lapply(CONTACT_RAW, symmetrise)

## Students per grade, and enrollment vector for each school
PER_GRADE  <- c(Elementary = 26L, Middle = 59L, High = 58L)
ENROLLMENT <- sapply(SCHOOLS, function(s)
  rep(PER_GRADE[[s]], nrow(SCHOOL_CMAT[[s]])), simplify = FALSE)

## R0 is hit through a per-school transmission rate (matrix is NOT rescaled):
##   p = R0 / (INF_DAYS * largest eigenvalue)
TRANS_RATE <- vapply(SCHOOL_CMAT,
                     function(S) R0 / (INF_DAYS * max(abs(eigen(S)$values))),
                     numeric(1))


## ---------------------------------------------------------------------
## 3. Simulator
##    duration = c(high, medium, low) quarantine days
##    returns one final outbreak size per simulation
##    (everyone who is no longer susceptible at the last day)
## ---------------------------------------------------------------------

simulator <- function(school, duration, vaccinated, nsims = N_SIMS, seed = SEED) {

  sizes <- ENROLLMENT[[school]]
  n     <- sum(sizes)

  m <- ModelMeaslesMixingRiskQuarantine(
    n                           = n,
    prevalence                  = 1 / n,          # one index case
    contact_matrix              = SCHOOL_CMAT[[school]],
    transmission_rate           = TRANS_RATE[[school]],
    prop_vaccinated             = vaccinated,
    detection_rate_quarantine   = 0.0,
    contact_tracing_days_window = 7,
    quarantine_period_high      = as.integer(duration[1]),
    quarantine_period_medium    = as.integer(duration[2]),
    quarantine_period_low       = as.integer(duration[3])
  )

  for (k in seq_along(sizes))
    m <- add_entity(m, entity("Grade", as.integer(sizes[k]), FALSE))

  run_multiple(m, ndays = N_DAYS, nsims = nsims, seed = seed,
               saver = make_saver("total_hist"),
               nthreads = N_THREADS, verbose = FALSE)

  ans <- run_multiple_get_results(m, freader = data.table::fread,
                                  nthreads = 1L)$total_hist

  ans[date == max(date)][
    state != "Susceptible" & state != "Susceptible Quarantine",
    .(total = sum(counts)), by = sim_num]$total
}


## ---------------------------------------------------------------------
## 4. Design helpers
## ---------------------------------------------------------------------

## Strategies written as "H/M/L" days -> table with q_high, q_med, q_low.
## Names (optional) become the scenario label.
##   strategy_table(c("No quarantine" = "0/0/0", "21/7/7"))
strategy_table <- function(strats) {
  labs <- if (is.null(names(strats))) strats else
    ifelse(names(strats) == "", strats, names(strats))
  d <- lapply(tstrsplit(strats, "/", fixed = TRUE), as.integer)
  data.table(strategy = labs, q_high = d[[1]], q_med = d[[2]], q_low = d[[3]])
}

## Every row of `design` at every vaccination level
cross_vacc <- function(design, vaccs) {
  out <- design[rep(seq_len(nrow(design)), times = length(vaccs))]
  out[, vacc := rep(vaccs, each = nrow(design))]
  out[]
}

## Full grid of quarantine days, keeping only high >= medium >= low
policy_grid <- function(vaccs, highs, meds, lows) {
  g <- as.data.table(expand.grid(vacc = vaccs, q_high = highs,
                                 q_med = meds, q_low = lows))
  g <- g[q_high >= q_med & q_med >= q_low]
  setorder(g, vacc, q_high, q_med, q_low)
  g[]
}

## Run every design row in every school. Returns one row per simulation.
##
## By default every scenario uses the same seed (SEED). Run k of one
## scenario then starts from the same random numbers as run k of every
## other scenario, which makes comparisons between scenarios less noisy.
##
## independent_seeds = TRUE gives every school x scenario its own seed,
## so scenarios are independent samples. Use this whenever a statistical
## test that assumes independence (e.g. Wilcoxon rank-sum) is applied.
simulate_design <- function(design, schools = SCHOOLS, independent_seeds = FALSE) {
  rbindlist(lapply(seq_along(schools), function(j) {
    s <- schools[j]
    rbindlist(lapply(seq_len(nrow(design)), function(i) {
      d    <- design[i]
      seed <- if (independent_seeds) SEED + 1000L * j + i else SEED
      message(sprintf("%-10s | %3d/%d | vacc %3.0f%% | H/M/L = %d/%d/%d | seed %d",
                      s, i, nrow(design), 100 * d$vacc,
                      d$q_high, d$q_med, d$q_low, seed))
      data.table(school = s, d, seed = seed,
                 size = simulator(s, c(d$q_high, d$q_med, d$q_low), d$vacc,
                                  seed = seed))
    }))
  }))
}

## Simulations are slow, so each analysis caches its raw output as
## results/<name>_n<N_SIMS>.csv. Delete the file to re-run.
run_or_load <- function(name, expr) {
  dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
  file <- file.path(RESULTS_DIR, sprintf("%s_n%d.csv", name, N_SIMS))
  if (file.exists(file)) {
    message("Loading cached simulations: ", file)
    dt <- fread(file)
  } else {
    dt <- expr                 # simulations run here (lazy evaluation)
    fwrite(dt, file)
    message("Saved simulations: ", file)
  }
  dt[, school := factor(school, levels = SCHOOLS)]
  dt[]
}


## When a .qmd is rendered to GitHub Markdown (.md), keep its figures in
## md_figures/ so they survive the HTML render's clean-up and can be
## committed alongside the .md
if (isTRUE(getOption("knitr.in.progress")) &&
    identical(knitr::pandoc_to(), "commonmark")) {
  doc <- sub("\\..*$", "", basename(knitr::current_input()))
  knitr::opts_chunk$set(fig.path = file.path("md_figures", paste0(doc, "-")))
}


## ---------------------------------------------------------------------
## 5. Small utilities
## ---------------------------------------------------------------------

## Wilson score interval for a binomial proportion x / n
wilson <- function(x, n, z = 1.96) {
  p   <- x / n
  d   <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(lo = ctr - hw, hi = ctr + hw)
}

## Difference of two independent proportions, p_new - p_ref, with
## Newcombe's hybrid score 95% CI (built from two Wilson intervals).
## x = number of runs with an outbreak, n = number of runs.
## Assumes the two scenarios are independent (separate seeds).
newcombe_diff <- function(x_new, n_new, x_ref, n_ref) {
  p1 <- x_new / n_new;  p2 <- x_ref / n_ref
  w1 <- wilson(x_new, n_new);  w2 <- wilson(x_ref, n_ref)
  d  <- p1 - p2
  list(diff = d,
       lo   = d - sqrt((p1 - w1$lo)^2 + (w2$hi - p2)^2),
       hi   = d + sqrt((w1$hi - p1)^2 + (p2 - w2$lo)^2))
}

pct      <- scales::percent
vacc_lab <- function(v) paste0(round(100 * v), "%")

SCHOOL_COLORS <- c(Elementary = "#E69F00", Middle = "#009E73", High = "#0072B2")

theme_measles <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.text       = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey95", colour = NA),
          plot.title       = element_text(face = "bold", size = base_size + 2),
          plot.subtitle    = element_text(colour = "grey35"),
          legend.position  = "bottom")
}

## Save a PNG copy to figures/ and return the plot so it also shows inline
save_fig <- function(p, file, width, height, dpi = 200) {
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(FIG_DIR, file), p, width = width, height = height,
         dpi = dpi, bg = "white")
  p
}
