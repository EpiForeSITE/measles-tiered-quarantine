# 10b. Outbreak risk and quarantine burden, one table per school and
vaccination

One table per school and vaccination level, one row per strategy (high /
medium / low-risk quarantine days). **Every number is for one outbreak**
(one run of 100 days starting with one infected student), computed in
each run and averaged over the runs, never added up across runs. Uses
the same runs as file 10 (same seeds): if file 10 has been rendered with
the same `N_RUNS`, nothing is simulated again.

<table>
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Column</th>
<th style="text-align: left;">Meaning</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">P(&gt;10), P(&gt;20)</td>
<td style="text-align: left;">Chance the outbreak ends with more than 10
/ 20 cases</td>
</tr>
<tr class="even">
<td style="text-align: left;">Benefit vs 7/0/0</td>
<td style="text-align: left;">P(&gt;10) with 7/0/0 minus P(&gt;10) with
this strategy, percentage points (95% CI, Newcombe)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Students sent home</td>
<td style="text-align: left;">Different students quarantined at least
once</td>
</tr>
<tr class="even">
<td style="text-align: left;">Rounds</td>
<td style="text-align: left;">Days on which a new group is sent
home</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Once / twice / 3+</td>
<td style="text-align: left;">Of the students sent home, how many times
each was sent home</td>
</tr>
<tr class="even">
<td style="text-align: left;">School days per student sent home</td>
<td style="text-align: left;">School days at home per student sent home
(all their stays added)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">School days per unvaccinated student</td>
<td style="text-align: left;">School days missed per unvaccinated
student (only they can be quarantined)</td>
</tr>
<tr class="even">
<td style="text-align: left;">vs 7/0/0: extra school days; cases
prevented</td>
<td style="text-align: left;">Extra school days per unvaccinated
student, and cases prevented per outbreak, compared with 7/0/0</td>
</tr>
</tbody>
</table>

School days = calendar days × 5/7: the 100 simulated days are about 71
school days, so no school-day column can be larger than 71. | Infected
when sent home | Share of stays where the student was already infected |

<div>

> **Note**
>
> Stays of only 1 day are not counted (a “0-day” tier still keeps a
> student home for 1 day). The low-risk tier is every other unvaccinated
> student, so any strategy with low-risk days sends every unvaccinated
> student home.

</div>

``` r
source("measles_model_dur.R")    # new durations (and RASH_DAYS if set); measles_model.R unchanged
```

    measles from: /uufs/chpc.utah.edu/common/home/u1567499/R_libs_dur/measles | Duration model: 1 (new durations) | Rash contagious days: 0 | R0: 12 | results in: results_dur 
    Students per school: Elementary 339, Middle 466, High 473 
       Elementary : K=57 1=66 2=48 3=63 4=67 5=38 
       Middle : 5=144 6=157 7=100 8=65 
       High : 9=125 10=149 11=108 12=91 

``` r
library(data.table)
library(ggplot2)
if (!exists("vacc_lab")) vacc_lab <- function(v) paste0(round(100 * v), "%")
if (!exists("N_THREADS")) N_THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

## ---- Settings ----------------------------------------------------------------
N_RUNS   <- as.integer(Sys.getenv("N_RUNS", unset = 1000))   # runs per school x vaccination x policy
SEED_BUR <- SEED + 77L                                       # same seeds as file 10 (shares its cache)
VACCS    <- c(0.60, 0.70, 0.80, 0.90)             # same as file 10
POLS     <- c("7/0/0", "21/0/0", "21/14/0", "21/14/7", "21/7/7", "21/21/21")   # same as file 10
REF_POL  <- "7/0/0"           # reference for benefit and cases prevented (minimal quarantine)
K        <- 10L               # outbreak = more than K cases
QST <- c("Quarantined Susceptible", "Quarantined Latent", "Quarantined Prodromal", "Quarantined Recovered")
INF <- c("Quarantined Latent", "Quarantined Prodromal")

COLS <- setNames(c("#C8531E", "#E8A33D", "#9ecae1", "#4292c6", "#7A6FD0", "#08306b"), POLS)
VL   <- vacc_lab(VACCS)
th <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title.position = "plot", legend.position = "bottom",
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))
save10 <- function(p, file, w, h) {
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(FIG_DIR, file), p, width = w, height = h, dpi = 200, bg = "white"); p
}
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

## The burden-by-tier counters need the patched package with "Tracked" parameters
local({
  so <- file.path(find.package("measles"), "libs", "measles.so")
  if (!length(grepRaw("Tracked ", readBin(so, "raw", file.info(so)$size), fixed = TRUE)))
    stop("This measles package has no burden-by-tier counters; install the newest .hpp: ", so)
})
```

``` r
## ---- Model with a daily state recorder --------------------------------------
rec <- new.env()
build_model <- function(school, vacc, d) {
  sizes <- ENROLLMENT[[school]]
  m <- ModelMeaslesMixingRiskQuarantine(
    n = sum(sizes), prevalence = 1 / sum(sizes),
    contact_matrix = SCHOOL_CMAT[[school]], transmission_rate = TRANS_RATE[[school]],
    prop_vaccinated = vacc, detection_rate_quarantine = 0.0, contact_tracing_days_window = 7,
    quarantine_period_high = d[1], quarantine_period_medium = d[2], quarantine_period_low = d[3])
  for (g in sizes) m <- add_entity(m, entity("Grade", as.integer(g), FALSE))
  f  <- function(m) rec$days[[length(rec$days) + 1L]] <- get_agents_states(m)
  ge <- globalevent_fun(f, "State recorder")
  rec$keep[[length(rec$keep) + 1L]] <- list(f, ge, m)    # keep R from garbage-collecting it
  m <- add_globalevent(m, ge)
  verbose_off(m)
  m
}

## ---- One run -> burden measures ---------------------------------------------
one_run <- function(m, seed, d) {
  rec$days <- list()
  run(m, ndays = N_DAYS, seed = seed)
  S <- do.call(rbind, rec$days); nd <- nrow(S)
  isQ <- matrix(S %in% QST, nd)
  ## drop 1-day stays ("0-day" tiers), keep stays of 2+ days
  for (i in which(colSums(isQ) > 0)) {
    r <- rle(isQ[, i]); r$values[r$values & r$lengths == 1] <- FALSE; isQ[, i] <- inverse.rle(r)
  }
  inf   <- matrix(S %in% INF, nd) & isQ
  entry <- isQ & rbind(FALSE, !isQ[-nd, , drop = FALSE])     # first day of each stay
  stays <- colSums(entry); home <- rowSums(isQ); qd <- which(home > 0)

  ## Burden by tier, exact: counted inside the model when each student is sent
  ## home and on each day at home (patched package), read back as parameters
  tp <- function(w) vapply(c("high", "medium", "low"), function(t) get_param(m, paste("Tracked", w, t)), 0)
  tiers <- data.table(tier = c("High", "Medium", "Low"), students = tp("students"), stays = tp("stays"),
                      days = tp("days"), inf = tp("infected stays"))
  list(row = data.table(
         size        = nrow(get_transmissions(m)),
         rounds      = sum(rowSums(entry) > 0),
         students_q  = sum(stays > 0), stays = sum(stays),
         once = sum(stays == 1), twice = sum(stays == 2), more = sum(stays >= 3),
         days_home   = sum(isQ), days_inf = sum(inf),
         stays_inf   = sum(entry & inf),
         peak        = max(home),
         weeks       = if (length(qd)) (max(qd) - min(qd) + 1) / 7 else 0,
         tracked     = sum(tiers$days)),          # days in tiers with more than 0 days
       home = home, tiers = tiers)
}

## ---- One scenario: all its runs, in parallel, cached --------------------------
burden_scenario <- function(school, vacc, pol, seeds) {
  f <- file.path(RESULTS_DIR, sprintf("10_burden2_%s_%d_%s_n%d_seed%d.rds", school, round(100 * vacc),
                                      gsub("/", "-", pol), length(seeds), SEED_BUR))
  if (file.exists(f)) return(readRDS(f))
  d <- as.integer(strsplit(pol, "/")[[1]])
  nc <- max(1L, min(N_THREADS, length(seeds)))
  ch <- split(seeds, ceiling(seq_along(seeds) / ceiling(length(seeds) / nc)))
  out <- unlist(parallel::mclapply(ch, function(s) {
    m <- build_model(school, vacc, d); lapply(s, function(x) one_run(m, x, d)) }, mc.cores = nc),
    recursive = FALSE)
  tag <- data.table(school = school, vacc = vacc, pol = pol)
  tt  <- rbindlist(lapply(seq_along(out), function(k) out[[k]]$tiers[, run := k]))
  res <- list(
    runs  = cbind(tag, rbindlist(lapply(out, `[[`, "row"))),
    home  = cbind(tag, data.table(day = seq_along(out[[1]]$home),
                                  home = Reduce(`+`, lapply(out, `[[`, "home")) / length(out))),
    tiers = cbind(tag, tt))                      # one row per run and tier
  saveRDS(res, f)
  res
}
```

``` r
scen <- CJ(school = SCHOOLS, vacc = VACCS, pol = POLS, sorted = FALSE)
set.seed(SEED_BUR)              # one set.seed, all seeds at once, no repeats
seeds <- split(sample.int(.Machine$integer.max, nrow(scen) * N_RUNS), rep(seq_len(nrow(scen)), each = N_RUNS))
R <- lapply(seq_len(nrow(scen)), function(k) {
  message(sprintf("%-10s %3.0f%% %s", scen$school[k], 100 * scen$vacc[k], scen$pol[k]))
  burden_scenario(scen$school[k], scen$vacc[k], scen$pol[k], seeds[[k]])
})
RUNS  <- rbindlist(lapply(R, `[[`, "runs"))
HOME  <- rbindlist(lapply(R, `[[`, "home"))
TIERS <- rbindlist(lapply(R, `[[`, "tiers"))
for (d in list(RUNS, HOME, TIERS))
  d[, `:=`(school = factor(school, SCHOOLS), pol = factor(pol, POLS), vl = factor(vacc_lab(vacc), VL))]
```

1000 independent runs per school × vaccination × policy, 72 settings in
total.

``` r
ci_mean <- function(v) { m <- mean(v); h <- 1.96 * sd(v) / sqrt(length(v)); c(m, m - h, m + h) }
fmt_ci  <- function(v, dig = 1) { x <- ci_mean(v); sprintf(paste0("%.", dig, "f (%.", dig, "f%.", dig, "f)"), x[1], x[2], x[3]) }

## Every measure is computed for EACH run (one outbreak) first, then averaged
## over runs. Shares and per-stay measures only exist in runs where someone
## was sent home; they are averaged over those runs (q_runs = how many).
RUNS[, `:=`(
  r_once        = fifelse(students_q > 0, once / students_q, NA_real_),
  r_twice       = fifelse(students_q > 0, twice / students_q, NA_real_),
  r_more        = fifelse(students_q > 0, more / students_q, NA_real_),
  r_per_stay    = fifelse(stays > 0, days_home / stays * 5 / 7, NA_real_),         # school days
  r_per_student = fifelse(students_q > 0, days_home / students_q * 5 / 7, NA_real_),
  r_targeted    = fifelse(stays > 0, stays_inf / stays, NA_real_),
  r_inf_days    = fifelse(days_home > 0, days_inf / days_home, NA_real_))]
mn <- function(v) mean(v, na.rm = TRUE)
## Unvaccinated students (only they can be quarantined) and school days
## (the model runs 7 days a week: calendar days x 5/7; 100 days = 71 school days)
RUNS[, nu := round(vapply(as.character(school), function(x) sum(ENROLLMENT[[x]]), 0) * (1 - vacc))]
RUNS[, `:=`(u_days = days_home * 5 / 7 / nu, u_healthy = (days_home - days_inf) * 5 / 7 / nu,
            u_inf = days_inf * 5 / 7 / nu)]

SUM <- RUNS[, .(
  x10 = sum(size > 10), x20 = sum(size > 20), n = .N,
  risk        = mean(size > K),
  cases       = mean(size),
  rounds      = mean(rounds),
  students_q  = mean(students_q),
  q_runs      = sum(students_q > 0),
  once        = mn(r_once),
  twice       = mn(r_twice),
  more        = mn(r_more),
  days_home   = mean(days_home),
  healthy     = mean(days_home - days_inf),
  infected    = mean(days_inf),
  inf_share   = mn(r_inf_days),
  per_stay    = mn(r_per_stay),
  per_student = mn(r_per_student),
  peak        = mean(peak),
  weeks       = mean(weeks),
  zero_tier   = mean(pmax(days_home - tracked, 0)),   # back-to-back 1-day stays of 0-day tiers
  targeted    = mn(r_targeted),
  ## text with 95% CI over runs, for the tables
  t_rounds = fmt_ci(rounds), t_students = fmt_ci(students_q, 0), t_days = fmt_ci(days_home, 0),
  t_peak = fmt_ci(peak, 0), t_weeks = fmt_ci(weeks),
  t_per_stay = fmt_ci(na.omit(r_per_stay)), t_per_student = fmt_ci(na.omit(r_per_student)),
  u_days_m = mean(u_days), u_healthy = mean(u_healthy), u_inf = mean(u_inf), t_u_days = fmt_ci(u_days)),
  by = .(school, vacc, vl, pol)]

## Easier scales for schools: school days (the model runs 7 days a week, so
## calendar days x 5/7), per unvaccinated student and per enrolled student
SUM[, `:=`(n_enrolled = vapply(as.character(school), function(x) sum(ENROLLMENT[[x]]), 0))]
SUM[, `:=`(n_unvacc     = round(n_enrolled * (1 - vacc)),
           school_days  = days_home * 5 / 7)]
SUM[, `:=`(per_unvacc   = school_days / n_unvacc,
           per_enrolled = school_days / n_enrolled)]

## Benefit vs the reference: P(>10) reference - P(>10) strategy, Newcombe 95% CI
if (!exists("wilson")) wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n; c <- (p + z^2 / (2 * n)) / d
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  list(lo = c - h, hi = c + h)
}
if (!exists("newcombe_diff")) newcombe_diff <- function(x_new, n_new, x_ref, n_ref) {
  p1 <- x_new / n_new; p2 <- x_ref / n_ref
  w1 <- wilson(x_new, n_new); w2 <- wilson(x_ref, n_ref); d <- p1 - p2
  list(diff = d, lo = d - sqrt((p1 - w1$lo)^2 + (w2$hi - p2)^2),
       hi = d + sqrt((w1$hi - p1)^2 + (p2 - w2$lo)^2))
}
rx <- SUM[pol == REF_POL, .(school, vacc, xr = x10, nr = n)]
SUM <- rx[SUM, on = .(school, vacc)]
bd  <- newcombe_diff(SUM$xr, SUM$nr, SUM$x10, SUM$n)      # ref - strategy = benefit
SUM[, `:=`(benefit = bd$diff, ben_lo = bd$lo, ben_hi = bd$hi)]

## Is it worth it? Extra days at home per case prevented, vs the reference policy
ref <- SUM[pol == REF_POL, .(school, vacc, ref_cases = cases, ref_days = days_home, ref_u = u_days_m)]
SUM <- ref[SUM, on = .(school, vacc)]
SUM[, prevented := ref_cases - cases]
SUM[, extra_u := u_days_m - ref_u]            # extra school days per unvaccinated student vs reference
SUM[, per_case := fifelse(pol == REF_POL, NA_real_,
                          fifelse(prevented > 0.1, (days_home - ref_days) / prevented, NA_real_))]
fwrite(SUM, file.path(RESULTS_DIR, sprintf("10b_summary_n%d.csv", N_RUNS)))

## Per tier, per outbreak (computed in each run, averaged over runs)
TIERS[, tier := factor(tier, c("High", "Medium", "Low"))]
TIERS[, nu := round(vapply(as.character(school), function(x) sum(ENROLLMENT[[x]]), 0) * (1 - vacc))]
TSUM <- TIERS[, .(students = mean(students), stays = mean(stays), days = mean(days),
                  days_u = mean(days * 5 / 7 / nu),
                  t_students = fmt_ci(students, 0), t_days_u = fmt_ci(days * 5 / 7 / nu),
                  per_student = mean(fifelse(students > 0, days / students * 5 / 7, NA_real_), na.rm = TRUE),
                  inf_share   = mean(fifelse(stays > 0, inf / stays, NA_real_), na.rm = TRUE)),
              by = .(school, vacc, vl, pol, tier)]
TSUM[is.nan(per_student), per_student := NA][is.nan(inf_share), inf_share := NA]
fwrite(TSUM, file.path(RESULTS_DIR, sprintf("10b_by_tier_n%d.csv", N_RUNS)))
```

## Tables for the slides

One table per school and vaccination level.

``` r
pc  <- function(x) sprintf("%.0f%%", 100 * x)
pts <- function(x) sprintf("%+.0f", 100 * x)
BT <- SUM[, .(
  School = school, Vaccinated = vl, Strategy = factor(pol, POLS),
  `P(>10)` = pc(x10 / n), `P(>20)` = pc(x20 / n),
  `Benefit vs 7/0/0 (pts)` = fifelse(pol == REF_POL, "reference",
                                     sprintf("%s (%s to %s)", pts(benefit), pts(ben_lo), pts(ben_hi))),
  `Students sent home` = sprintf("%.0f", students_q),
  Rounds = sprintf("%.1f", rounds),
  `Once / twice / 3+` = sprintf("%s / %s / %s", pc(once), pc(twice), pc(more)),
  `School days per student sent home` = sprintf("%.0f", per_student),
  `School days per unvaccinated student` = sprintf("%.1f", u_days_m),
  `vs 7/0/0: extra school days; cases prevented` = fifelse(pol == REF_POL, "reference",
      sprintf("%+.1f; %s", extra_u, fifelse(prevented > 0.1, sprintf("%.1f", prevented), "none"))),
  `Infected when sent home` = pc(targeted))]
setorder(BT, School, Vaccinated, Strategy)
fwrite(BT, file.path(RESULTS_DIR, sprintf("10b_slide_table_n%d.csv", N_RUNS)))
```

<div class="panel-tabset">

### Elementary

**Elementary school, 60% vaccinated** (339 students, 136 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">97%</td>
<td style="text-align: right;">95%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">130</td>
<td style="text-align: right;">9.1</td>
<td style="text-align: right;">9% / 14% / 77%</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">16.8</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">13%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">67%</td>
<td style="text-align: right;">50%</td>
<td style="text-align: right;">+30 (+27 to +33)</td>
<td style="text-align: right;">115</td>
<td style="text-align: right;">4.5</td>
<td style="text-align: right;">44% / 22% / 34%</td>
<td style="text-align: right;">14</td>
<td style="text-align: right;">11.7</td>
<td style="text-align: right;">-5.1; 37.1</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">48%</td>
<td style="text-align: right;">27%</td>
<td style="text-align: right;">+49 (+46 to +52)</td>
<td style="text-align: right;">111</td>
<td style="text-align: right;">4.0</td>
<td style="text-align: right;">59% / 21% / 20%</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">10.9</td>
<td style="text-align: right;">-5.9; 44.1</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">38%</td>
<td style="text-align: right;">19%</td>
<td style="text-align: right;">+59 (+56 to +62)</td>
<td style="text-align: right;">133</td>
<td style="text-align: right;">3.2</td>
<td style="text-align: right;">23% / 42% / 35%</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">17.5</td>
<td style="text-align: right;">+0.7; 46.5</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">61%</td>
<td style="text-align: right;">37%</td>
<td style="text-align: right;">+36 (+33 to +39)</td>
<td style="text-align: right;">134</td>
<td style="text-align: right;">3.6</td>
<td style="text-align: right;">16% / 36% / 48%</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">20.8</td>
<td style="text-align: right;">+4.0; 40.8</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">31%</td>
<td style="text-align: right;">11%</td>
<td style="text-align: right;">+66 (+62 to +69)</td>
<td style="text-align: right;">133</td>
<td style="text-align: right;">1.2</td>
<td style="text-align: right;">85% / 14% / 1%</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">18.0</td>
<td style="text-align: right;">+1.2; 48.8</td>
<td style="text-align: right;">5%</td>
</tr>
</tbody>
</table>

**Elementary school, 70% vaccinated** (339 students, 102 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">91%</td>
<td style="text-align: right;">86%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">92</td>
<td style="text-align: right;">8.3</td>
<td style="text-align: right;">16% / 16% / 68%</td>
<td style="text-align: right;">14</td>
<td style="text-align: right;">12.6</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">13%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">53%</td>
<td style="text-align: right;">30%</td>
<td style="text-align: right;">+38 (+34 to +41)</td>
<td style="text-align: right;">78</td>
<td style="text-align: right;">3.7</td>
<td style="text-align: right;">58% / 21% / 22%</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">8.9</td>
<td style="text-align: right;">-3.7; 25.9</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">34%</td>
<td style="text-align: right;">14%</td>
<td style="text-align: right;">+56 (+53 to +60)</td>
<td style="text-align: right;">73</td>
<td style="text-align: right;">3.3</td>
<td style="text-align: right;">70% / 19% / 11%</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">8.4</td>
<td style="text-align: right;">-4.2; 30.8</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">27%</td>
<td style="text-align: right;">9%</td>
<td style="text-align: right;">+63 (+60 to +67)</td>
<td style="text-align: right;">99</td>
<td style="text-align: right;">3.0</td>
<td style="text-align: right;">25% / 46% / 29%</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">15.8</td>
<td style="text-align: right;">+3.2; 32.1</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">44%</td>
<td style="text-align: right;">22%</td>
<td style="text-align: right;">+46 (+43 to +50)</td>
<td style="text-align: right;">99</td>
<td style="text-align: right;">3.3</td>
<td style="text-align: right;">18% / 37% / 44%</td>
<td style="text-align: right;">19</td>
<td style="text-align: right;">18.3</td>
<td style="text-align: right;">+5.7; 28.3</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">22%</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">+69 (+65 to +72)</td>
<td style="text-align: right;">99</td>
<td style="text-align: right;">1.1</td>
<td style="text-align: right;">87% / 12% / 1%</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">17.3</td>
<td style="text-align: right;">+4.7; 33.2</td>
<td style="text-align: right;">5%</td>
</tr>
</tbody>
</table>

**Elementary school, 80% vaccinated** (339 students, 68 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">74%</td>
<td style="text-align: right;">57%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">51</td>
<td style="text-align: right;">6.2</td>
<td style="text-align: right;">31% / 24% / 45%</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">7.3</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">12%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">31%</td>
<td style="text-align: right;">11%</td>
<td style="text-align: right;">+43 (+39 to +47)</td>
<td style="text-align: right;">38</td>
<td style="text-align: right;">2.6</td>
<td style="text-align: right;">77% / 16% / 8%</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">6.0</td>
<td style="text-align: right;">-1.3; 13.2</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">17%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">+57 (+53 to +61)</td>
<td style="text-align: right;">40</td>
<td style="text-align: right;">2.5</td>
<td style="text-align: right;">82% / 14% / 4%</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">6.4</td>
<td style="text-align: right;">-0.9; 15.6</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">14%</td>
<td style="text-align: right;">2%</td>
<td style="text-align: right;">+60 (+57 to +64)</td>
<td style="text-align: right;">64</td>
<td style="text-align: right;">2.6</td>
<td style="text-align: right;">27% / 48% / 25%</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">13.9</td>
<td style="text-align: right;">+6.5; 15.9</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">27%</td>
<td style="text-align: right;">7%</td>
<td style="text-align: right;">+48 (+44 to +51)</td>
<td style="text-align: right;">64</td>
<td style="text-align: right;">2.8</td>
<td style="text-align: right;">22% / 40% / 38%</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">15.4</td>
<td style="text-align: right;">+8.0; 14.0</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">9%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+66 (+62 to +69)</td>
<td style="text-align: right;">65</td>
<td style="text-align: right;">1.1</td>
<td style="text-align: right;">89% / 11% / 0%</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">16.2</td>
<td style="text-align: right;">+8.8; 16.7</td>
<td style="text-align: right;">4%</td>
</tr>
</tbody>
</table>

**Elementary school, 90% vaccinated** (339 students, 34 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">23%</td>
<td style="text-align: right;">2%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">2.7</td>
<td style="text-align: right;">65% / 24% / 11%</td>
<td style="text-align: right;">7</td>
<td style="text-align: right;">2.4</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">12%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">7%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+15 (+12 to +18)</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">1.6</td>
<td style="text-align: right;">92% / 7% / 1%</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">3.7</td>
<td style="text-align: right;">+1.3; 2.3</td>
<td style="text-align: right;">12%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+19 (+16 to +22)</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">1.6</td>
<td style="text-align: right;">93% / 6% / 1%</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">4.4</td>
<td style="text-align: right;">+1.9; 2.9</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+18 (+16 to +21)</td>
<td style="text-align: right;">30</td>
<td style="text-align: right;">1.9</td>
<td style="text-align: right;">43% / 42% / 15%</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">10.7</td>
<td style="text-align: right;">+8.2; 2.7</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+18 (+15 to +21)</td>
<td style="text-align: right;">30</td>
<td style="text-align: right;">1.9</td>
<td style="text-align: right;">42% / 39% / 19%</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">10.4</td>
<td style="text-align: right;">+8.0; 2.7</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">2%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+21 (+18 to +24)</td>
<td style="text-align: right;">30</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">92% / 7% / 0%</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">14.5</td>
<td style="text-align: right;">+12.1; 3.2</td>
<td style="text-align: right;">4%</td>
</tr>
</tbody>
</table>

### Middle

**Middle school, 60% vaccinated** (466 students, 186 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">77%</td>
<td style="text-align: right;">52%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">155</td>
<td style="text-align: right;">5.4</td>
<td style="text-align: right;">22% / 21% / 57%</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">9.2</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">16%</td>
<td style="text-align: right;">2%</td>
<td style="text-align: right;">+60 (+57 to +64)</td>
<td style="text-align: right;">118</td>
<td style="text-align: right;">1.7</td>
<td style="text-align: right;">82% / 13% / 5%</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">5.7</td>
<td style="text-align: right;">-3.5; 14.9</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">16%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">+61 (+57 to +64)</td>
<td style="text-align: right;">121</td>
<td style="text-align: right;">1.8</td>
<td style="text-align: right;">82% / 13% / 5%</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">5.7</td>
<td style="text-align: right;">-3.5; 14.6</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">15%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">+61 (+58 to +65)</td>
<td style="text-align: right;">180</td>
<td style="text-align: right;">2.2</td>
<td style="text-align: right;">35% / 49% / 17%</td>
<td style="text-align: right;">14</td>
<td style="text-align: right;">13.2</td>
<td style="text-align: right;">+4.1; 14.6</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">16%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">+61 (+57 to +64)</td>
<td style="text-align: right;">182</td>
<td style="text-align: right;">2.2</td>
<td style="text-align: right;">31% / 51% / 18%</td>
<td style="text-align: right;">14</td>
<td style="text-align: right;">13.5</td>
<td style="text-align: right;">+4.4; 14.7</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">15%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">+62 (+58 to +65)</td>
<td style="text-align: right;">180</td>
<td style="text-align: right;">1.1</td>
<td style="text-align: right;">91% / 8% / 1%</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">16.6</td>
<td style="text-align: right;">+7.5; 15.0</td>
<td style="text-align: right;">2%</td>
</tr>
</tbody>
</table>

**Middle school, 70% vaccinated** (466 students, 140 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">57%</td>
<td style="text-align: right;">26%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">97</td>
<td style="text-align: right;">4.2</td>
<td style="text-align: right;">34% / 23% / 42%</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">6.3</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">9%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">+48 (+44 to +52)</td>
<td style="text-align: right;">69</td>
<td style="text-align: right;">1.5</td>
<td style="text-align: right;">87% / 10% / 3%</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">4.9</td>
<td style="text-align: right;">-1.5; 8.6</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">9%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">+48 (+44 to +52)</td>
<td style="text-align: right;">73</td>
<td style="text-align: right;">1.5</td>
<td style="text-align: right;">89% / 9% / 3%</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">4.9</td>
<td style="text-align: right;">-1.4; 8.7</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">7%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+50 (+47 to +53)</td>
<td style="text-align: right;">132</td>
<td style="text-align: right;">2.1</td>
<td style="text-align: right;">36% / 53% / 12%</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">11.9</td>
<td style="text-align: right;">+5.6; 8.8</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">10%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">+48 (+44 to +51)</td>
<td style="text-align: right;">132</td>
<td style="text-align: right;">2.0</td>
<td style="text-align: right;">36% / 51% / 13%</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">12.2</td>
<td style="text-align: right;">+5.8; 8.6</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">8%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">+49 (+45 to +52)</td>
<td style="text-align: right;">132</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">93% / 7% / 0%</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">15.9</td>
<td style="text-align: right;">+9.6; 8.7</td>
<td style="text-align: right;">2%</td>
</tr>
</tbody>
</table>

**Middle school, 80% vaccinated** (466 students, 93 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">32%</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">46</td>
<td style="text-align: right;">2.9</td>
<td style="text-align: right;">47% / 24% / 29%</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">4.0</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+28 (+25 to +31)</td>
<td style="text-align: right;">37</td>
<td style="text-align: right;">1.2</td>
<td style="text-align: right;">93% / 6% / 1%</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">4.3</td>
<td style="text-align: right;">+0.2; 4.2</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+28 (+25 to +31)</td>
<td style="text-align: right;">37</td>
<td style="text-align: right;">1.2</td>
<td style="text-align: right;">94% / 5% / 1%</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">4.2</td>
<td style="text-align: right;">+0.2; 4.3</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+28 (+25 to +31)</td>
<td style="text-align: right;">84</td>
<td style="text-align: right;">1.8</td>
<td style="text-align: right;">44% / 47% / 9%</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">10.6</td>
<td style="text-align: right;">+6.6; 4.2</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+29 (+26 to +32)</td>
<td style="text-align: right;">84</td>
<td style="text-align: right;">1.7</td>
<td style="text-align: right;">44% / 47% / 9%</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">10.5</td>
<td style="text-align: right;">+6.5; 4.3</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+28 (+25 to +31)</td>
<td style="text-align: right;">84</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">93% / 7% / 1%</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">14.9</td>
<td style="text-align: right;">+10.9; 4.1</td>
<td style="text-align: right;">2%</td>
</tr>
</tbody>
</table>

**Middle school, 90% vaccinated** (466 students, 47 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">1.6</td>
<td style="text-align: right;">62% / 23% / 14%</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">2.0</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+4 (+3 to +5)</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">0.9</td>
<td style="text-align: right;">96% / 4% / 0%</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">3.5</td>
<td style="text-align: right;">+1.4; 0.9</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+4 (+2 to +5)</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">96% / 4% / 0%</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">3.7</td>
<td style="text-align: right;">+1.7; 0.9</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+4 (+3 to +6)</td>
<td style="text-align: right;">39</td>
<td style="text-align: right;">1.4</td>
<td style="text-align: right;">59% / 37% / 4%</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">8.5</td>
<td style="text-align: right;">+6.5; 0.9</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+4 (+3 to +5)</td>
<td style="text-align: right;">38</td>
<td style="text-align: right;">1.4</td>
<td style="text-align: right;">59% / 37% / 5%</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">8.5</td>
<td style="text-align: right;">+6.5; 0.9</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+4 (+3 to +5)</td>
<td style="text-align: right;">40</td>
<td style="text-align: right;">0.9</td>
<td style="text-align: right;">96% / 4% / 0%</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">13.3</td>
<td style="text-align: right;">+11.3; 0.9</td>
<td style="text-align: right;">2%</td>
</tr>
</tbody>
</table>

### High

**High school, 60% vaccinated** (473 students, 189 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">98%</td>
<td style="text-align: right;">97%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">183</td>
<td style="text-align: right;">8.7</td>
<td style="text-align: right;">8% / 14% / 77%</td>
<td style="text-align: right;">23</td>
<td style="text-align: right;">22.2</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">12%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">84%</td>
<td style="text-align: right;">66%</td>
<td style="text-align: right;">+14 (+11 to +16)</td>
<td style="text-align: right;">177</td>
<td style="text-align: right;">5.1</td>
<td style="text-align: right;">36% / 28% / 35%</td>
<td style="text-align: right;">20</td>
<td style="text-align: right;">18.7</td>
<td style="text-align: right;">-3.5; 40.0</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">57%</td>
<td style="text-align: right;">37%</td>
<td style="text-align: right;">+40 (+37 to +44)</td>
<td style="text-align: right;">174</td>
<td style="text-align: right;">4.7</td>
<td style="text-align: right;">58% / 24% / 18%</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">17.2</td>
<td style="text-align: right;">-5.0; 50.6</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">47%</td>
<td style="text-align: right;">29%</td>
<td style="text-align: right;">+51 (+48 to +54)</td>
<td style="text-align: right;">187</td>
<td style="text-align: right;">3.5</td>
<td style="text-align: right;">28% / 39% / 33%</td>
<td style="text-align: right;">21</td>
<td style="text-align: right;">20.6</td>
<td style="text-align: right;">-1.6; 53.4</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">73%</td>
<td style="text-align: right;">53%</td>
<td style="text-align: right;">+24 (+21 to +27)</td>
<td style="text-align: right;">187</td>
<td style="text-align: right;">4.1</td>
<td style="text-align: right;">16% / 29% / 55%</td>
<td style="text-align: right;">26</td>
<td style="text-align: right;">25.6</td>
<td style="text-align: right;">+3.4; 45.0</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">34%</td>
<td style="text-align: right;">17%</td>
<td style="text-align: right;">+64 (+60 to +67)</td>
<td style="text-align: right;">188</td>
<td style="text-align: right;">1.2</td>
<td style="text-align: right;">82% / 16% / 1%</td>
<td style="text-align: right;">20</td>
<td style="text-align: right;">19.4</td>
<td style="text-align: right;">-2.8; 57.0</td>
<td style="text-align: right;">4%</td>
</tr>
</tbody>
</table>

**High school, 70% vaccinated** (473 students, 142 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">95%</td>
<td style="text-align: right;">94%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">133</td>
<td style="text-align: right;">8.6</td>
<td style="text-align: right;">10% / 13% / 77%</td>
<td style="text-align: right;">20</td>
<td style="text-align: right;">19.0</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">11%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">73%</td>
<td style="text-align: right;">50%</td>
<td style="text-align: right;">+22 (+19 to +25)</td>
<td style="text-align: right;">122</td>
<td style="text-align: right;">4.6</td>
<td style="text-align: right;">45% / 28% / 28%</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">16.1</td>
<td style="text-align: right;">-2.9; 30.2</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">42%</td>
<td style="text-align: right;">19%</td>
<td style="text-align: right;">+53 (+50 to +56)</td>
<td style="text-align: right;">119</td>
<td style="text-align: right;">3.9</td>
<td style="text-align: right;">66% / 23% / 11%</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">14.4</td>
<td style="text-align: right;">-4.6; 38.6</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">36%</td>
<td style="text-align: right;">17%</td>
<td style="text-align: right;">+59 (+56 to +62)</td>
<td style="text-align: right;">139</td>
<td style="text-align: right;">3.4</td>
<td style="text-align: right;">26% / 40% / 34%</td>
<td style="text-align: right;">20</td>
<td style="text-align: right;">19.4</td>
<td style="text-align: right;">+0.4; 39.6</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">60%</td>
<td style="text-align: right;">34%</td>
<td style="text-align: right;">+34 (+31 to +38)</td>
<td style="text-align: right;">138</td>
<td style="text-align: right;">3.8</td>
<td style="text-align: right;">17% / 29% / 54%</td>
<td style="text-align: right;">24</td>
<td style="text-align: right;">22.9</td>
<td style="text-align: right;">+3.8; 34.2</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">26%</td>
<td style="text-align: right;">8%</td>
<td style="text-align: right;">+68 (+65 to +71)</td>
<td style="text-align: right;">139</td>
<td style="text-align: right;">1.2</td>
<td style="text-align: right;">84% / 15% / 1%</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">18.0</td>
<td style="text-align: right;">-1.1; 42.0</td>
<td style="text-align: right;">4%</td>
</tr>
</tbody>
</table>

**High school, 80% vaccinated** (473 students, 95 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">82%</td>
<td style="text-align: right;">75%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">79</td>
<td style="text-align: right;">7.2</td>
<td style="text-align: right;">20% / 16% / 64%</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">12.8</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">48%</td>
<td style="text-align: right;">18%</td>
<td style="text-align: right;">+34 (+30 to +38)</td>
<td style="text-align: right;">70</td>
<td style="text-align: right;">3.3</td>
<td style="text-align: right;">65% / 23% / 12%</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">11.7</td>
<td style="text-align: right;">-1.1; 17.6</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">28%</td>
<td style="text-align: right;">8%</td>
<td style="text-align: right;">+54 (+50 to +58)</td>
<td style="text-align: right;">71</td>
<td style="text-align: right;">3.2</td>
<td style="text-align: right;">74% / 19% / 6%</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">11.8</td>
<td style="text-align: right;">-1.1; 21.1</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">23%</td>
<td style="text-align: right;">7%</td>
<td style="text-align: right;">+60 (+56 to +63)</td>
<td style="text-align: right;">91</td>
<td style="text-align: right;">2.9</td>
<td style="text-align: right;">32% / 40% / 28%</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">16.7</td>
<td style="text-align: right;">+3.8; 21.8</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">39%</td>
<td style="text-align: right;">14%</td>
<td style="text-align: right;">+43 (+39 to +47)</td>
<td style="text-align: right;">91</td>
<td style="text-align: right;">3.3</td>
<td style="text-align: right;">20% / 33% / 47%</td>
<td style="text-align: right;">20</td>
<td style="text-align: right;">19.5</td>
<td style="text-align: right;">+6.7; 18.8</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">13%</td>
<td style="text-align: right;">2%</td>
<td style="text-align: right;">+69 (+65 to +72)</td>
<td style="text-align: right;">89</td>
<td style="text-align: right;">1.1</td>
<td style="text-align: right;">87% / 12% / 1%</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">16.4</td>
<td style="text-align: right;">+3.5; 23.5</td>
<td style="text-align: right;">3%</td>
</tr>
</tbody>
</table>

**High school, 90% vaccinated** (473 students, 47 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 3%" />
<col style="width: 10%" />
<col style="width: 8%" />
<col style="width: 3%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">P(&gt;10)</th>
<th style="text-align: right;">P(&gt;20)</th>
<th style="text-align: right;">Benefit vs 7/0/0 (pts)</th>
<th style="text-align: right;">Students sent home</th>
<th style="text-align: right;">Rounds</th>
<th style="text-align: right;">Once / twice / 3+</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">School days per unvaccinated student</th>
<th style="text-align: right;">vs 7/0/0: extra school days; cases
prevented</th>
<th style="text-align: right;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">7/0/0</td>
<td style="text-align: right;">40%</td>
<td style="text-align: right;">12%</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">27</td>
<td style="text-align: right;">3.7</td>
<td style="text-align: right;">52% / 25% / 23%</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">5.4</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">12%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+27 (+23 to +31)</td>
<td style="text-align: right;">23</td>
<td style="text-align: right;">2.0</td>
<td style="text-align: right;">87% / 11% / 2%</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">7.3</td>
<td style="text-align: right;">+1.9; 4.2</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/14/0</td>
<td style="text-align: right;">8%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+32 (+28 to +35)</td>
<td style="text-align: right;">27</td>
<td style="text-align: right;">2.0</td>
<td style="text-align: right;">88% / 10% / 2%</td>
<td style="text-align: right;">14</td>
<td style="text-align: right;">8.2</td>
<td style="text-align: right;">+2.8; 5.1</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/7</td>
<td style="text-align: right;">6%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+34 (+30 to +37)</td>
<td style="text-align: right;">42</td>
<td style="text-align: right;">2.1</td>
<td style="text-align: right;">44% / 40% / 16%</td>
<td style="text-align: right;">14</td>
<td style="text-align: right;">12.8</td>
<td style="text-align: right;">+7.4; 5.3</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">11%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">+28 (+25 to +32)</td>
<td style="text-align: right;">42</td>
<td style="text-align: right;">2.2</td>
<td style="text-align: right;">39% / 34% / 27%</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">13.3</td>
<td style="text-align: right;">+7.9; 4.4</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">+36 (+33 to +39)</td>
<td style="text-align: right;">42</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">89% / 10% / 1%</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">15.3</td>
<td style="text-align: right;">+10.0; 5.7</td>
<td style="text-align: right;">3%</td>
</tr>
</tbody>
</table>

</div>

## Days at home by risk tier

Exact counts from the model, per outbreak (averaged over runs).

``` r
TT <- TSUM[students >= 0.05, .(School = school, Vaccinated = vl, Strategy = factor(pol, POLS), Tier = tier,
            Students = sprintf("%.1f", students), `School days per student` = sprintf("%.0f", per_student),
            `School days per unvaccinated student` = sprintf("%.1f", days_u), `Infected when sent home` = pc(inf_share))]
setorder(TT, School, Vaccinated, Strategy, Tier)
fwrite(TT, file.path(RESULTS_DIR, sprintf("10b_by_tier_table_n%d.csv", N_RUNS)))
knitr::kable(TT)
```

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 8%" />
<col style="width: 6%" />
<col style="width: 5%" />
<col style="width: 6%" />
<col style="width: 18%" />
<col style="width: 28%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">School</th>
<th style="text-align: left;">Vaccinated</th>
<th style="text-align: left;">Strategy</th>
<th style="text-align: left;">Tier</th>
<th style="text-align: left;">Students</th>
<th style="text-align: left;">School days per student</th>
<th style="text-align: left;">School days per unvaccinated student</th>
<th style="text-align: left;">Infected when sent home</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">105.3</td>
<td style="text-align: left;">13</td>
<td style="text-align: left;">10.0</td>
<td style="text-align: left;">23%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">72.6</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">8.6</td>
<td style="text-align: left;">16%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">61.6</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">7.2</td>
<td style="text-align: left;">12%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">22.4</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.7</td>
<td style="text-align: left;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">41.8</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">4.8</td>
<td style="text-align: left;">14%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">15.6</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.2</td>
<td style="text-align: left;">12%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">112.4</td>
<td style="text-align: left;">14</td>
<td style="text-align: left;">11.5</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">49.7</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">5.6</td>
<td style="text-align: left;">17%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">19.2</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">0.7</td>
<td style="text-align: left;">11%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">120.0</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">14.4</td>
<td style="text-align: left;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">25.5</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">2.9</td>
<td style="text-align: left;">17%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">13.7</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">1.5</td>
<td style="text-align: left;">11%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">102.2</td>
<td style="text-align: left;">18</td>
<td style="text-align: left;">13.6</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">68.1</td>
<td style="text-align: left;">12</td>
<td style="text-align: left;">7.9</td>
<td style="text-align: left;">22%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">45.5</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">6.9</td>
<td style="text-align: left;">17%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">37.3</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">5.7</td>
<td style="text-align: left;">13%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">14.8</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.5</td>
<td style="text-align: left;">10%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">28.4</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">4.3</td>
<td style="text-align: left;">14%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">11.5</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.1</td>
<td style="text-align: left;">10%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">82.6</td>
<td style="text-align: left;">13</td>
<td style="text-align: left;">10.4</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">32.7</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">4.8</td>
<td style="text-align: left;">17%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">14.3</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">0.7</td>
<td style="text-align: left;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">87.8</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">12.7</td>
<td style="text-align: left;">2%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">18.6</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">2.8</td>
<td style="text-align: left;">17%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">10.0</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">1.4</td>
<td style="text-align: left;">11%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">75.5</td>
<td style="text-align: left;">17</td>
<td style="text-align: left;">13.1</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">33.7</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">5.1</td>
<td style="text-align: left;">20%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">22.6</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">5.1</td>
<td style="text-align: left;">16%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">20.0</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">4.5</td>
<td style="text-align: left;">13%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">8.9</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.3</td>
<td style="text-align: left;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">16.8</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">3.7</td>
<td style="text-align: left;">14%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">7.6</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.1</td>
<td style="text-align: left;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">53.1</td>
<td style="text-align: left;">11</td>
<td style="text-align: left;">9.0</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">18.6</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">4.1</td>
<td style="text-align: left;">16%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">9.2</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">0.7</td>
<td style="text-align: left;">10%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">55.9</td>
<td style="text-align: left;">13</td>
<td style="text-align: left;">10.5</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">11.7</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">2.6</td>
<td style="text-align: left;">15%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">6.5</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">1.4</td>
<td style="text-align: left;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">49.1</td>
<td style="text-align: left;">17</td>
<td style="text-align: left;">12.2</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">8.8</td>
<td style="text-align: left;">8</td>
<td style="text-align: left;">2.1</td>
<td style="text-align: left;">15%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">7.8</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">3.4</td>
<td style="text-align: left;">15%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">7.2</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">3.1</td>
<td style="text-align: left;">14%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">3.6</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.1</td>
<td style="text-align: left;">10%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">6.9</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">3.0</td>
<td style="text-align: left;">13%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">3.6</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.1</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">23.8</td>
<td style="text-align: left;">9</td>
<td style="text-align: left;">6.6</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">6.7</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">2.9</td>
<td style="text-align: left;">14%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">3.8</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">0.6</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">24.7</td>
<td style="text-align: left;">9</td>
<td style="text-align: left;">6.9</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">5.4</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">2.4</td>
<td style="text-align: left;">15%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">2.9</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">1.3</td>
<td style="text-align: left;">11%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">22.7</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">10.9</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">59.5</td>
<td style="text-align: left;">18</td>
<td style="text-align: left;">5.9</td>
<td style="text-align: left;">14%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">53.2</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.7</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">52.1</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.6</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.9</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">50.6</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.5</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.8</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">131.7</td>
<td style="text-align: left;">12</td>
<td style="text-align: left;">8.7</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">52.0</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.6</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.8</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">134.3</td>
<td style="text-align: left;">12</td>
<td style="text-align: left;">8.9</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">48.4</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.3</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.7</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">0.1</td>
<td style="text-align: left;">4%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">131.5</td>
<td style="text-align: left;">17</td>
<td style="text-align: left;">12.3</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">40.5</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.6</td>
<td style="text-align: left;">12%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">36.7</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.2</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">37.4</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.3</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.7</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">36.2</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.1</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.6</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">96.3</td>
<td style="text-align: left;">11</td>
<td style="text-align: left;">7.7</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">37.6</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.3</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.6</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">96.9</td>
<td style="text-align: left;">11</td>
<td style="text-align: left;">7.8</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">36.1</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.2</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.5</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">0.1</td>
<td style="text-align: left;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">95.9</td>
<td style="text-align: left;">17</td>
<td style="text-align: left;">11.7</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">24.2</td>
<td style="text-align: left;">13</td>
<td style="text-align: left;">3.3</td>
<td style="text-align: left;">10%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">23.0</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">3.9</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">23.3</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">3.9</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.4</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">23.2</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">3.9</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.3</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">4%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">61.2</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">6.6</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">22.9</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">3.9</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.4</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">2%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">61.6</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">6.6</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">22.8</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">3.9</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.3</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">0.1</td>
<td style="text-align: left;">4%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">61.4</td>
<td style="text-align: left;">17</td>
<td style="text-align: left;">11.0</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">10.5</td>
<td style="text-align: left;">8</td>
<td style="text-align: left;">1.9</td>
<td style="text-align: left;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">10.4</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">3.4</td>
<td style="text-align: left;">7%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">11.0</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">3.6</td>
<td style="text-align: left;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.2</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">3%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">10.7</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">3.5</td>
<td style="text-align: left;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.2</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">28.3</td>
<td style="text-align: left;">8</td>
<td style="text-align: left;">4.9</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">10.5</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">3.4</td>
<td style="text-align: left;">7%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.2</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">28.1</td>
<td style="text-align: left;">8</td>
<td style="text-align: left;">5.0</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">10.5</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">3.4</td>
<td style="text-align: left;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">0.2</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">0.0</td>
<td style="text-align: left;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">28.9</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">9.8</td>
<td style="text-align: left;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">179.1</td>
<td style="text-align: left;">18</td>
<td style="text-align: left;">17.2</td>
<td style="text-align: left;">12%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">160.9</td>
<td style="text-align: left;">19</td>
<td style="text-align: left;">16.1</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">143.0</td>
<td style="text-align: left;">17</td>
<td style="text-align: left;">13.5</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">43.7</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">2.4</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">104.7</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">9.2</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">31.8</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.7</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">147.4</td>
<td style="text-align: left;">12</td>
<td style="text-align: left;">9.7</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">124.0</td>
<td style="text-align: left;">17</td>
<td style="text-align: left;">11.1</td>
<td style="text-align: left;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">41.8</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">1.2</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">166.6</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">13.2</td>
<td style="text-align: left;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">58.6</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">5.0</td>
<td style="text-align: left;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">27.8</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">2.2</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">126.9</td>
<td style="text-align: left;">18</td>
<td style="text-align: left;">12.2</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">129.1</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">14.8</td>
<td style="text-align: left;">12%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">109.8</td>
<td style="text-align: left;">18</td>
<td style="text-align: left;">14.2</td>
<td style="text-align: left;">7%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">94.0</td>
<td style="text-align: left;">17</td>
<td style="text-align: left;">11.3</td>
<td style="text-align: left;">4%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">30.1</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">2.2</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">73.3</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">8.4</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">24.4</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.7</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">109.7</td>
<td style="text-align: left;">12</td>
<td style="text-align: left;">9.2</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">83.6</td>
<td style="text-align: left;">17</td>
<td style="text-align: left;">9.8</td>
<td style="text-align: left;">7%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">30.5</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">1.1</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">120.7</td>
<td style="text-align: left;">14</td>
<td style="text-align: left;">11.8</td>
<td style="text-align: left;">2%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">41.6</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.6</td>
<td style="text-align: left;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">21.1</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">2.2</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">90.8</td>
<td style="text-align: left;">17</td>
<td style="text-align: left;">11.1</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">73.3</td>
<td style="text-align: left;">13</td>
<td style="text-align: left;">10.4</td>
<td style="text-align: left;">10%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">60.0</td>
<td style="text-align: left;">17</td>
<td style="text-align: left;">10.8</td>
<td style="text-align: left;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">53.0</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">9.2</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">18.7</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">2.0</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">43.5</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">7.3</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">15.6</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.7</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">68.8</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">7.7</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">48.9</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">8.4</td>
<td style="text-align: left;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">19.9</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">1.1</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">76.7</td>
<td style="text-align: left;">12</td>
<td style="text-align: left;">9.9</td>
<td style="text-align: left;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">25.7</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.2</td>
<td style="text-align: left;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">13.0</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">2.0</td>
<td style="text-align: left;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">57.4</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">10.1</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">24.1</td>
<td style="text-align: left;">9</td>
<td style="text-align: left;">4.9</td>
<td style="text-align: left;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">20.4</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">7.0</td>
<td style="text-align: left;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">19.0</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">6.3</td>
<td style="text-align: left;">4%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">8.0</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.8</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">16.8</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">5.5</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">7.0</td>
<td style="text-align: left;">10</td>
<td style="text-align: left;">1.5</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">29.7</td>
<td style="text-align: left;">9</td>
<td style="text-align: left;">5.7</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">17.5</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">5.8</td>
<td style="text-align: left;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">8.5</td>
<td style="text-align: left;">5</td>
<td style="text-align: left;">1.0</td>
<td style="text-align: left;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">32.1</td>
<td style="text-align: left;">9</td>
<td style="text-align: left;">6.4</td>
<td style="text-align: left;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: left;">12.0</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">4.0</td>
<td style="text-align: left;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: left;">6.3</td>
<td style="text-align: left;">15</td>
<td style="text-align: left;">2.0</td>
<td style="text-align: left;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: left;">26.9</td>
<td style="text-align: left;">16</td>
<td style="text-align: left;">9.4</td>
<td style="text-align: left;">1%</td>
</tr>
</tbody>
</table>
