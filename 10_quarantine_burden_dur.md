# 10. What quarantine means for a school: rounds, students, days at home

Every number here is **per outbreak** (one simulation that starts with
one infected student and runs 100 days), averaged over independent runs.

<table>
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question a school asks</th>
<th style="text-align: left;">Measure</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">How many students will be sent home?</td>
<td style="text-align: left;">Students quarantined at least once</td>
</tr>
<tr class="even">
<td style="text-align: left;">How often will we have to send a group
home?</td>
<td style="text-align: left;">Rounds: days on which a new group is sent
home</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Will the same child be sent home
repeatedly?</td>
<td style="text-align: left;">Students sent home once, twice, 3+
times</td>
</tr>
<tr class="even">
<td style="text-align: left;">How many school days are lost?</td>
<td style="text-align: left;">School days an unvaccinated student
misses, healthy vs infected days</td>
</tr>
<tr class="odd">
<td style="text-align: left;">How long is each child out?</td>
<td style="text-align: left;">School days per stay at home, and per
student sent home</td>
</tr>
<tr class="even">
<td style="text-align: left;">How many children are out at once?</td>
<td style="text-align: left;">Most students at home on one day</td>
</tr>
<tr class="odd">
<td style="text-align: left;">How long does the disruption last?</td>
<td style="text-align: left;">Weeks from the first to the last day
anyone is at home</td>
</tr>
<tr class="even">
<td style="text-align: left;">Is it worth it?</td>
<td style="text-align: left;">Extra school days per unvaccinated
student, and cases prevented, compared with 7/0/0</td>
</tr>
<tr class="odd">
<td style="text-align: left;">How well targeted is it?</td>
<td style="text-align: left;">Share of stays at home where the student
was infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">How many students and days in each
tier?</td>
<td style="text-align: left;">Students sent home and school days each,
per tier (high, medium, low)</td>
</tr>
</tbody>
</table>

<div>

> **Note**
>
> **What counts as a quarantine.** A “0-day” tier still keeps a student
> home for 1 day in the model (release rule: days since sent home ≥ tier
> days). Those 1-day stays are left out, so only stays of 2 days or more
> count. **Low-risk tier = every other unvaccinated student**, so any
> policy with low-risk days sends every unvaccinated student home.

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
SEED_BUR <- SEED + 77L                                       # own seeds, independent of 01/03/09
VACCS    <- c(0.60, 0.70, 0.80, 0.90)
POLS     <- c("7/0/0", "21/0/0", "21/14/0", "21/14/7", "21/7/7", "21/21/21")
REF_POL  <- "7/0/0"           # reference for "cases prevented" (minimal quarantine)
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

## Is it worth it? Extra days at home per case prevented, vs the reference policy
ref <- SUM[pol == REF_POL, .(school, vacc, ref_cases = cases, ref_days = days_home, ref_u = u_days_m)]
SUM <- ref[SUM, on = .(school, vacc)]
SUM[, prevented := ref_cases - cases]
SUM[, extra_u := u_days_m - ref_u]            # extra school days per unvaccinated student vs reference
SUM[, per_case := fifelse(pol == REF_POL, NA_real_,
                          fifelse(prevented > 0.1, (days_home - ref_days) / prevented, NA_real_))]
fwrite(SUM, file.path(RESULTS_DIR, sprintf("10_burden_summary_n%d.csv", N_RUNS)))

## Per tier, per outbreak (computed in each run, averaged over runs)
TIERS[, tier := factor(tier, c("High", "Medium", "Low"))]
TIERS[, nu := round(vapply(as.character(school), function(x) sum(ENROLLMENT[[x]]), 0) * (1 - vacc))]
TSUM <- TIERS[, .(students = mean(students), stays = mean(stays), days = mean(days),
                  days_u = mean(days * 5 / 7 / nu),             # school days per unvaccinated student
                  t_students = fmt_ci(students, 0), t_days_u = fmt_ci(days * 5 / 7 / nu),
                  per_student = mean(fifelse(students > 0, days / students * 5 / 7, NA_real_), na.rm = TRUE),
                  inf_share   = mean(fifelse(stays > 0, inf / stays, NA_real_), na.rm = TRUE)),
              by = .(school, vacc, vl, pol, tier)]
TSUM[is.nan(per_student), per_student := NA][is.nan(inf_share), inf_share := NA]
fwrite(TSUM, file.path(RESULTS_DIR, sprintf("10_burden_by_tier_n%d.csv", N_RUNS)))
```

## The trade-off in one figure

Each point is one policy. **Down = safer, left = fewer school days
missed.** The x-axis is the school days an unvaccinated student misses
to quarantine in one outbreak, on average: all students’ days at home
divided by the number of unvaccinated students (only they can be
quarantined). The 100 simulated days are about 71 school days, so this
can never exceed 71. Dividing by a fixed number per school does not
change which policy costs more. A policy is **worth its cost** when no
cheaper policy comes within 2 percentage points of its risk; the blue
line joins those policies. Hollow points cost more than a policy that is
about as safe.

``` r
TOL <- 0.02     # "about as safe" = risk at most 2 points higher (as MIN_GAIN in file 03)
tradeoff_plot <- function(v) {
  TR <- RUNS[vacc == v, {
    nu  <- round(sum(ENROLLMENT[[as.character(school)]]) * (1 - v))   # unvaccinated students
    sd_ <- days_home * 5 / 7 / nu; p <- mean(size > K); n <- .N      # school days per unvaccinated student
    z <- 1.96; d <- 1 + z^2 / n; c <- (p + z^2 / (2 * n)) / d
    h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d          # Wilson interval
    .(cost = mean(sd_), cost_lo = mean(sd_) - z * sd(sd_) / sqrt(n), cost_hi = mean(sd_) + z * sd(sd_) / sqrt(n),
      risk = p, risk_lo = c - h, risk_hi = c + h)
  }, by = .(school, pol)]
  TR[, worth := vapply(seq_len(.N), function(i) !any(cost < cost[i] & risk <= risk[i] + TOL), TRUE), by = school]
  TR[, status := factor(fifelse(worth, "Worth it: no cheaper policy is as safe",
                                "Not worth it: a cheaper policy is about as safe"),
                        c("Worth it: no cheaper policy is as safe", "Not worth it: a cheaper policy is about as safe"))]
  TR[, rng := diff(range(cost)), by = school]
  TR[, below := vapply(seq_len(.N), function(i)
    any(abs(cost - cost[i]) < 0.12 * rng[i] & risk > risk[i] & risk - risk[i] < 0.09), TRUE), by = school]
  setorder(TR, school, cost)
  ggplot(TR, aes(cost, risk)) +
    geom_line(data = TR[worth == TRUE], aes(group = school), colour = "#3A7DBA", linewidth = 1) +
    geom_errorbar(aes(ymin = risk_lo, ymax = risk_hi), width = 0, colour = "grey55", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = cost_lo, xmax = cost_hi), height = 0, colour = "grey55", linewidth = 0.4) +
    geom_point(aes(shape = status, fill = status), size = 3.6, colour = "#1F4E79", stroke = 0.9) +
    geom_text(aes(label = pol, colour = status, vjust = ifelse(below, 1.9, -1.05)),
              size = 3.5, fontface = "bold", show.legend = FALSE) +
    facet_wrap(~ school, nrow = 1, scales = "free_x") +
    scale_shape_manual(values = c(21, 21), name = NULL, drop = FALSE) +
    scale_fill_manual(values = c("#3A7DBA", "white"), name = NULL, drop = FALSE) +
    scale_colour_manual(values = c("#1F4E79", "grey50")) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25)) +
    coord_cartesian(ylim = c(-0.07, 1.04), clip = "off") +      # room for labels under points near 0
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
    labs(x = "School days missed per unvaccinated student in one outbreak (quarantine)",
         y = sprintf("Chance of more than %d cases", K),
         title = "The quarantine trade-off: protection bought per school day lost",
         subtitle = sprintf(paste0("%s vaccinated. Each point is one policy (high / medium / low-risk quarantine days). ",
                                   "Down = safer, left = fewer school days missed.\n",
                                   "School days = calendar days x 5/7 (the 100 simulated days = 71 school days). Bars = 95%% CI.\n",
                                   "Blue line = policies worth their cost (no cheaper policy is within %.0f points of their risk)."),
                            vacc_lab(v), 100 * TOL)) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title.position = "plot",
          plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey30"), panel.spacing.x = unit(1.6, "lines"))
}
```

<div class="panel-tabset">

### 60% vaccinated

### 70% vaccinated

![](md_figures/10_quarantine_burden_dur-tradeoff-1.png)

### 80% vaccinated

![](md_figures/10_quarantine_burden_dur-tradeoff-2.png)

### 90% vaccinated

![](md_figures/10_quarantine_burden_dur-tradeoff-3.png)![](md_figures/10_quarantine_burden_dur-tradeoff-4.png)

</div>

## School cards: the nine answers

One table per school and vaccination level. Columns are the policies
(high / medium / low-risk quarantine days). **Every number is for one
outbreak**: it is computed in each run and then averaged over the runs
(never added up across runs). Numbers in brackets are 95% confidence
intervals of that average. Shares and days per stay are averaged over
the runs in which someone was sent home.

<div class="panel-tabset">

### Elementary

**Elementary school, 60% vaccinated** (339 students, 136 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 22%" />
<col style="width: 13%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">97%</td>
<td style="text-align: right;">67%</td>
<td style="text-align: right;">48%</td>
<td style="text-align: right;">38%</td>
<td style="text-align: right;">61%</td>
<td style="text-align: right;">31%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">130 (129131)</td>
<td style="text-align: right;">115 (113118)</td>
<td style="text-align: right;">111 (109114)</td>
<td style="text-align: right;">133 (132134)</td>
<td style="text-align: right;">134 (133134)</td>
<td style="text-align: right;">133 (132134)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">9.1 (8.99.2)</td>
<td style="text-align: right;">4.5 (4.44.7)</td>
<td style="text-align: right;">4.0 (3.94.2)</td>
<td style="text-align: right;">3.2 (3.13.3)</td>
<td style="text-align: right;">3.6 (3.43.7)</td>
<td style="text-align: right;">1.2 (1.11.2)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">9% / 14% / 77%</td>
<td style="text-align: right;">44% / 22% / 34%</td>
<td style="text-align: right;">59% / 21% / 20%</td>
<td style="text-align: right;">23% / 42% / 35%</td>
<td style="text-align: right;">16% / 36% / 48%</td>
<td style="text-align: right;">85% / 14% / 1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">16.8 (16.517.1); healthy 14.3, infected
2.5</td>
<td style="text-align: right;">11.7 (11.312.1); healthy 10.8, infected
0.9</td>
<td style="text-align: right;">10.9 (10.511.3); healthy 10.2, infected
0.7</td>
<td style="text-align: right;">17.5 (17.017.9); healthy 16.8, infected
0.6</td>
<td style="text-align: right;">20.8 (20.321.3); healthy 19.9, infected
0.9</td>
<td style="text-align: right;">18.0 (17.618.4); healthy 17.5, infected
0.5</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">4.3 (4.34.4) / 17.4 (17.117.7)</td>
<td style="text-align: right;">7.3 (7.17.5) / 14.0 (13.714.4)</td>
<td style="text-align: right;">8.7 (8.58.9) / 13.5 (13.213.8)</td>
<td style="text-align: right;">7.9 (7.78.1) / 17.9 (17.418.3)</td>
<td style="text-align: right;">7.9 (7.78.1) / 21.2 (20.621.7)</td>
<td style="text-align: right;">16.0 (15.816.3) / 18.3 (17.918.7)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">127 (126128)</td>
<td style="text-align: right;">114 (111116)</td>
<td style="text-align: right;">109 (106111)</td>
<td style="text-align: right;">133 (132134)</td>
<td style="text-align: right;">134 (133134)</td>
<td style="text-align: right;">133 (132134)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">10.4 (10.310.5)</td>
<td style="text-align: right;">7.2 (7.07.5)</td>
<td style="text-align: right;">6.5 (6.36.7)</td>
<td style="text-align: right;">5.6 (5.55.8)</td>
<td style="text-align: right;">6.6 (6.46.8)</td>
<td style="text-align: right;">3.9 (3.84.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">-5.1 days; 37.1 cases prevented</td>
<td style="text-align: right;">-5.9 days; 44.1 cases prevented</td>
<td style="text-align: right;">+0.7 days; 46.5 cases prevented</td>
<td style="text-align: right;">+4.0 days; 40.8 cases prevented</td>
<td style="text-align: right;">+1.2 days; 48.8 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">13%</td>
<td style="text-align: right;">7%</td>
<td style="text-align: right;">6%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">105.3 students, 13 school days each; 23%
infected</td>
<td style="text-align: right;">72.6 students, 16 school days each; 16%
infected</td>
<td style="text-align: right;">61.6 students, 16 school days each; 12%
infected</td>
<td style="text-align: right;">41.8 students, 15 school days each; 14%
infected</td>
<td style="text-align: right;">49.7 students, 15 school days each; 17%
infected</td>
<td style="text-align: right;">25.5 students, 15 school days each; 17%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">22.4 students, 10 school days each; 10%
infected</td>
<td style="text-align: right;">15.6 students, 10 school days each; 12%
infected</td>
<td style="text-align: right;">19.2 students, 5 school days each; 11%
infected</td>
<td style="text-align: right;">13.7 students, 15 school days each; 11%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">112.4 students, 14 school days each; 1%
infected</td>
<td style="text-align: right;">120.0 students, 16 school days each; 2%
infected</td>
<td style="text-align: right;">102.2 students, 18 school days each; 1%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">6.8</td>
<td style="text-align: right;">3.1</td>
<td style="text-align: right;">2.0</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0.0</td>
</tr>
</tbody>
</table>

**Elementary school, 70% vaccinated** (339 students, 102 unvaccinated)

<table>
<colgroup>
<col style="width: 22%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">91%</td>
<td style="text-align: right;">53%</td>
<td style="text-align: right;">34%</td>
<td style="text-align: right;">27%</td>
<td style="text-align: right;">44%</td>
<td style="text-align: right;">22%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">92 (9093)</td>
<td style="text-align: right;">78 (7680)</td>
<td style="text-align: right;">73 (7175)</td>
<td style="text-align: right;">99 (98100)</td>
<td style="text-align: right;">99 (98100)</td>
<td style="text-align: right;">99 (98100)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">8.3 (8.28.5)</td>
<td style="text-align: right;">3.7 (3.63.8)</td>
<td style="text-align: right;">3.3 (3.13.4)</td>
<td style="text-align: right;">3.0 (2.93.1)</td>
<td style="text-align: right;">3.3 (3.23.4)</td>
<td style="text-align: right;">1.1 (1.11.2)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">16% / 16% / 68%</td>
<td style="text-align: right;">58% / 21% / 22%</td>
<td style="text-align: right;">70% / 19% / 11%</td>
<td style="text-align: right;">25% / 46% / 29%</td>
<td style="text-align: right;">18% / 37% / 44%</td>
<td style="text-align: right;">87% / 12% / 1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">12.6 (12.312.9); healthy 10.7, infected
1.9</td>
<td style="text-align: right;">8.9 (8.69.3); healthy 8.2, infected
0.8</td>
<td style="text-align: right;">8.4 (8.18.7); healthy 7.8, infected
0.6</td>
<td style="text-align: right;">15.8 (15.416.3); healthy 15.3, infected
0.5</td>
<td style="text-align: right;">18.3 (17.818.8); healthy 17.5, infected
0.8</td>
<td style="text-align: right;">17.3 (16.917.7); healthy 16.9, infected
0.5</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">3.9 (3.84.0) / 13.6 (13.413.9)</td>
<td style="text-align: right;">8.1 (7.88.4) / 12.4 (12.112.7)</td>
<td style="text-align: right;">9.4 (9.29.7) / 12.4 (12.212.7)</td>
<td style="text-align: right;">7.5 (7.47.7) / 16.4 (16.016.8)</td>
<td style="text-align: right;">7.3 (7.27.5) / 18.8 (18.319.3)</td>
<td style="text-align: right;">15.9 (15.616.1) / 17.9 (17.518.3)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">88 (8790)</td>
<td style="text-align: right;">76 (7478)</td>
<td style="text-align: right;">71 (6873)</td>
<td style="text-align: right;">99 (98100)</td>
<td style="text-align: right;">99 (98100)</td>
<td style="text-align: right;">99 (98100)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">9.8 (9.610.0)</td>
<td style="text-align: right;">6.7 (6.46.9)</td>
<td style="text-align: right;">5.7 (5.55.8)</td>
<td style="text-align: right;">5.3 (5.15.4)</td>
<td style="text-align: right;">6.0 (5.86.1)</td>
<td style="text-align: right;">3.8 (3.73.9)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">-3.7 days; 25.9 cases prevented</td>
<td style="text-align: right;">-4.2 days; 30.8 cases prevented</td>
<td style="text-align: right;">+3.2 days; 32.1 cases prevented</td>
<td style="text-align: right;">+5.7 days; 28.3 cases prevented</td>
<td style="text-align: right;">+4.7 days; 33.2 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">13%</td>
<td style="text-align: right;">8%</td>
<td style="text-align: right;">7%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">68.1 students, 12 school days each; 22%
infected</td>
<td style="text-align: right;">45.5 students, 15 school days each; 17%
infected</td>
<td style="text-align: right;">37.3 students, 15 school days each; 13%
infected</td>
<td style="text-align: right;">28.4 students, 15 school days each; 14%
infected</td>
<td style="text-align: right;">32.7 students, 15 school days each; 17%
infected</td>
<td style="text-align: right;">18.6 students, 15 school days each; 17%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">14.8 students, 10 school days each; 10%
infected</td>
<td style="text-align: right;">11.5 students, 10 school days each; 10%
infected</td>
<td style="text-align: right;">14.3 students, 5 school days each; 10%
infected</td>
<td style="text-align: right;">10.0 students, 15 school days each; 11%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">82.6 students, 13 school days each; 1%
infected</td>
<td style="text-align: right;">87.8 students, 15 school days each; 2%
infected</td>
<td style="text-align: right;">75.5 students, 17 school days each; 1%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">4.7</td>
<td style="text-align: right;">2.0</td>
<td style="text-align: right;">1.2</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.0</td>
</tr>
</tbody>
</table>

**Elementary school, 80% vaccinated** (339 students, 68 unvaccinated)

<table>
<colgroup>
<col style="width: 22%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">74%</td>
<td style="text-align: right;">31%</td>
<td style="text-align: right;">17%</td>
<td style="text-align: right;">14%</td>
<td style="text-align: right;">27%</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">51 (4952)</td>
<td style="text-align: right;">38 (3740)</td>
<td style="text-align: right;">40 (3841)</td>
<td style="text-align: right;">64 (6465)</td>
<td style="text-align: right;">64 (6365)</td>
<td style="text-align: right;">65 (6466)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">6.2 (6.06.4)</td>
<td style="text-align: right;">2.6 (2.52.7)</td>
<td style="text-align: right;">2.5 (2.32.6)</td>
<td style="text-align: right;">2.6 (2.62.7)</td>
<td style="text-align: right;">2.8 (2.72.9)</td>
<td style="text-align: right;">1.1 (1.11.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">31% / 24% / 45%</td>
<td style="text-align: right;">77% / 16% / 8%</td>
<td style="text-align: right;">82% / 14% / 4%</td>
<td style="text-align: right;">27% / 48% / 25%</td>
<td style="text-align: right;">22% / 40% / 38%</td>
<td style="text-align: right;">89% / 11% / 0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">7.3 (7.17.6); healthy 6.2, infected
1.1</td>
<td style="text-align: right;">6.0 (5.86.3); healthy 5.5, infected
0.5</td>
<td style="text-align: right;">6.4 (6.26.7); healthy 5.9, infected
0.5</td>
<td style="text-align: right;">13.9 (13.514.3); healthy 13.4, infected
0.5</td>
<td style="text-align: right;">15.4 (14.915.9); healthy 14.8, infected
0.6</td>
<td style="text-align: right;">16.2 (15.816.6); healthy 15.8, infected
0.4</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">3.9 (3.83.9) / 9.6 (9.39.8)</td>
<td style="text-align: right;">10.3 (10.010.6) / 12.3 (12.112.6)</td>
<td style="text-align: right;">10.4 (10.210.7) / 12.1 (11.912.4)</td>
<td style="text-align: right;">7.0 (6.97.1) / 14.7 (14.315.0)</td>
<td style="text-align: right;">6.6 (6.56.7) / 16.3 (15.816.7)</td>
<td style="text-align: right;">15.4 (15.215.6) / 17.0 (16.717.3)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">47 (4649)</td>
<td style="text-align: right;">37 (3539)</td>
<td style="text-align: right;">38 (3639)</td>
<td style="text-align: right;">64 (6465)</td>
<td style="text-align: right;">64 (6365)</td>
<td style="text-align: right;">65 (6466)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">8.6 (8.48.8)</td>
<td style="text-align: right;">5.6 (5.45.8)</td>
<td style="text-align: right;">4.9 (4.75.1)</td>
<td style="text-align: right;">4.7 (4.64.9)</td>
<td style="text-align: right;">5.3 (5.15.5)</td>
<td style="text-align: right;">3.5 (3.43.6)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">-1.3 days; 13.2 cases prevented</td>
<td style="text-align: right;">-0.9 days; 15.6 cases prevented</td>
<td style="text-align: right;">+6.5 days; 15.9 cases prevented</td>
<td style="text-align: right;">+8.0 days; 14.0 cases prevented</td>
<td style="text-align: right;">+8.8 days; 16.7 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">12%</td>
<td style="text-align: right;">9%</td>
<td style="text-align: right;">7%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">33.7 students, 10 school days each; 20%
infected</td>
<td style="text-align: right;">22.6 students, 15 school days each; 16%
infected</td>
<td style="text-align: right;">20.0 students, 15 school days each; 13%
infected</td>
<td style="text-align: right;">16.8 students, 15 school days each; 14%
infected</td>
<td style="text-align: right;">18.6 students, 15 school days each; 16%
infected</td>
<td style="text-align: right;">11.7 students, 15 school days each; 15%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">8.9 students, 10 school days each; 10%
infected</td>
<td style="text-align: right;">7.6 students, 10 school days each; 10%
infected</td>
<td style="text-align: right;">9.2 students, 5 school days each; 10%
infected</td>
<td style="text-align: right;">6.5 students, 15 school days each; 10%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">53.1 students, 11 school days each; 1%
infected</td>
<td style="text-align: right;">55.9 students, 13 school days each; 1%
infected</td>
<td style="text-align: right;">49.1 students, 17 school days each; 1%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">2.2</td>
<td style="text-align: right;">0.9</td>
<td style="text-align: right;">0.6</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">—</td>
</tr>
</tbody>
</table>

**Elementary school, 90% vaccinated** (339 students, 34 unvaccinated)

<table>
<colgroup>
<col style="width: 22%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">23%</td>
<td style="text-align: right;">7%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">13 (1214)</td>
<td style="text-align: right;">11 (1112)</td>
<td style="text-align: right;">12 (1213)</td>
<td style="text-align: right;">30 (2931)</td>
<td style="text-align: right;">30 (2931)</td>
<td style="text-align: right;">30 (3031)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">2.7 (2.62.9)</td>
<td style="text-align: right;">1.6 (1.51.6)</td>
<td style="text-align: right;">1.6 (1.51.6)</td>
<td style="text-align: right;">1.9 (1.82.0)</td>
<td style="text-align: right;">1.9 (1.82.0)</td>
<td style="text-align: right;">1.0 (0.91.0)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">65% / 24% / 11%</td>
<td style="text-align: right;">92% / 7% / 1%</td>
<td style="text-align: right;">93% / 6% / 1%</td>
<td style="text-align: right;">43% / 42% / 15%</td>
<td style="text-align: right;">42% / 39% / 19%</td>
<td style="text-align: right;">92% / 7% / 0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">2.4 (2.32.6); healthy 2.1, infected
0.3</td>
<td style="text-align: right;">3.7 (3.53.9); healthy 3.3, infected
0.3</td>
<td style="text-align: right;">4.4 (4.24.6); healthy 4.0, infected
0.4</td>
<td style="text-align: right;">10.7 (10.311.0); healthy 10.3, infected
0.4</td>
<td style="text-align: right;">10.4 (10.010.8); healthy 10.1, infected
0.4</td>
<td style="text-align: right;">14.5 (14.114.9); healthy 14.2, infected
0.3</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">4.6 (4.54.6) / 6.7 (6.56.9)</td>
<td style="text-align: right;">12.6 (12.312.8) / 13.4 (13.113.6)</td>
<td style="text-align: right;">12.1 (11.912.3) / 12.8 (12.613.0)</td>
<td style="text-align: right;">6.9 (6.86.9) / 12.1 (11.812.4)</td>
<td style="text-align: right;">6.4 (6.36.5) / 11.8 (11.512.2)</td>
<td style="text-align: right;">15.2 (15.015.3) / 16.4 (16.116.7)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">11 (1112)</td>
<td style="text-align: right;">11 (1011)</td>
<td style="text-align: right;">12 (1112)</td>
<td style="text-align: right;">30 (2931)</td>
<td style="text-align: right;">30 (2931)</td>
<td style="text-align: right;">30 (3031)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">4.5 (4.24.7)</td>
<td style="text-align: right;">4.1 (3.94.3)</td>
<td style="text-align: right;">3.8 (3.63.9)</td>
<td style="text-align: right;">3.9 (3.84.1)</td>
<td style="text-align: right;">3.9 (3.74.1)</td>
<td style="text-align: right;">3.1 (3.03.2)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">+1.3 days; 2.3 cases prevented</td>
<td style="text-align: right;">+1.9 days; 2.9 cases prevented</td>
<td style="text-align: right;">+8.2 days; 2.7 cases prevented</td>
<td style="text-align: right;">+8.0 days; 2.7 cases prevented</td>
<td style="text-align: right;">+12.1 days; 3.2 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">12%</td>
<td style="text-align: right;">12%</td>
<td style="text-align: right;">10%</td>
<td style="text-align: right;">2%</td>
<td style="text-align: right;">2%</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">8.8 students, 8 school days each; 15%
infected</td>
<td style="text-align: right;">7.8 students, 15 school days each; 15%
infected</td>
<td style="text-align: right;">7.2 students, 15 school days each; 14%
infected</td>
<td style="text-align: right;">6.9 students, 15 school days each; 13%
infected</td>
<td style="text-align: right;">6.7 students, 15 school days each; 14%
infected</td>
<td style="text-align: right;">5.4 students, 15 school days each; 15%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">3.6 students, 10 school days each; 10%
infected</td>
<td style="text-align: right;">3.6 students, 10 school days each; 9%
infected</td>
<td style="text-align: right;">3.8 students, 5 school days each; 8%
infected</td>
<td style="text-align: right;">2.9 students, 15 school days each; 11%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">23.8 students, 9 school days each; 1%
infected</td>
<td style="text-align: right;">24.7 students, 9 school days each; 1%
infected</td>
<td style="text-align: right;">22.7 students, 16 school days each; 1%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">0.3</td>
<td style="text-align: right;">0.2</td>
<td style="text-align: right;">0.2</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">—</td>
</tr>
</tbody>
</table>

### Middle

**Middle school, 60% vaccinated** (466 students, 186 unvaccinated)

<table>
<colgroup>
<col style="width: 22%" />
<col style="width: 13%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 13%" />
<col style="width: 13%" />
<col style="width: 13%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">77%</td>
<td style="text-align: right;">16%</td>
<td style="text-align: right;">16%</td>
<td style="text-align: right;">15%</td>
<td style="text-align: right;">16%</td>
<td style="text-align: right;">15%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">155 (151159)</td>
<td style="text-align: right;">118 (114122)</td>
<td style="text-align: right;">121 (117125)</td>
<td style="text-align: right;">180 (178182)</td>
<td style="text-align: right;">182 (180183)</td>
<td style="text-align: right;">180 (178182)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">5.4 (5.25.5)</td>
<td style="text-align: right;">1.7 (1.71.8)</td>
<td style="text-align: right;">1.8 (1.81.9)</td>
<td style="text-align: right;">2.2 (2.22.3)</td>
<td style="text-align: right;">2.2 (2.22.3)</td>
<td style="text-align: right;">1.1 (1.01.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">22% / 21% / 57%</td>
<td style="text-align: right;">82% / 13% / 5%</td>
<td style="text-align: right;">82% / 13% / 5%</td>
<td style="text-align: right;">35% / 49% / 17%</td>
<td style="text-align: right;">31% / 51% / 18%</td>
<td style="text-align: right;">91% / 8% / 1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">9.2 (8.99.5); healthy 8.5, infected
0.7</td>
<td style="text-align: right;">5.7 (5.55.8); healthy 5.5, infected
0.2</td>
<td style="text-align: right;">5.7 (5.55.9); healthy 5.5, infected
0.2</td>
<td style="text-align: right;">13.2 (12.913.5); healthy 13.0, infected
0.2</td>
<td style="text-align: right;">13.5 (13.213.8); healthy 13.3, infected
0.2</td>
<td style="text-align: right;">16.6 (16.217.0); healthy 16.4, infected
0.2</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">3.9 (3.84.0) / 11.1 (10.811.3)</td>
<td style="text-align: right;">10.1 (9.710.4) / 11.3 (10.911.6)</td>
<td style="text-align: right;">9.9 (9.510.2) / 11.0 (10.711.4)</td>
<td style="text-align: right;">7.6 (7.47.7) / 13.7 (13.414.0)</td>
<td style="text-align: right;">7.3 (7.27.5) / 13.9 (13.614.1)</td>
<td style="text-align: right;">15.7 (15.515.9) / 17.1 (16.817.5)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">151 (147155)</td>
<td style="text-align: right;">117 (112121)</td>
<td style="text-align: right;">120 (116124)</td>
<td style="text-align: right;">180 (178182)</td>
<td style="text-align: right;">182 (180183)</td>
<td style="text-align: right;">180 (178182)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">7.8 (7.68.0)</td>
<td style="text-align: right;">3.6 (3.53.7)</td>
<td style="text-align: right;">3.6 (3.53.7)</td>
<td style="text-align: right;">3.7 (3.63.9)</td>
<td style="text-align: right;">3.8 (3.73.9)</td>
<td style="text-align: right;">3.4 (3.43.5)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">-3.5 days; 14.9 cases prevented</td>
<td style="text-align: right;">-3.5 days; 14.6 cases prevented</td>
<td style="text-align: right;">+4.1 days; 14.6 cases prevented</td>
<td style="text-align: right;">+4.4 days; 14.7 cases prevented</td>
<td style="text-align: right;">+7.5 days; 15.0 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">59.5 students, 18 school days each; 14%
infected</td>
<td style="text-align: right;">53.2 students, 16 school days each; 8%
infected</td>
<td style="text-align: right;">52.1 students, 16 school days each; 9%
infected</td>
<td style="text-align: right;">50.6 students, 16 school days each; 9%
infected</td>
<td style="text-align: right;">52.0 students, 16 school days each; 8%
infected</td>
<td style="text-align: right;">48.4 students, 16 school days each; 9%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">0.9 students, 10 school days each; 5%
infected</td>
<td style="text-align: right;">0.8 students, 10 school days each; 5%
infected</td>
<td style="text-align: right;">0.8 students, 5 school days each; 5%
infected</td>
<td style="text-align: right;">0.7 students, 15 school days each; 4%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">131.7 students, 12 school days each; 0%
infected</td>
<td style="text-align: right;">134.3 students, 12 school days each; 0%
infected</td>
<td style="text-align: right;">131.5 students, 17 school days each; 0%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">3.3</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">1.1</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.0</td>
</tr>
</tbody>
</table>

**Middle school, 70% vaccinated** (466 students, 140 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 22%" />
<col style="width: 13%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">57%</td>
<td style="text-align: right;">9%</td>
<td style="text-align: right;">9%</td>
<td style="text-align: right;">7%</td>
<td style="text-align: right;">10%</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">97 (93100)</td>
<td style="text-align: right;">69 (6673)</td>
<td style="text-align: right;">73 (7076)</td>
<td style="text-align: right;">132 (130134)</td>
<td style="text-align: right;">132 (130134)</td>
<td style="text-align: right;">132 (130134)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">4.2 (4.04.3)</td>
<td style="text-align: right;">1.5 (1.41.5)</td>
<td style="text-align: right;">1.5 (1.51.6)</td>
<td style="text-align: right;">2.1 (2.02.1)</td>
<td style="text-align: right;">2.0 (2.02.1)</td>
<td style="text-align: right;">1.0 (1.01.0)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">34% / 23% / 42%</td>
<td style="text-align: right;">87% / 10% / 3%</td>
<td style="text-align: right;">89% / 9% / 3%</td>
<td style="text-align: right;">36% / 53% / 12%</td>
<td style="text-align: right;">36% / 51% / 13%</td>
<td style="text-align: right;">93% / 7% / 0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">6.3 (6.16.6); healthy 5.8, infected
0.5</td>
<td style="text-align: right;">4.9 (4.75.0); healthy 4.7, infected
0.2</td>
<td style="text-align: right;">4.9 (4.85.1); healthy 4.7, infected
0.2</td>
<td style="text-align: right;">11.9 (11.612.2); healthy 11.7, infected
0.2</td>
<td style="text-align: right;">12.2 (11.812.5); healthy 12.0, infected
0.2</td>
<td style="text-align: right;">15.9 (15.516.3); healthy 15.7, infected
0.2</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">4.1 (4.04.2) / 9.8 (9.510.1)</td>
<td style="text-align: right;">11.5 (11.211.8) / 12.5 (12.212.9)</td>
<td style="text-align: right;">11.4 (11.011.7) / 12.2 (11.912.6)</td>
<td style="text-align: right;">7.2 (7.17.3) / 12.7 (12.412.9)</td>
<td style="text-align: right;">7.3 (7.27.4) / 12.9 (12.713.2)</td>
<td style="text-align: right;">15.6 (15.415.8) / 16.8 (16.517.2)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">93 (9097)</td>
<td style="text-align: right;">69 (6672)</td>
<td style="text-align: right;">72 (6975)</td>
<td style="text-align: right;">132 (130134)</td>
<td style="text-align: right;">132 (130134)</td>
<td style="text-align: right;">132 (130134)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">6.7 (6.46.9)</td>
<td style="text-align: right;">3.4 (3.33.5)</td>
<td style="text-align: right;">3.4 (3.33.5)</td>
<td style="text-align: right;">3.4 (3.33.5)</td>
<td style="text-align: right;">3.6 (3.43.7)</td>
<td style="text-align: right;">3.3 (3.23.4)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">-1.5 days; 8.6 cases prevented</td>
<td style="text-align: right;">-1.4 days; 8.7 cases prevented</td>
<td style="text-align: right;">+5.6 days; 8.8 cases prevented</td>
<td style="text-align: right;">+5.8 days; 8.6 cases prevented</td>
<td style="text-align: right;">+9.6 days; 8.7 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">40.5 students, 16 school days each; 12%
infected</td>
<td style="text-align: right;">36.7 students, 16 school days each; 8%
infected</td>
<td style="text-align: right;">37.4 students, 16 school days each; 8%
infected</td>
<td style="text-align: right;">36.2 students, 16 school days each; 8%
infected</td>
<td style="text-align: right;">37.6 students, 16 school days each; 8%
infected</td>
<td style="text-align: right;">36.1 students, 16 school days each; 8%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">0.7 students, 10 school days each; 5%
infected</td>
<td style="text-align: right;">0.6 students, 10 school days each; 6%
infected</td>
<td style="text-align: right;">0.6 students, 5 school days each; 5%
infected</td>
<td style="text-align: right;">0.5 students, 15 school days each; 6%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">96.3 students, 11 school days each; 0%
infected</td>
<td style="text-align: right;">96.9 students, 11 school days each; 0%
infected</td>
<td style="text-align: right;">95.9 students, 17 school days each; 0%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">1.7</td>
<td style="text-align: right;">0.6</td>
<td style="text-align: right;">0.6</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.0</td>
</tr>
</tbody>
</table>

**Middle school, 80% vaccinated** (466 students, 93 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 22%" />
<col style="width: 13%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">32%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">46 (4448)</td>
<td style="text-align: right;">37 (3539)</td>
<td style="text-align: right;">37 (3539)</td>
<td style="text-align: right;">84 (8286)</td>
<td style="text-align: right;">84 (8286)</td>
<td style="text-align: right;">84 (8386)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">2.9 (2.83.1)</td>
<td style="text-align: right;">1.2 (1.11.2)</td>
<td style="text-align: right;">1.2 (1.21.3)</td>
<td style="text-align: right;">1.8 (1.71.8)</td>
<td style="text-align: right;">1.7 (1.71.8)</td>
<td style="text-align: right;">1.0 (1.01.0)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">47% / 24% / 29%</td>
<td style="text-align: right;">93% / 6% / 1%</td>
<td style="text-align: right;">94% / 5% / 1%</td>
<td style="text-align: right;">44% / 47% / 9%</td>
<td style="text-align: right;">44% / 47% / 9%</td>
<td style="text-align: right;">93% / 7% / 1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">4.0 (3.84.2); healthy 3.7, infected
0.3</td>
<td style="text-align: right;">4.3 (4.14.4); healthy 4.1, infected
0.2</td>
<td style="text-align: right;">4.2 (4.14.4); healthy 4.1, infected
0.2</td>
<td style="text-align: right;">10.6 (10.310.9); healthy 10.4, infected
0.2</td>
<td style="text-align: right;">10.5 (10.210.8); healthy 10.4, infected
0.2</td>
<td style="text-align: right;">14.9 (14.515.3); healthy 14.7, infected
0.2</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">4.5 (4.44.6) / 9.0 (8.69.3)</td>
<td style="text-align: right;">12.9 (12.513.2) / 13.5 (13.113.8)</td>
<td style="text-align: right;">12.7 (12.413.0) / 13.1 (12.813.4)</td>
<td style="text-align: right;">7.2 (7.17.3) / 11.7 (11.511.9)</td>
<td style="text-align: right;">7.1 (7.07.2) / 11.7 (11.411.9)</td>
<td style="text-align: right;">15.3 (15.115.4) / 16.4 (16.116.7)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">45 (4247)</td>
<td style="text-align: right;">37 (3539)</td>
<td style="text-align: right;">37 (3539)</td>
<td style="text-align: right;">84 (8286)</td>
<td style="text-align: right;">84 (8286)</td>
<td style="text-align: right;">84 (8386)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">5.0 (4.85.3)</td>
<td style="text-align: right;">3.1 (3.03.2)</td>
<td style="text-align: right;">3.0 (2.93.1)</td>
<td style="text-align: right;">3.2 (3.13.3)</td>
<td style="text-align: right;">3.2 (3.13.3)</td>
<td style="text-align: right;">3.1 (3.03.2)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">+0.2 days; 4.2 cases prevented</td>
<td style="text-align: right;">+0.2 days; 4.3 cases prevented</td>
<td style="text-align: right;">+6.6 days; 4.2 cases prevented</td>
<td style="text-align: right;">+6.5 days; 4.3 cases prevented</td>
<td style="text-align: right;">+10.9 days; 4.1 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">24.2 students, 13 school days each; 10%
infected</td>
<td style="text-align: right;">23.0 students, 16 school days each; 8%
infected</td>
<td style="text-align: right;">23.3 students, 15 school days each; 8%
infected</td>
<td style="text-align: right;">23.2 students, 16 school days each; 8%
infected</td>
<td style="text-align: right;">22.9 students, 16 school days each; 8%
infected</td>
<td style="text-align: right;">22.8 students, 16 school days each; 8%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">0.4 students, 10 school days each; 6%
infected</td>
<td style="text-align: right;">0.3 students, 10 school days each; 4%
infected</td>
<td style="text-align: right;">0.4 students, 5 school days each; 2%
infected</td>
<td style="text-align: right;">0.3 students, 15 school days each; 4%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">61.2 students, 10 school days each; 0%
infected</td>
<td style="text-align: right;">61.6 students, 10 school days each; 0%
infected</td>
<td style="text-align: right;">61.4 students, 17 school days each; 0%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">0.7</td>
<td style="text-align: right;">0.3</td>
<td style="text-align: right;">0.3</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">0.0</td>
</tr>
</tbody>
</table>

**Middle school, 90% vaccinated** (466 students, 47 unvaccinated)

<table>
<colgroup>
<col style="width: 22%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">0%</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">13 (1214)</td>
<td style="text-align: right;">12 (1213)</td>
<td style="text-align: right;">13 (1314)</td>
<td style="text-align: right;">39 (3840)</td>
<td style="text-align: right;">38 (3740)</td>
<td style="text-align: right;">40 (3841)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">1.6 (1.51.7)</td>
<td style="text-align: right;">0.9 (0.90.9)</td>
<td style="text-align: right;">1.0 (0.91.0)</td>
<td style="text-align: right;">1.4 (1.31.4)</td>
<td style="text-align: right;">1.4 (1.31.4)</td>
<td style="text-align: right;">0.9 (0.90.9)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">62% / 23% / 14%</td>
<td style="text-align: right;">96% / 4% / 0%</td>
<td style="text-align: right;">96% / 4% / 0%</td>
<td style="text-align: right;">59% / 37% / 4%</td>
<td style="text-align: right;">59% / 37% / 5%</td>
<td style="text-align: right;">96% / 4% / 0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">2.0 (1.92.1); healthy 1.9, infected
0.1</td>
<td style="text-align: right;">3.5 (3.33.6); healthy 3.3, infected
0.1</td>
<td style="text-align: right;">3.7 (3.63.9); healthy 3.6, infected
0.1</td>
<td style="text-align: right;">8.5 (8.28.8); healthy 8.3, infected
0.1</td>
<td style="text-align: right;">8.5 (8.28.8); healthy 8.3, infected
0.1</td>
<td style="text-align: right;">13.3 (12.913.7); healthy 13.2, infected
0.1</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">4.9 (4.84.9) / 7.7 (7.58.0)</td>
<td style="text-align: right;">14.0 (13.914.2) / 14.6 (14.314.9)</td>
<td style="text-align: right;">14.1 (13.914.3) / 14.6 (14.314.8)</td>
<td style="text-align: right;">7.3 (7.27.3) / 10.3 (10.110.5)</td>
<td style="text-align: right;">7.2 (7.17.2) / 10.3 (10.110.6)</td>
<td style="text-align: right;">15.1 (15.015.3) / 15.8 (15.616.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">13 (1213)</td>
<td style="text-align: right;">12 (1213)</td>
<td style="text-align: right;">13 (1214)</td>
<td style="text-align: right;">39 (3840)</td>
<td style="text-align: right;">38 (3740)</td>
<td style="text-align: right;">40 (3841)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">2.6 (2.52.8)</td>
<td style="text-align: right;">2.7 (2.62.8)</td>
<td style="text-align: right;">2.8 (2.72.9)</td>
<td style="text-align: right;">2.8 (2.72.9)</td>
<td style="text-align: right;">2.8 (2.72.9)</td>
<td style="text-align: right;">2.8 (2.72.9)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">+1.4 days; 0.9 cases prevented</td>
<td style="text-align: right;">+1.7 days; 0.9 cases prevented</td>
<td style="text-align: right;">+6.5 days; 0.9 cases prevented</td>
<td style="text-align: right;">+6.5 days; 0.9 cases prevented</td>
<td style="text-align: right;">+11.3 days; 0.9 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">6%</td>
<td style="text-align: right;">6%</td>
<td style="text-align: right;">6%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">1%</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">10.5 students, 8 school days each; 7%
infected</td>
<td style="text-align: right;">10.4 students, 15 school days each; 7%
infected</td>
<td style="text-align: right;">11.0 students, 15 school days each; 6%
infected</td>
<td style="text-align: right;">10.7 students, 15 school days each; 7%
infected</td>
<td style="text-align: right;">10.5 students, 15 school days each; 7%
infected</td>
<td style="text-align: right;">10.5 students, 15 school days each; 7%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">0.2 students, 10 school days each; 3%
infected</td>
<td style="text-align: right;">0.2 students, 10 school days each; 5%
infected</td>
<td style="text-align: right;">0.2 students, 5 school days each; 3%
infected</td>
<td style="text-align: right;">0.2 students, 15 school days each; 4%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">28.3 students, 8 school days each; 0%
infected</td>
<td style="text-align: right;">28.1 students, 8 school days each; 0%
infected</td>
<td style="text-align: right;">28.9 students, 16 school days each; 0%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
</tr>
</tbody>
</table>

### High

**High school, 60% vaccinated** (473 students, 189 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 22%" />
<col style="width: 13%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">98%</td>
<td style="text-align: right;">84%</td>
<td style="text-align: right;">57%</td>
<td style="text-align: right;">47%</td>
<td style="text-align: right;">73%</td>
<td style="text-align: right;">34%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">183 (181184)</td>
<td style="text-align: right;">177 (174179)</td>
<td style="text-align: right;">174 (172176)</td>
<td style="text-align: right;">187 (186188)</td>
<td style="text-align: right;">187 (186188)</td>
<td style="text-align: right;">188 (187188)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">8.7 (8.58.8)</td>
<td style="text-align: right;">5.1 (5.05.3)</td>
<td style="text-align: right;">4.7 (4.54.9)</td>
<td style="text-align: right;">3.5 (3.43.6)</td>
<td style="text-align: right;">4.1 (4.04.2)</td>
<td style="text-align: right;">1.2 (1.21.2)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">8% / 14% / 77%</td>
<td style="text-align: right;">36% / 28% / 35%</td>
<td style="text-align: right;">58% / 24% / 18%</td>
<td style="text-align: right;">28% / 39% / 33%</td>
<td style="text-align: right;">16% / 29% / 55%</td>
<td style="text-align: right;">82% / 16% / 1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">22.2 (21.922.5); healthy 20.0, infected
2.2</td>
<td style="text-align: right;">18.7 (18.319.2); healthy 17.9, infected
0.9</td>
<td style="text-align: right;">17.2 (16.817.7); healthy 16.6, infected
0.6</td>
<td style="text-align: right;">20.6 (20.121.2); healthy 20.1, infected
0.6</td>
<td style="text-align: right;">25.6 (25.026.1); healthy 24.7, infected
0.9</td>
<td style="text-align: right;">19.4 (18.919.8); healthy 18.9, infected
0.4</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">5.9 (5.96.0) / 22.8 (22.523.0)</td>
<td style="text-align: right;">9.7 (9.59.9) / 19.8 (19.520.2)</td>
<td style="text-align: right;">11.6 (11.411.8) / 18.4 (18.018.8)</td>
<td style="text-align: right;">9.5 (9.39.6) / 20.8 (20.321.4)</td>
<td style="text-align: right;">9.2 (9.09.4) / 25.8 (25.226.4)</td>
<td style="text-align: right;">16.6 (16.316.9) / 19.5 (19.120.0)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">180 (178181)</td>
<td style="text-align: right;">175 (173177)</td>
<td style="text-align: right;">171 (168173)</td>
<td style="text-align: right;">187 (186188)</td>
<td style="text-align: right;">187 (186188)</td>
<td style="text-align: right;">188 (187188)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">10.3 (10.110.4)</td>
<td style="text-align: right;">8.3 (8.18.5)</td>
<td style="text-align: right;">7.2 (7.17.4)</td>
<td style="text-align: right;">6.3 (6.16.5)</td>
<td style="text-align: right;">7.6 (7.47.8)</td>
<td style="text-align: right;">4.2 (4.14.4)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">-3.5 days; 40.0 cases prevented</td>
<td style="text-align: right;">-5.0 days; 50.6 cases prevented</td>
<td style="text-align: right;">-1.6 days; 53.4 cases prevented</td>
<td style="text-align: right;">+3.4 days; 45.0 cases prevented</td>
<td style="text-align: right;">-2.8 days; 57.0 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">12%</td>
<td style="text-align: right;">6%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">179.1 students, 18 school days each; 12%
infected</td>
<td style="text-align: right;">160.9 students, 19 school days each; 8%
infected</td>
<td style="text-align: right;">143.0 students, 17 school days each; 5%
infected</td>
<td style="text-align: right;">104.7 students, 16 school days each; 5%
infected</td>
<td style="text-align: right;">124.0 students, 17 school days each; 7%
infected</td>
<td style="text-align: right;">58.6 students, 16 school days each; 6%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">43.7 students, 10 school days each; 9%
infected</td>
<td style="text-align: right;">31.8 students, 10 school days each; 9%
infected</td>
<td style="text-align: right;">41.8 students, 5 school days each; 9%
infected</td>
<td style="text-align: right;">27.8 students, 15 school days each; 9%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">147.4 students, 12 school days each; 1%
infected</td>
<td style="text-align: right;">166.6 students, 15 school days each; 3%
infected</td>
<td style="text-align: right;">126.9 students, 18 school days each; 1%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">5.0</td>
<td style="text-align: right;">2.7</td>
<td style="text-align: right;">1.3</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0.0</td>
</tr>
</tbody>
</table>

**High school, 70% vaccinated** (473 students, 142 unvaccinated)

<table>
<colgroup>
<col style="width: 22%" />
<col style="width: 13%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">95%</td>
<td style="text-align: right;">73%</td>
<td style="text-align: right;">42%</td>
<td style="text-align: right;">36%</td>
<td style="text-align: right;">60%</td>
<td style="text-align: right;">26%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">133 (132135)</td>
<td style="text-align: right;">122 (120125)</td>
<td style="text-align: right;">119 (117121)</td>
<td style="text-align: right;">139 (139140)</td>
<td style="text-align: right;">138 (137139)</td>
<td style="text-align: right;">139 (138140)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">8.6 (8.48.8)</td>
<td style="text-align: right;">4.6 (4.44.7)</td>
<td style="text-align: right;">3.9 (3.84.1)</td>
<td style="text-align: right;">3.4 (3.33.5)</td>
<td style="text-align: right;">3.8 (3.74.0)</td>
<td style="text-align: right;">1.2 (1.11.2)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">10% / 13% / 77%</td>
<td style="text-align: right;">45% / 28% / 28%</td>
<td style="text-align: right;">66% / 23% / 11%</td>
<td style="text-align: right;">26% / 40% / 34%</td>
<td style="text-align: right;">17% / 29% / 54%</td>
<td style="text-align: right;">84% / 15% / 1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">19.0 (18.719.4); healthy 17.2, infected
1.9</td>
<td style="text-align: right;">16.1 (15.716.5); healthy 15.4, infected
0.7</td>
<td style="text-align: right;">14.4 (14.014.9); healthy 13.9, infected
0.5</td>
<td style="text-align: right;">19.4 (18.919.9); healthy 18.9, infected
0.5</td>
<td style="text-align: right;">22.9 (22.323.5); healthy 22.1, infected
0.8</td>
<td style="text-align: right;">18.0 (17.618.4); healthy 17.6, infected
0.4</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">5.2 (5.25.3) / 19.9 (19.720.2)</td>
<td style="text-align: right;">10.1 (9.910.3) / 18.4 (18.118.7)</td>
<td style="text-align: right;">11.8 (11.612.0) / 16.8 (16.517.2)</td>
<td style="text-align: right;">8.7 (8.68.9) / 19.8 (19.320.3)</td>
<td style="text-align: right;">8.4 (8.38.5) / 23.5 (22.924.1)</td>
<td style="text-align: right;">15.8 (15.616.1) / 18.4 (18.018.8)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">130 (128131)</td>
<td style="text-align: right;">120 (118123)</td>
<td style="text-align: right;">116 (114118)</td>
<td style="text-align: right;">139 (139140)</td>
<td style="text-align: right;">138 (137139)</td>
<td style="text-align: right;">139 (138140)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">10.2 (10.010.3)</td>
<td style="text-align: right;">7.9 (7.78.1)</td>
<td style="text-align: right;">6.5 (6.36.7)</td>
<td style="text-align: right;">6.2 (6.06.3)</td>
<td style="text-align: right;">7.1 (6.97.3)</td>
<td style="text-align: right;">3.9 (3.84.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">-2.9 days; 30.2 cases prevented</td>
<td style="text-align: right;">-4.6 days; 38.6 cases prevented</td>
<td style="text-align: right;">+0.4 days; 39.6 cases prevented</td>
<td style="text-align: right;">+3.8 days; 34.2 cases prevented</td>
<td style="text-align: right;">-1.1 days; 42.0 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">11%</td>
<td style="text-align: right;">6%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">129.1 students, 16 school days each; 12%
infected</td>
<td style="text-align: right;">109.8 students, 18 school days each; 7%
infected</td>
<td style="text-align: right;">94.0 students, 17 school days each; 4%
infected</td>
<td style="text-align: right;">73.3 students, 16 school days each; 5%
infected</td>
<td style="text-align: right;">83.6 students, 17 school days each; 7%
infected</td>
<td style="text-align: right;">41.6 students, 16 school days each; 6%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">30.1 students, 10 school days each; 9%
infected</td>
<td style="text-align: right;">24.4 students, 10 school days each; 9%
infected</td>
<td style="text-align: right;">30.5 students, 5 school days each; 8%
infected</td>
<td style="text-align: right;">21.1 students, 15 school days each; 9%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">109.7 students, 12 school days each; 1%
infected</td>
<td style="text-align: right;">120.7 students, 14 school days each; 2%
infected</td>
<td style="text-align: right;">90.8 students, 17 school days each; 1%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">4.2</td>
<td style="text-align: right;">1.9</td>
<td style="text-align: right;">0.9</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0.0</td>
</tr>
</tbody>
</table>

**High school, 80% vaccinated** (473 students, 95 unvaccinated)

<table style="width:100%;">
<colgroup>
<col style="width: 22%" />
<col style="width: 13%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">82%</td>
<td style="text-align: right;">48%</td>
<td style="text-align: right;">28%</td>
<td style="text-align: right;">23%</td>
<td style="text-align: right;">39%</td>
<td style="text-align: right;">13%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">79 (7780)</td>
<td style="text-align: right;">70 (6872)</td>
<td style="text-align: right;">71 (6972)</td>
<td style="text-align: right;">91 (8992)</td>
<td style="text-align: right;">91 (8992)</td>
<td style="text-align: right;">89 (8890)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">7.2 (7.07.4)</td>
<td style="text-align: right;">3.3 (3.23.4)</td>
<td style="text-align: right;">3.2 (3.13.4)</td>
<td style="text-align: right;">2.9 (2.83.0)</td>
<td style="text-align: right;">3.3 (3.23.4)</td>
<td style="text-align: right;">1.1 (1.11.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">20% / 16% / 64%</td>
<td style="text-align: right;">65% / 23% / 12%</td>
<td style="text-align: right;">74% / 19% / 6%</td>
<td style="text-align: right;">32% / 40% / 28%</td>
<td style="text-align: right;">20% / 33% / 47%</td>
<td style="text-align: right;">87% / 12% / 1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">12.8 (12.413.2); healthy 11.6, infected
1.2</td>
<td style="text-align: right;">11.7 (11.312.1); healthy 11.2, infected
0.5</td>
<td style="text-align: right;">11.8 (11.412.2); healthy 11.3, infected
0.5</td>
<td style="text-align: right;">16.7 (16.217.1); healthy 16.2, infected
0.4</td>
<td style="text-align: right;">19.5 (18.920.1); healthy 18.9, infected
0.6</td>
<td style="text-align: right;">16.4 (15.916.8); healthy 16.0, infected
0.3</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">4.6 (4.64.7) / 14.6 (14.314.9)</td>
<td style="text-align: right;">11.4 (11.211.6) / 15.9 (15.616.2)</td>
<td style="text-align: right;">12.0 (11.912.2) / 15.5 (15.215.8)</td>
<td style="text-align: right;">8.4 (8.38.5) / 17.5 (17.017.9)</td>
<td style="text-align: right;">7.8 (7.77.9) / 20.4 (19.921.0)</td>
<td style="text-align: right;">15.5 (15.315.7) / 17.5 (17.117.8)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">75 (7377)</td>
<td style="text-align: right;">67 (6669)</td>
<td style="text-align: right;">67 (6669)</td>
<td style="text-align: right;">91 (8992)</td>
<td style="text-align: right;">91 (8992)</td>
<td style="text-align: right;">89 (8890)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">9.0 (8.89.2)</td>
<td style="text-align: right;">6.7 (6.56.9)</td>
<td style="text-align: right;">5.9 (5.76.1)</td>
<td style="text-align: right;">5.5 (5.35.6)</td>
<td style="text-align: right;">6.4 (6.26.6)</td>
<td style="text-align: right;">3.6 (3.43.7)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">-1.1 days; 17.6 cases prevented</td>
<td style="text-align: right;">-1.1 days; 21.1 cases prevented</td>
<td style="text-align: right;">+3.8 days; 21.8 cases prevented</td>
<td style="text-align: right;">+6.7 days; 18.8 cases prevented</td>
<td style="text-align: right;">+3.5 days; 23.5 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">9%</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">3%</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">73.3 students, 13 school days each; 10%
infected</td>
<td style="text-align: right;">60.0 students, 17 school days each; 7%
infected</td>
<td style="text-align: right;">53.0 students, 16 school days each; 5%
infected</td>
<td style="text-align: right;">43.5 students, 16 school days each; 5%
infected</td>
<td style="text-align: right;">48.9 students, 16 school days each; 6%
infected</td>
<td style="text-align: right;">25.7 students, 16 school days each; 6%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">18.7 students, 10 school days each; 9%
infected</td>
<td style="text-align: right;">15.6 students, 10 school days each; 9%
infected</td>
<td style="text-align: right;">19.9 students, 5 school days each; 8%
infected</td>
<td style="text-align: right;">13.0 students, 15 school days each; 9%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">68.8 students, 10 school days each; 1%
infected</td>
<td style="text-align: right;">76.7 students, 12 school days each; 2%
infected</td>
<td style="text-align: right;">57.4 students, 16 school days each; 1%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">2.4</td>
<td style="text-align: right;">1.0</td>
<td style="text-align: right;">0.5</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.1</td>
<td style="text-align: right;">0.0</td>
</tr>
</tbody>
</table>

**High school, 90% vaccinated** (473 students, 47 unvaccinated)

<table>
<colgroup>
<col style="width: 22%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Question</th>
<th style="text-align: right;">7/0/0</th>
<th style="text-align: right;">21/0/0</th>
<th style="text-align: right;">21/14/0</th>
<th style="text-align: right;">21/14/7</th>
<th style="text-align: right;">21/7/7</th>
<th style="text-align: right;">21/21/21</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Outbreak risk (more than 10 cases)</td>
<td style="text-align: right;">40%</td>
<td style="text-align: right;">12%</td>
<td style="text-align: right;">8%</td>
<td style="text-align: right;">6%</td>
<td style="text-align: right;">11%</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">1. Students sent home at least once (per
outbreak)</td>
<td style="text-align: right;">27 (2628)</td>
<td style="text-align: right;">23 (2224)</td>
<td style="text-align: right;">27 (2628)</td>
<td style="text-align: right;">42 (4143)</td>
<td style="text-align: right;">42 (4143)</td>
<td style="text-align: right;">42 (4143)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">2. Rounds per outbreak (days a new group
is sent home)</td>
<td style="text-align: right;">3.7 (3.53.8)</td>
<td style="text-align: right;">2.0 (1.92.1)</td>
<td style="text-align: right;">2.0 (1.92.1)</td>
<td style="text-align: right;">2.1 (2.02.2)</td>
<td style="text-align: right;">2.2 (2.12.3)</td>
<td style="text-align: right;">1.0 (1.01.0)</td>
</tr>
<tr class="even">
<td style="text-align: left;">3. Of those sent home: once / twice / 3+
times</td>
<td style="text-align: right;">52% / 25% / 23%</td>
<td style="text-align: right;">87% / 11% / 2%</td>
<td style="text-align: right;">88% / 10% / 2%</td>
<td style="text-align: right;">44% / 40% / 16%</td>
<td style="text-align: right;">39% / 34% / 27%</td>
<td style="text-align: right;">89% / 10% / 1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">4. School days missed per unvaccinated
student (healthy + infected days)</td>
<td style="text-align: right;">5.4 (5.15.7); healthy 5.0, infected
0.4</td>
<td style="text-align: right;">7.3 (6.97.6); healthy 7.0, infected
0.3</td>
<td style="text-align: right;">8.2 (7.98.5); healthy 7.9, infected
0.3</td>
<td style="text-align: right;">12.8 (12.313.2); healthy 12.5, infected
0.3</td>
<td style="text-align: right;">13.3 (12.813.8); healthy 12.9, infected
0.4</td>
<td style="text-align: right;">15.3 (14.915.8); healthy 15.1, infected
0.3</td>
</tr>
<tr class="even">
<td style="text-align: left;">5. School days at home per stay / per
student sent home</td>
<td style="text-align: right;">4.7 (4.74.8) / 8.5 (8.28.7)</td>
<td style="text-align: right;">13.4 (13.213.6) / 14.9 (14.715.1)</td>
<td style="text-align: right;">12.7 (12.512.8) / 14.1 (13.914.3)</td>
<td style="text-align: right;">8.1 (8.08.1) / 14.3 (13.914.7)</td>
<td style="text-align: right;">7.4 (7.37.4) / 14.9 (14.415.3)</td>
<td style="text-align: right;">15.3 (15.215.5) / 17.1 (16.717.5)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">6. Most students at home on one day</td>
<td style="text-align: right;">23 (2224)</td>
<td style="text-align: right;">22 (2123)</td>
<td style="text-align: right;">25 (2426)</td>
<td style="text-align: right;">42 (4143)</td>
<td style="text-align: right;">42 (4143)</td>
<td style="text-align: right;">42 (4143)</td>
</tr>
<tr class="even">
<td style="text-align: left;">7. Weeks from first to last
quarantine</td>
<td style="text-align: right;">5.8 (5.56.1)</td>
<td style="text-align: right;">5.0 (4.85.2)</td>
<td style="text-align: right;">4.5 (4.34.6)</td>
<td style="text-align: right;">4.3 (4.14.4)</td>
<td style="text-align: right;">4.6 (4.44.8)</td>
<td style="text-align: right;">3.3 (3.13.4)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">8. Compared with 7/0/0: extra school days
per unvaccinated student; cases prevented</td>
<td style="text-align: right;">reference</td>
<td style="text-align: right;">+1.9 days; 4.2 cases prevented</td>
<td style="text-align: right;">+2.8 days; 5.1 cases prevented</td>
<td style="text-align: right;">+7.4 days; 5.3 cases prevented</td>
<td style="text-align: right;">+7.9 days; 4.4 cases prevented</td>
<td style="text-align: right;">+10.0 days; 5.7 cases prevented</td>
</tr>
<tr class="even">
<td style="text-align: left;">9. Share of stays where the student was
infected</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">5%</td>
<td style="text-align: right;">4%</td>
<td style="text-align: right;">2%</td>
<td style="text-align: right;">2%</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">24.1 students, 9 school days each; 6%
infected</td>
<td style="text-align: right;">20.4 students, 16 school days each; 6%
infected</td>
<td style="text-align: right;">19.0 students, 15 school days each; 4%
infected</td>
<td style="text-align: right;">16.8 students, 15 school days each; 5%
infected</td>
<td style="text-align: right;">17.5 students, 16 school days each; 6%
infected</td>
<td style="text-align: right;">12.0 students, 16 school days each; 5%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Medium-risk tier: students sent home,
school days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">8.0 students, 10 school days each; 8%
infected</td>
<td style="text-align: right;">7.0 students, 10 school days each; 8%
infected</td>
<td style="text-align: right;">8.5 students, 5 school days each; 7%
infected</td>
<td style="text-align: right;">6.3 students, 15 school days each; 8%
infected</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Low-risk tier: students sent home, school
days each</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">—</td>
<td style="text-align: right;">29.7 students, 9 school days each; 1%
infected</td>
<td style="text-align: right;">32.1 students, 9 school days each; 1%
infected</td>
<td style="text-align: right;">26.9 students, 16 school days each; 1%
infected</td>
</tr>
<tr class="even">
<td style="text-align: left;">Tiers with 0 days: back-to-back 1-day
stays (school days per unvaccinated student)</td>
<td style="text-align: right;">0.5</td>
<td style="text-align: right;">0.3</td>
<td style="text-align: right;">0.2</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.0</td>
<td style="text-align: right;">0.0</td>
</tr>
</tbody>
</table>

</div>

## 1. What each policy costs and what it buys

``` r
p <- ggplot(SUM, aes(u_healthy, risk, colour = pol)) +
  geom_line(aes(group = 1), colour = "grey80", linewidth = 0.5) +   # joined from fewest to most days
  geom_point(size = 3) +
  geom_text(aes(label = pol), vjust = -1, size = 2.8, show.legend = FALSE) +
  facet_grid(vl ~ school, scales = "free_x") +
  scale_colour_manual(values = COLS, name = "Quarantine days (high / medium / low)") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  labs(x = "Healthy school days missed per unvaccinated student, per outbreak", y = sprintf("Chance of more than %d cases", K),
       title = "What each policy costs (healthy days at home) and what it buys (lower outbreak risk)",
       subtitle = "Lower = safer; further right = more school missed by students who were never infected") + th
save10(p, "10_fig1_cost_vs_benefit.png", 12, 10)
```

![](md_figures/10_quarantine_burden_dur-fig1-1.png)

## 2. How many students are at home on each day

``` r
p <- ggplot(HOME, aes(day, home, colour = pol)) + geom_line(linewidth = 0.8) +
  facet_grid(vl ~ school, scales = "free_y") +
  scale_colour_manual(values = COLS, name = "Quarantine days (high / medium / low)") +
  labs(x = "Day since the first case arrived", y = "Students at home (average over runs)",
       title = "How many students are at home on each day",
       subtitle = "Average over all runs, including those where the outbreak died out early") + th
save10(p, "10_fig2_home_by_day.png", 12, 10)
```

![](md_figures/10_quarantine_burden_dur-fig2-1.png)

## 3. Who is at home: healthy or infected

``` r
d3 <- melt(SUM, id.vars = c("school", "vl", "pol"), measure.vars = c("u_inf", "u_healthy"),
           variable.name = "who", value.name = "days")
d3[, who := factor(c(u_inf = "Infected (latent or contagious)", u_healthy = "Healthy")[as.character(who)],
                   c("Healthy", "Infected (latent or contagious)"))]
lab3 <- SUM[, .(school, vl, pol, days_home = u_days_m, lab = sprintf("%.1f (%.0f%% infected)", u_days_m,
                                                          100 * inf_share))]
p <- ggplot(d3, aes(days, pol, fill = who)) + geom_col(width = 0.65) +
  geom_text(data = lab3, aes(days_home, pol, label = lab), hjust = -0.08, size = 2.6, inherit.aes = FALSE) +
  facet_grid(vl ~ school, scales = "free_x") +
  scale_fill_manual(values = c("Healthy" = "#B4B2A9", "Infected (latent or contagious)" = "#C8531E"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.6))) +
  labs(x = "School days missed per unvaccinated student, per outbreak", y = NULL,
       title = "Who is at home: most quarantine days fall on students who were never infected",
       subtitle = "Average over runs; label = share of those days spent by students who were infected") + th
save10(p, "10_fig3_who_is_home.png", 12, 10)
```

![](md_figures/10_quarantine_burden_dur-fig3-1.png)

## 4. How many times a child is sent home

``` r
d4 <- melt(SUM, id.vars = c("school", "vl", "pol"), measure.vars = c("once", "twice", "more"),
           variable.name = "times", value.name = "share")
d4[, times := factor(c(once = "Once", twice = "Twice", more = "3+ times")[as.character(times)],
                     c("Once", "Twice", "3+ times"))]
p <- ggplot(d4, aes(pol, share, fill = times)) + geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
  facet_grid(vl ~ school) +
  scale_fill_manual(values = c("Once" = "#deebf7", "Twice" = "#6baed6", "3+ times" = "#08306b"), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Quarantine days (high / medium / low)", y = "Share of students sent home",
       title = "How many times a student is sent home during one outbreak") + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save10(p, "10_fig4_times_sent_home.png", 12, 10)
```

![](md_figures/10_quarantine_burden_dur-fig4-1.png)

## 5. Burden by risk tier

Counted inside the model (exact): each time a student is sent home, in
the tier they were put in, and every day they spend at home. Tiers with
0 days are not counted. Per outbreak, averaged over runs.

The tier days add up to row 4 of the school cards, apart from one extra:
a “0-day” tier still keeps a student home for 1 day, and when cases are
detected on consecutive days these 1-day stays run into each other.
Those days are shown on their own in the last row of each school card.

``` r
p <- ggplot(TSUM, aes(pol, days_u, fill = tier)) + geom_col(width = 0.7) +
  facet_grid(vl ~ school, scales = "free_y") +
  scale_fill_manual(values = c(High = "#C8531E", Medium = "#E8A33D", Low = "#B4B2A9"), name = "Tier") +
  labs(x = "Quarantine days (high / medium / low)", y = "School days missed per unvaccinated student",
       title = "School days missed by risk tier, per outbreak") + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save10(p, "10_fig5_days_by_tier.png", 12, 10)
```

![](md_figures/10_quarantine_burden_dur-fig5-1.png)

<div class="panel-tabset">

### 60% vaccinated

<table>
<colgroup>
<col style="width: 6%" />
<col style="width: 5%" />
<col style="width: 3%" />
<col style="width: 15%" />
<col style="width: 9%" />
<col style="width: 26%" />
<col style="width: 19%" />
<col style="width: 14%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">School</th>
<th style="text-align: left;">Policy</th>
<th style="text-align: left;">Tier</th>
<th style="text-align: right;">Students sent home (95% CI)</th>
<th style="text-align: right;">Times sent home</th>
<th style="text-align: right;">School days per unvaccinated student (95%
CI)</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">Sent home while infected</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">105 (104107)</td>
<td style="text-align: right;">293.9</td>
<td style="text-align: right;">10.0 (9.810.2)</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">23%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">73 (7075)</td>
<td style="text-align: right;">87.2</td>
<td style="text-align: right;">8.6 (8.38.9)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">16%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">62 (6064)</td>
<td style="text-align: right;">71.3</td>
<td style="text-align: right;">7.2 (6.97.5)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">12%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">22 (2123)</td>
<td style="text-align: right;">23.5</td>
<td style="text-align: right;">1.7 (1.61.8)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">42 (4143)</td>
<td style="text-align: right;">46.4</td>
<td style="text-align: right;">4.8 (4.64.9)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">14%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">16 (1516)</td>
<td style="text-align: right;">16.0</td>
<td style="text-align: right;">1.2 (1.11.2)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">12%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">112 (111113)</td>
<td style="text-align: right;">316.2</td>
<td style="text-align: right;">11.5 (11.211.8)</td>
<td style="text-align: right;">14</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">50 (4851)</td>
<td style="text-align: right;">55.9</td>
<td style="text-align: right;">5.6 (5.45.7)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">17%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">19 (1820)</td>
<td style="text-align: right;">20.2</td>
<td style="text-align: right;">0.7 (0.70.8)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">11%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">120 (119121)</td>
<td style="text-align: right;">399.3</td>
<td style="text-align: right;">14.4 (14.114.8)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">25 (2526)</td>
<td style="text-align: right;">28.0</td>
<td style="text-align: right;">2.9 (2.83.0)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">17%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">14 (1314)</td>
<td style="text-align: right;">13.8</td>
<td style="text-align: right;">1.5 (1.41.5)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">11%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">102 (101103)</td>
<td style="text-align: right;">124.5</td>
<td style="text-align: right;">13.6 (13.313.9)</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">60 (5861)</td>
<td style="text-align: right;">228.9</td>
<td style="text-align: right;">5.9 (5.76.1)</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">14%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">53 (5255)</td>
<td style="text-align: right;">59.9</td>
<td style="text-align: right;">4.7 (4.54.8)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">52 (5153)</td>
<td style="text-align: right;">58.7</td>
<td style="text-align: right;">4.6 (4.44.7)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">1 (11)</td>
<td style="text-align: right;">0.9</td>
<td style="text-align: right;">0.0 (0.00.1)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">51 (4952)</td>
<td style="text-align: right;">57.6</td>
<td style="text-align: right;">4.5 (4.44.6)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">1 (11)</td>
<td style="text-align: right;">0.8</td>
<td style="text-align: right;">0.0 (0.00.0)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">132 (130133)</td>
<td style="text-align: right;">323.3</td>
<td style="text-align: right;">8.7 (8.58.9)</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">52 (5153)</td>
<td style="text-align: right;">59.0</td>
<td style="text-align: right;">4.6 (4.54.8)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">1 (11)</td>
<td style="text-align: right;">0.8</td>
<td style="text-align: right;">0.0 (0.00.0)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">134 (133136)</td>
<td style="text-align: right;">331.2</td>
<td style="text-align: right;">8.9 (8.79.1)</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">48 (4749)</td>
<td style="text-align: right;">54.5</td>
<td style="text-align: right;">4.3 (4.14.4)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">1 (11)</td>
<td style="text-align: right;">0.7</td>
<td style="text-align: right;">0.1 (0.10.1)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">131 (130133)</td>
<td style="text-align: right;">152.9</td>
<td style="text-align: right;">12.3 (12.012.6)</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">179 (178181)</td>
<td style="text-align: right;">688.2</td>
<td style="text-align: right;">17.2 (16.917.4)</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">12%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">161 (158163)</td>
<td style="text-align: right;">223.0</td>
<td style="text-align: right;">16.1 (15.716.4)</td>
<td style="text-align: right;">19</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">143 (140146)</td>
<td style="text-align: right;">180.7</td>
<td style="text-align: right;">13.5 (13.213.9)</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">44 (4245)</td>
<td style="text-align: right;">46.2</td>
<td style="text-align: right;">2.4 (2.32.5)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">105 (102107)</td>
<td style="text-align: right;">122.9</td>
<td style="text-align: right;">9.2 (9.09.5)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">32 (3133)</td>
<td style="text-align: right;">32.8</td>
<td style="text-align: right;">1.7 (1.71.8)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">147 (146149)</td>
<td style="text-align: right;">370.1</td>
<td style="text-align: right;">9.7 (9.49.9)</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">124 (121127)</td>
<td style="text-align: right;">150.9</td>
<td style="text-align: right;">11.1 (10.811.4)</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">42 (4143)</td>
<td style="text-align: right;">45.1</td>
<td style="text-align: right;">1.2 (1.11.2)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">167 (165168)</td>
<td style="text-align: right;">508.9</td>
<td style="text-align: right;">13.2 (12.913.5)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">59 (5760)</td>
<td style="text-align: right;">65.5</td>
<td style="text-align: right;">5.0 (4.95.2)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">28 (2729)</td>
<td style="text-align: right;">28.2</td>
<td style="text-align: right;">2.2 (2.12.2)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">127 (125129)</td>
<td style="text-align: right;">155.5</td>
<td style="text-align: right;">12.2 (11.812.5)</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">1%</td>
</tr>
</tbody>
</table>

### 70% vaccinated

<table>
<colgroup>
<col style="width: 6%" />
<col style="width: 5%" />
<col style="width: 3%" />
<col style="width: 15%" />
<col style="width: 9%" />
<col style="width: 26%" />
<col style="width: 19%" />
<col style="width: 14%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">School</th>
<th style="text-align: left;">Policy</th>
<th style="text-align: left;">Tier</th>
<th style="text-align: right;">Students sent home (95% CI)</th>
<th style="text-align: right;">Times sent home</th>
<th style="text-align: right;">School days per unvaccinated student (95%
CI)</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">Sent home while infected</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">68 (6770)</td>
<td style="text-align: right;">173.1</td>
<td style="text-align: right;">7.9 (7.78.1)</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">22%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">45 (4447)</td>
<td style="text-align: right;">53.0</td>
<td style="text-align: right;">6.9 (6.77.2)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">17%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">37 (3639)</td>
<td style="text-align: right;">41.7</td>
<td style="text-align: right;">5.7 (5.55.9)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">13%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">15 (1416)</td>
<td style="text-align: right;">15.4</td>
<td style="text-align: right;">1.5 (1.41.6)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">28 (2729)</td>
<td style="text-align: right;">31.2</td>
<td style="text-align: right;">4.3 (4.14.4)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">14%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">12 (1112)</td>
<td style="text-align: right;">11.8</td>
<td style="text-align: right;">1.1 (1.11.2)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">83 (8283)</td>
<td style="text-align: right;">213.7</td>
<td style="text-align: right;">10.4 (10.110.7)</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">33 (3234)</td>
<td style="text-align: right;">36.2</td>
<td style="text-align: right;">4.8 (4.75.0)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">17%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">14 (1415)</td>
<td style="text-align: right;">15.0</td>
<td style="text-align: right;">0.7 (0.70.8)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">88 (8789)</td>
<td style="text-align: right;">262.9</td>
<td style="text-align: right;">12.7 (12.413.0)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">19 (1819)</td>
<td style="text-align: right;">20.2</td>
<td style="text-align: right;">2.8 (2.72.9)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">17%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">10 (1010)</td>
<td style="text-align: right;">10.1</td>
<td style="text-align: right;">1.4 (1.41.5)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">11%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">75 (7576)</td>
<td style="text-align: right;">89.9</td>
<td style="text-align: right;">13.1 (12.813.4)</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">40 (3942)</td>
<td style="text-align: right;">133.9</td>
<td style="text-align: right;">4.6 (4.44.8)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">12%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">37 (3638)</td>
<td style="text-align: right;">40.9</td>
<td style="text-align: right;">4.2 (4.14.4)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">37 (3639)</td>
<td style="text-align: right;">41.2</td>
<td style="text-align: right;">4.3 (4.14.4)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">1 (11)</td>
<td style="text-align: right;">0.7</td>
<td style="text-align: right;">0.0 (0.00.1)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">36 (3537)</td>
<td style="text-align: right;">39.7</td>
<td style="text-align: right;">4.1 (4.04.3)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">1 (01)</td>
<td style="text-align: right;">0.6</td>
<td style="text-align: right;">0.0 (0.00.0)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">96 (9598)</td>
<td style="text-align: right;">216.7</td>
<td style="text-align: right;">7.7 (7.57.9)</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">38 (3739)</td>
<td style="text-align: right;">41.5</td>
<td style="text-align: right;">4.3 (4.24.5)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">1 (11)</td>
<td style="text-align: right;">0.6</td>
<td style="text-align: right;">0.0 (0.00.0)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">97 (9598)</td>
<td style="text-align: right;">219.1</td>
<td style="text-align: right;">7.8 (7.68.0)</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">36 (3537)</td>
<td style="text-align: right;">39.9</td>
<td style="text-align: right;">4.2 (4.04.3)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">1 (01)</td>
<td style="text-align: right;">0.5</td>
<td style="text-align: right;">0.1 (0.00.1)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">96 (9497)</td>
<td style="text-align: right;">109.2</td>
<td style="text-align: right;">11.7 (11.412.0)</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">129 (128131)</td>
<td style="text-align: right;">446.4</td>
<td style="text-align: right;">14.8 (14.615.1)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">12%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">110 (107112)</td>
<td style="text-align: right;">147.6</td>
<td style="text-align: right;">14.2 (13.814.6)</td>
<td style="text-align: right;">18</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">94 (9296)</td>
<td style="text-align: right;">113.4</td>
<td style="text-align: right;">11.3 (11.011.7)</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">30 (2931)</td>
<td style="text-align: right;">31.7</td>
<td style="text-align: right;">2.2 (2.12.3)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">73 (7175)</td>
<td style="text-align: right;">84.0</td>
<td style="text-align: right;">8.4 (8.28.7)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">24 (2425)</td>
<td style="text-align: right;">25.3</td>
<td style="text-align: right;">1.7 (1.71.8)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">110 (108111)</td>
<td style="text-align: right;">264.0</td>
<td style="text-align: right;">9.2 (8.99.4)</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">84 (8186)</td>
<td style="text-align: right;">100.6</td>
<td style="text-align: right;">9.8 (9.510.1)</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">30 (3031)</td>
<td style="text-align: right;">32.9</td>
<td style="text-align: right;">1.1 (1.11.2)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">121 (119122)</td>
<td style="text-align: right;">342.2</td>
<td style="text-align: right;">11.8 (11.512.1)</td>
<td style="text-align: right;">14</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">42 (4143)</td>
<td style="text-align: right;">45.3</td>
<td style="text-align: right;">4.6 (4.54.8)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">21 (2122)</td>
<td style="text-align: right;">21.4</td>
<td style="text-align: right;">2.2 (2.12.2)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">91 (8992)</td>
<td style="text-align: right;">106.8</td>
<td style="text-align: right;">11.1 (10.811.4)</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">1%</td>
</tr>
</tbody>
</table>

### 80% vaccinated

<table>
<colgroup>
<col style="width: 6%" />
<col style="width: 5%" />
<col style="width: 3%" />
<col style="width: 15%" />
<col style="width: 9%" />
<col style="width: 26%" />
<col style="width: 19%" />
<col style="width: 14%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">School</th>
<th style="text-align: left;">Policy</th>
<th style="text-align: left;">Tier</th>
<th style="text-align: right;">Students sent home (95% CI)</th>
<th style="text-align: right;">Times sent home</th>
<th style="text-align: right;">School days per unvaccinated student (95%
CI)</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">Sent home while infected</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">34 (3335)</td>
<td style="text-align: right;">75.1</td>
<td style="text-align: right;">5.1 (4.95.3)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">20%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">23 (2223)</td>
<td style="text-align: right;">25.6</td>
<td style="text-align: right;">5.1 (4.95.3)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">16%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">20 (1921)</td>
<td style="text-align: right;">21.6</td>
<td style="text-align: right;">4.5 (4.34.6)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">13%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">9 (89)</td>
<td style="text-align: right;">9.2</td>
<td style="text-align: right;">1.3 (1.31.4)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">17 (1617)</td>
<td style="text-align: right;">18.1</td>
<td style="text-align: right;">3.7 (3.63.9)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">14%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">8 (78)</td>
<td style="text-align: right;">7.8</td>
<td style="text-align: right;">1.1 (1.11.2)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">53 (5254)</td>
<td style="text-align: right;">123.3</td>
<td style="text-align: right;">9.0 (8.79.3)</td>
<td style="text-align: right;">11</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">19 (1819)</td>
<td style="text-align: right;">20.3</td>
<td style="text-align: right;">4.1 (4.04.3)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">16%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">9 (910)</td>
<td style="text-align: right;">9.7</td>
<td style="text-align: right;">0.7 (0.70.7)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">56 (5557)</td>
<td style="text-align: right;">144.8</td>
<td style="text-align: right;">10.5 (10.210.9)</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">12 (1112)</td>
<td style="text-align: right;">12.5</td>
<td style="text-align: right;">2.6 (2.52.7)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">15%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">6 (67)</td>
<td style="text-align: right;">6.6</td>
<td style="text-align: right;">1.4 (1.31.5)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">49 (4850)</td>
<td style="text-align: right;">55.7</td>
<td style="text-align: right;">12.2 (11.912.5)</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">24 (2325)</td>
<td style="text-align: right;">63.2</td>
<td style="text-align: right;">3.3 (3.13.5)</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">23 (2224)</td>
<td style="text-align: right;">25.1</td>
<td style="text-align: right;">3.9 (3.84.1)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">23 (2324)</td>
<td style="text-align: right;">24.7</td>
<td style="text-align: right;">3.9 (3.74.0)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">0 (00)</td>
<td style="text-align: right;">0.4</td>
<td style="text-align: right;">0.0 (0.00.0)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">23 (2324)</td>
<td style="text-align: right;">25.0</td>
<td style="text-align: right;">3.9 (3.84.1)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">0 (00)</td>
<td style="text-align: right;">0.3</td>
<td style="text-align: right;">0.0 (0.00.0)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">61 (6063)</td>
<td style="text-align: right;">122.9</td>
<td style="text-align: right;">6.6 (6.46.8)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">23 (2224)</td>
<td style="text-align: right;">24.9</td>
<td style="text-align: right;">3.9 (3.84.0)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">0 (00)</td>
<td style="text-align: right;">0.4</td>
<td style="text-align: right;">0.0 (0.00.0)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">62 (6063)</td>
<td style="text-align: right;">123.0</td>
<td style="text-align: right;">6.6 (6.46.8)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">23 (2223)</td>
<td style="text-align: right;">24.9</td>
<td style="text-align: right;">3.9 (3.84.0)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">0 (00)</td>
<td style="text-align: right;">0.3</td>
<td style="text-align: right;">0.1 (0.00.1)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">61 (6063)</td>
<td style="text-align: right;">68.1</td>
<td style="text-align: right;">11.0 (10.611.3)</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">73 (7275)</td>
<td style="text-align: right;">208.9</td>
<td style="text-align: right;">10.4 (10.110.7)</td>
<td style="text-align: right;">13</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">60 (5862)</td>
<td style="text-align: right;">73.8</td>
<td style="text-align: right;">10.8 (10.411.1)</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">53 (5155)</td>
<td style="text-align: right;">61.2</td>
<td style="text-align: right;">9.2 (8.99.5)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">19 (1819)</td>
<td style="text-align: right;">19.7</td>
<td style="text-align: right;">2.0 (2.02.1)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">43 (4245)</td>
<td style="text-align: right;">48.2</td>
<td style="text-align: right;">7.3 (7.07.5)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">16 (1516)</td>
<td style="text-align: right;">16.1</td>
<td style="text-align: right;">1.7 (1.61.7)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">69 (6870)</td>
<td style="text-align: right;">147.4</td>
<td style="text-align: right;">7.7 (7.47.9)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">49 (4750)</td>
<td style="text-align: right;">57.1</td>
<td style="text-align: right;">8.4 (8.18.7)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">20 (1921)</td>
<td style="text-align: right;">21.5</td>
<td style="text-align: right;">1.1 (1.11.2)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">77 (7578)</td>
<td style="text-align: right;">191.4</td>
<td style="text-align: right;">9.9 (9.610.2)</td>
<td style="text-align: right;">12</td>
<td style="text-align: right;">2%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">26 (2526)</td>
<td style="text-align: right;">27.4</td>
<td style="text-align: right;">4.2 (4.14.3)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">13 (1313)</td>
<td style="text-align: right;">13.2</td>
<td style="text-align: right;">2.0 (2.02.1)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">57 (5659)</td>
<td style="text-align: right;">64.8</td>
<td style="text-align: right;">10.1 (9.810.4)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">1%</td>
</tr>
</tbody>
</table>

### 90% vaccinated

<table>
<colgroup>
<col style="width: 6%" />
<col style="width: 5%" />
<col style="width: 3%" />
<col style="width: 15%" />
<col style="width: 9%" />
<col style="width: 26%" />
<col style="width: 19%" />
<col style="width: 14%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">School</th>
<th style="text-align: left;">Policy</th>
<th style="text-align: left;">Tier</th>
<th style="text-align: right;">Students sent home (95% CI)</th>
<th style="text-align: right;">Times sent home</th>
<th style="text-align: right;">School days per unvaccinated student (95%
CI)</th>
<th style="text-align: right;">School days per student sent home</th>
<th style="text-align: right;">Sent home while infected</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">9 (89)</td>
<td style="text-align: right;">14.7</td>
<td style="text-align: right;">2.1 (2.02.2)</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">15%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">8 (78)</td>
<td style="text-align: right;">8.4</td>
<td style="text-align: right;">3.4 (3.33.6)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">15%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">7 (77)</td>
<td style="text-align: right;">7.6</td>
<td style="text-align: right;">3.1 (3.03.3)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">14%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">4 (34)</td>
<td style="text-align: right;">3.7</td>
<td style="text-align: right;">1.1 (1.01.1)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">10%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">7 (77)</td>
<td style="text-align: right;">7.2</td>
<td style="text-align: right;">3.0 (2.93.1)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">13%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">4 (34)</td>
<td style="text-align: right;">3.7</td>
<td style="text-align: right;">1.1 (1.01.1)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">9%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">24 (2324)</td>
<td style="text-align: right;">44.7</td>
<td style="text-align: right;">6.6 (6.36.8)</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">7 (67)</td>
<td style="text-align: right;">7.1</td>
<td style="text-align: right;">2.9 (2.83.0)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">14%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">4 (44)</td>
<td style="text-align: right;">3.9</td>
<td style="text-align: right;">0.6 (0.50.6)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">25 (2425)</td>
<td style="text-align: right;">47.4</td>
<td style="text-align: right;">6.9 (6.77.2)</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">5 (56)</td>
<td style="text-align: right;">5.6</td>
<td style="text-align: right;">2.4 (2.32.4)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">15%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">3 (33)</td>
<td style="text-align: right;">3.0</td>
<td style="text-align: right;">1.3 (1.21.3)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">11%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">23 (2223)</td>
<td style="text-align: right;">24.8</td>
<td style="text-align: right;">10.9 (10.611.2)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">11 (1011)</td>
<td style="text-align: right;">18.2</td>
<td style="text-align: right;">1.9 (1.82.0)</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">10 (1011)</td>
<td style="text-align: right;">10.9</td>
<td style="text-align: right;">3.4 (3.23.5)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">11 (1111)</td>
<td style="text-align: right;">11.5</td>
<td style="text-align: right;">3.6 (3.53.7)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">0 (00)</td>
<td style="text-align: right;">0.2</td>
<td style="text-align: right;">0.0 (0.00.0)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">11 (1011)</td>
<td style="text-align: right;">11.2</td>
<td style="text-align: right;">3.5 (3.43.6)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">0 (00)</td>
<td style="text-align: right;">0.2</td>
<td style="text-align: right;">0.0 (0.00.0)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">28 (2729)</td>
<td style="text-align: right;">46.5</td>
<td style="text-align: right;">4.9 (4.75.1)</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">10 (1011)</td>
<td style="text-align: right;">11.0</td>
<td style="text-align: right;">3.4 (3.33.6)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">0 (00)</td>
<td style="text-align: right;">0.2</td>
<td style="text-align: right;">0.0 (0.00.0)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">3%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">28 (2729)</td>
<td style="text-align: right;">47.1</td>
<td style="text-align: right;">5.0 (4.85.2)</td>
<td style="text-align: right;">8</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">10 (1011)</td>
<td style="text-align: right;">11.0</td>
<td style="text-align: right;">3.4 (3.33.6)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">0 (00)</td>
<td style="text-align: right;">0.2</td>
<td style="text-align: right;">0.0 (0.00.1)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">29 (2830)</td>
<td style="text-align: right;">30.8</td>
<td style="text-align: right;">9.8 (9.510.1)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">7/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">24 (2325)</td>
<td style="text-align: right;">47.4</td>
<td style="text-align: right;">4.9 (4.65.1)</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/0/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">20 (2021)</td>
<td style="text-align: right;">23.3</td>
<td style="text-align: right;">7.0 (6.77.3)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">19 (1820)</td>
<td style="text-align: right;">20.4</td>
<td style="text-align: right;">6.3 (6.06.5)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">4%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/0</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">8 (88)</td>
<td style="text-align: right;">8.3</td>
<td style="text-align: right;">1.8 (1.71.8)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">17 (1617)</td>
<td style="text-align: right;">17.9</td>
<td style="text-align: right;">5.5 (5.35.7)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">7 (77)</td>
<td style="text-align: right;">7.2</td>
<td style="text-align: right;">1.5 (1.41.6)</td>
<td style="text-align: right;">10</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/14/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">30 (2930)</td>
<td style="text-align: right;">54.2</td>
<td style="text-align: right;">5.7 (5.55.9)</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">17 (1718)</td>
<td style="text-align: right;">19.2</td>
<td style="text-align: right;">5.8 (5.66.1)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">6%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">9 (89)</td>
<td style="text-align: right;">9.1</td>
<td style="text-align: right;">1.0 (0.91.0)</td>
<td style="text-align: right;">5</td>
<td style="text-align: right;">7%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/7/7</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">32 (3133)</td>
<td style="text-align: right;">61.1</td>
<td style="text-align: right;">6.4 (6.26.7)</td>
<td style="text-align: right;">9</td>
<td style="text-align: right;">1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">High</td>
<td style="text-align: right;">12 (1212)</td>
<td style="text-align: right;">12.7</td>
<td style="text-align: right;">4.0 (3.84.1)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Medium</td>
<td style="text-align: right;">6 (66)</td>
<td style="text-align: right;">6.3</td>
<td style="text-align: right;">2.0 (1.92.0)</td>
<td style="text-align: right;">15</td>
<td style="text-align: right;">8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: left;">21/21/21</td>
<td style="text-align: left;">Low</td>
<td style="text-align: right;">27 (2628)</td>
<td style="text-align: right;">29.8</td>
<td style="text-align: right;">9.4 (9.19.7)</td>
<td style="text-align: right;">16</td>
<td style="text-align: right;">1%</td>
</tr>
</tbody>
</table>

</div>
