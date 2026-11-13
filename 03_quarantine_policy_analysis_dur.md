# 3. Which quarantine policy keeps outbreaks small?

This is the main analysis. Every combination of school, vaccination
coverage and **high / medium / low-risk quarantine days** (with high ≥
medium ≥ low) is simulated once, and all figures are computed from that
single set of runs.

<table>
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Figure</th>
<th style="text-align: left;">Question</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">1</td>
<td style="text-align: left;">How risky is each school with minimal
quarantine?</td>
</tr>
<tr class="even">
<td style="text-align: left;">2</td>
<td style="text-align: left;">What happens to outbreak risk when one
factor changes, and is the change clear?</td>
</tr>
<tr class="odd">
<td style="text-align: left;">3</td>
<td style="text-align: left;">Which quarantine combinations keep risk
low?</td>
</tr>
<tr class="even">
<td style="text-align: left;">4</td>
<td style="text-align: left;">Given high- and medium-risk days, how many
low-risk days are worth it?</td>
</tr>
<tr class="odd">
<td style="text-align: left;">5</td>
<td style="text-align: left;">How big do outbreaks get, with
high-risk-only and with optimal quarantine?</td>
</tr>
<tr class="even">
<td style="text-align: left;">6</td>
<td style="text-align: left;">How does risk change with low-risk days,
across thresholds?</td>
</tr>
</tbody>
</table>

Plain-language conclusions, computed from the results, are at the end
and are also written to `results/conclusions.txt`.

<div>

> **Note**
>
> Model code lives in `measles_model.R`. Simulations are cached in
> `results/`, and PNG copies of every figure are saved to `figures/`.

</div>

## Simulations

``` r
design <- policy_grid(VACCS, HIGHS, MEDS, LOWS)
## Every school x setting gets its own seed, so settings are independent
## samples (needed for the confidence intervals in section 4)
## The cache name includes the vaccination levels, so changing VACCS re-simulates
VACC_TAG <- paste(round(100 * VACCS), collapse = "-")
sims   <- run_or_load(sprintf("policy_grid_indep_seeds_v%s", VACC_TAG),
                      simulate_design(design, independent_seeds = TRUE))

## Runs with an outbreak (> K_MAIN cases, x) out of all runs (n), and the
## probability of an outbreak P = x / n, at every design point
cnt  <- sims[, .(x = sum(size > K_MAIN), n = .N),
             by = .(school, vacc, q_high, q_med, q_low)]
summ <- cnt[, .(school, vacc, q_high, q_med, q_low, P = x / n)]
summ[, vacc_lab := factor(vacc_lab(vacc), levels = VACC_LABS)]
```

136 quarantine/vaccination settings per school × 3 schools × 4,000
simulations = 1,632,000 runs.

<div>

> **What “outbreak risk” means in this report**
>
> Each simulation starts with one infected student. An **outbreak** is a
> simulation that ends with **more than 10 total cases**. The **outbreak
> risk** of a setting is the share of its 4,000 simulations that ended
> this way, i.e. P(final outbreak size \> 10). Figures 1 to 5 and the
> conclusions use this definition; Figure 6 repeats it for several other
> thresholds.

</div>

## 1. The problem

With minimal quarantine (only high-risk contacts quarantined), how
likely is an outbreak?

``` r
d1 <- summ[q_med == 0 & q_low == 0]

p1 <- ggplot(d1, aes(factor(q_high, levels = HIGHS), vacc_lab, fill = P)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = pct(P, accuracy = 1), colour = P > 0.5),
            size = 4, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradient(low = "#fff5eb", high = "#b30000", limits = c(0, 1),
                      labels = pct, name = ylab_risk) +
  facet_wrap(~ school) +
  labs(x = "High-risk quarantine (days)", y = "Students vaccinated",
       title = "1. Outbreak risk with only high-risk contacts quarantined",
       subtitle = "Darker = more dangerous") +
  theme(panel.grid = element_blank())

fig(p1, "fig1_baseline_risk.png", 10, 4.5)
```

![](md_figures/03_quarantine_policy_analysis_dur-fig1-1.png)

## 2. What drives outbreak risk?

Each panel changes **one** factor while the others stay fixed, and shows
the outbreak risk before (grey dot) and after (arrow). The fixed values
are named in each panel’s title. Vaccination is shown step by step (60%
→ 70% → 80% → 90%) at one fixed quarantine policy, with one arrow per
step; each quarantine tier is compared at the lowest and highest
vaccination level. Dots are shaded by vaccination coverage.

Each risk comes from a finite number of simulations, so two settings can
differ a little just by chance. Arrows are therefore coloured by whether
the change is **clear**: blue if the risk is clearly lower, red if
clearly higher, and grey if the 95% confidence interval of the change
includes 0 (the change can’t be told apart from simulation noise). The
table below the figure gives the numbers.

``` r
## The comparisons shown in Figure 2. Every setting used must be in the
## simulated grid (and satisfy high >= medium >= low); otherwise the code
## stops with a message naming the missing setting.
FIG2A <- list(
  vacc   = list(title = "Increasing vaccination",
                steps = VACCS,            # every level, one arrow per step
                fixed = c(q_high = 14L, q_med = 7L, q_low = 0L)),
  q_high = list(title = "Increasing high-risk quarantine",
                from = 14L, to = 21L,     # with medium = 14, high can only be 14 or 21
                fixed = c(q_med = 14L, q_low = 0L), vaccs = range(VACCS)),
  q_med  = list(title = "Increasing medium-risk quarantine",
                from = 0L, to = 14L,
                fixed = c(q_high = 14L, q_low = 0L), vaccs = range(VACCS)),
  q_low  = list(title = "Increasing low-risk quarantine",
                from = 0L, to = 7L,
                fixed = c(q_high = 14L, q_med = 7L), vaccs = range(VACCS)))

TIER_NAMES <- c(q_high = "high", q_med = "medium", q_low = "low")

## Outbreak counts (x out of n runs) of one simulated setting
lookup_xn <- function(s, v, h, m, l) {
  r <- cnt[.(s, v, as.integer(h), as.integer(m), as.integer(l)),
           on = .(school, vacc, q_high, q_med, q_low)]
  if (nrow(r) != 1 || is.na(r$x))
    stop(sprintf("Figure 2: %s, %s vaccinated, %d/%d/%d days was not simulated",
                 s, pct(v), h, m, l))
  list(x = r$x, n = r$n)
}

## "Increasing vaccination: 60% \u2192 70% \u2192 80% \u2192 90% (quarantine 14/7/0 days)"
panel_title <- function(var, sp) {
  change <- if (var == "vacc") paste(pct(sp$steps), collapse = " \u2192 ")
            else sprintf("%d \u2192 %d days", sp$from, sp$to)
  fixed <- if (var == "vacc")
    sprintf("quarantine %s days", paste(sp$fixed, collapse = "/"))
  else paste(sprintf("%s %d", TIER_NAMES[names(sp$fixed)], sp$fixed), collapse = ", ")
  sprintf("%s: %s  (%s%s)", sp$title, change, fixed, if (var == "vacc") "" else " days")
}

## One row per comparison. Vaccination: one row per step (60 -> 70, 70 -> 80,
## ...), all drawn on the school's line. Tiers: one row per school x vaccination.
d2a <- rbindlist(lapply(names(FIG2A), function(var) {
  sp <- FIG2A[[var]]
  q_at <- function(val, v) {
    x <- c(sp$fixed, if (var != "vacc") setNames(val, var))
    list(v = v, h = x[["q_high"]], m = x[["q_med"]], l = x[["q_low"]])
  }
  one <- function(s, a, b, row, label, step, nsteps) {
    ca <- lookup_xn(s, a$v, a$h, a$m, a$l)
    cb <- lookup_xn(s, b$v, b$h, b$m, b$l)
    data.table(factor = var, school = s, row = row, label = label,
               v_before = a$v, v_after = b$v, step = step, nsteps = nsteps,
               x_before = ca$x, n_before = ca$n, x_after = cb$x, n_after = cb$n,
               before = ca$x / ca$n, after = cb$x / cb$n)
  }
  if (var == "vacc") {
    lv <- sp$steps; ns <- length(lv) - 1L
    rbindlist(lapply(SCHOOLS, function(s) rbindlist(lapply(seq_len(ns), function(k)
      one(s, q_at(NA, lv[k]), q_at(NA, lv[k + 1]), row = s,
          label = sprintf("%s \u00b7 %s \u2192 %s", s, pct(lv[k]), pct(lv[k + 1])),
          step = k, nsteps = ns)))))
  } else {
    rbindlist(lapply(SCHOOLS, function(s) rbindlist(lapply(sp$vaccs, function(v) {
      lab <- sprintf("%s \u00b7 %s vaccinated", s, pct(v))
      one(s, q_at(sp$from, v), q_at(sp$to, v), row = lab, label = lab, step = 1L, nsteps = 1L)
    }))))
  }
}))

d2a[, factor_title := factor(vapply(factor, function(v) panel_title(v, FIG2A[[v]]), ""),
                             levels = vapply(names(FIG2A), function(v) panel_title(v, FIG2A[[v]]), ""))]
d2a[, case := factor(paste(factor, row, sep = "|"),
                     levels = rev(unique(paste(factor, row, sep = "|"))))]
## Change = risk after - risk before, with a 95% Newcombe CI
## (settings use separate seeds, so the two risks are independent)
chg <- newcombe_diff(d2a$x_after, d2a$n_after, d2a$x_before, d2a$n_before)
d2a[, `:=`(change = chg$diff, change_lo = chg$lo, change_hi = chg$hi)]
d2a[, direction := factor(fcase(change_hi < 0, "Clearly lower",
                                change_lo > 0, "Clearly higher",
                                default = "No clear change"),
                          levels = c("Clearly lower", "No clear change", "Clearly higher"))]
```

``` r
## Dots shaded by vaccination coverage (light = lowest, dark = highest)
VACC_FILL <- setNames(colorRampPalette(c("#deebf7", "#9ecae1", "#4292c6", "#08306b"))(length(VACCS)), VACC_LABS)
dots2 <- rbind(d2a[, .(factor_title, case, x = before, v = v_before)],
               d2a[factor == "vacc" & step == nsteps, .(factor_title, case, x = after, v = v_after)])
dots2[, vl := factor(vacc_lab(v), levels = VACC_LABS)]
first <- d2a[step == 1]          # risk label at the start of each row
last_ <- d2a[step == nsteps]     # risk label at the end of each row

p2a <- ggplot(d2a, aes(y = case)) +
  geom_segment(aes(x = before, xend = after, yend = case, colour = direction),
               linewidth = 1.2,
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  geom_point(data = dots2, aes(x = x, fill = vl), shape = 21, size = 3.2,
             colour = "grey25", stroke = 0.4) +
  ## labels sit outside the row's ends so they never overlap
  geom_text(data = first, aes(x = before, label = pct(before, accuracy = 1),
                hjust = ifelse(after <= before, -0.4, 1.4)),
            colour = "grey30", size = 3.3) +
  geom_text(data = last_, aes(x = after, label = pct(after, accuracy = 1), colour = direction,
                hjust = ifelse(after <= before, 1.3, -0.3)),
            size = 3.3, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = VACC_FILL, name = "Vaccinated", drop = FALSE) +
  scale_colour_manual(values = c("Clearly lower"   = "#2171b5",
                                 "No clear change" = "grey60",
                                 "Clearly higher"  = "#b30000"),
                      name = "Change in risk (95% CI)", drop = FALSE) +
  scale_x_continuous(labels = pct, limits = c(0, 1)) +
  scale_y_discrete(labels = function(x) sub("^[^|]*\\|", "", x)) +
  facet_wrap(~ factor_title, ncol = 1, scales = "free_y") +
  labs(x = ylab_risk, y = NULL,
       title = "2. What happens when one factor changes?",
       subtitle = paste0("Dot = before (shaded by vaccination), arrow = after. Everything not named in the panel title stays fixed.\n",
                         "Vaccination: one arrow per step. Grey arrow = change within simulation noise (95% CI includes 0)")) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 3.5)),
         colour = guide_legend(order = 2)) +
  theme(panel.grid.major.y = element_blank(),
        strip.text  = element_text(hjust = 0, size = 11),
        axis.text.y = element_text(size = 9.5),
        legend.box  = "vertical")

fig(p2a, "fig2_factor_effects.png", 11, 10)
```

![](md_figures/03_quarantine_policy_analysis_dur-fig2a-1.png)

**Small benefits, with 95% confidence intervals.** Large changes are
clear from the figure, so this table only lists the comparisons where
the **benefit is under 5 percentage points**, where simulation noise
matters most. Benefit = risk before − risk after, in percentage points:
positive = lower risk, negative = higher risk. *Clearly lower / higher*
means the whole interval is above / below 0; *No clear change* means the
interval includes 0. Intervals reflect simulation noise only, not
uncertainty in the model’s assumptions.

``` r
pts_signed <- function(v) sprintf("%+.1f", 100 * v)

## Benefit = before - after (the change with its sign flipped)
d2a[, `:=`(benefit = -change, benefit_lo = -change_hi, benefit_hi = -change_lo)]
small <- d2a[benefit < FIG2_TABLE_BENEFIT]

if (nrow(small) == 0) {
  cat(sprintf("\n\nEvery comparison has a benefit of %s points or more.\n\n",
              round(100 * FIG2_TABLE_BENEFIT)))
} else {
  for (v in names(FIG2A)) {
    d <- small[factor == v]
    cat(sprintf("\n\n*%s*\n\n", panel_title(v, FIG2A[[v]])))
    if (nrow(d) == 0) {
      cat(sprintf("All comparisons in this panel have a benefit of %s points or more.\n\n",
                  round(100 * FIG2_TABLE_BENEFIT)))
      next
    }
    print(knitr::kable(
      d[, .(Setting = label,
            `Risk before` = pct(before, accuracy = 0.1),
            `Risk after`  = pct(after,  accuracy = 0.1),
            `Benefit, points (95% CI)` = sprintf("%s (%s to %s)", pts_signed(benefit),
                                                 pts_signed(benefit_lo), pts_signed(benefit_hi)),
            Verdict = as.character(direction))],
      align = c("l", "r", "r", "r", "l")))
  }
}
```

*Increasing vaccination: 60% → 70% → 80% → 90% (quarantine 14/7/0 days)*

All comparisons in this panel have a benefit of 5 points or more.

*Increasing high-risk quarantine: 14 → 21 days (medium 14, low 0 days)*

<table>
<colgroup>
<col style="width: 31%" />
<col style="width: 13%" />
<col style="width: 12%" />
<col style="width: 27%" />
<col style="width: 15%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Setting</th>
<th style="text-align: right;">Risk before</th>
<th style="text-align: right;">Risk after</th>
<th style="text-align: right;">Benefit, points (95% CI)</th>
<th style="text-align: left;">Verdict</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Elementary · 90% vaccinated</td>
<td style="text-align: right;">6.2%</td>
<td style="text-align: right;">3.9%</td>
<td style="text-align: right;">+2.2 (+1.3 to +3.2)</td>
<td style="text-align: left;">Clearly lower</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle · 90% vaccinated</td>
<td style="text-align: right;">0.8%</td>
<td style="text-align: right;">0.3%</td>
<td style="text-align: right;">+0.5 (+0.2 to +0.8)</td>
<td style="text-align: left;">Clearly lower</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High · 90% vaccinated</td>
<td style="text-align: right;">10.5%</td>
<td style="text-align: right;">7.3%</td>
<td style="text-align: right;">+3.2 (+2.0 to +4.5)</td>
<td style="text-align: left;">Clearly lower</td>
</tr>
</tbody>
</table>

*Increasing medium-risk quarantine: 0 → 14 days (high 14, low 0 days)*

<table>
<colgroup>
<col style="width: 27%" />
<col style="width: 13%" />
<col style="width: 12%" />
<col style="width: 28%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Setting</th>
<th style="text-align: right;">Risk before</th>
<th style="text-align: right;">Risk after</th>
<th style="text-align: right;">Benefit, points (95% CI)</th>
<th style="text-align: left;">Verdict</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Middle · 60% vaccinated</td>
<td style="text-align: right;">30.2%</td>
<td style="text-align: right;">30.0%</td>
<td style="text-align: right;">+0.1 (-1.9 to +2.1)</td>
<td style="text-align: left;">No clear change</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle · 90% vaccinated</td>
<td style="text-align: right;">1.0%</td>
<td style="text-align: right;">0.8%</td>
<td style="text-align: right;">+0.3 (-0.1 to +0.7)</td>
<td style="text-align: left;">No clear change</td>
</tr>
</tbody>
</table>

*Increasing low-risk quarantine: 0 → 7 days (high 14, medium 7 days)*

<table>
<colgroup>
<col style="width: 30%" />
<col style="width: 13%" />
<col style="width: 11%" />
<col style="width: 27%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Setting</th>
<th style="text-align: right;">Risk before</th>
<th style="text-align: right;">Risk after</th>
<th style="text-align: right;">Benefit, points (95% CI)</th>
<th style="text-align: left;">Verdict</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Elementary · 90% vaccinated</td>
<td style="text-align: right;">9.7%</td>
<td style="text-align: right;">7.0%</td>
<td style="text-align: right;">+2.7 (+1.5 to +3.9)</td>
<td style="text-align: left;">Clearly lower</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle · 60% vaccinated</td>
<td style="text-align: right;">31.2%</td>
<td style="text-align: right;">29.3%</td>
<td style="text-align: right;">+1.9 (-0.1 to +3.9)</td>
<td style="text-align: left;">No clear change</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle · 90% vaccinated</td>
<td style="text-align: right;">0.8%</td>
<td style="text-align: right;">0.8%</td>
<td style="text-align: right;">+0.1 (-0.3 to +0.5)</td>
<td style="text-align: left;">No clear change</td>
</tr>
<tr class="even">
<td style="text-align: left;">High · 90% vaccinated</td>
<td style="text-align: right;">19.8%</td>
<td style="text-align: right;">16.2%</td>
<td style="text-align: right;">+3.6 (+1.9 to +5.3)</td>
<td style="text-align: left;">Clearly lower</td>
</tr>
</tbody>
</table>

## 3. Which quarantine combinations keep risk low?

Every quarantine combination, coloured by its outbreak risk on a
continuous scale: **light = low risk, dark = high risk**. The number in
each square is its risk. Squares with a **black outline** have a risk of
10% or less. Columns are medium / low-risk days, grouped by high-risk
days; low-risk days are shown in weekly steps to keep the figure
readable.

``` r
## Shortest policy per target (used in the conclusions)
d3 <- rbindlist(lapply(TARGETS, function(tg) {
  summ[, {
    ok <- .SD[P <= tg][order(q_high + q_med + q_low, q_low, q_med)]
    if (nrow(ok) == 0) .(days = NA_integer_, label = "Not possible")
    else .(days  = ok$q_high[1] + ok$q_med[1] + ok$q_low[1],
           label = sprintf("%d / %d / %d", ok$q_high[1], ok$q_med[1], ok$q_low[1]))
  }, by = .(school, vacc_lab)][, target := tg]
}))

d3p <- summ[q_low %in% FIG3_LOWS]
ML_LEVELS <- unique(d3p[order(q_med, q_low), sprintf("%d/%d", q_med, q_low)])
d3p[, ml := factor(sprintf("%d/%d", q_med, q_low), levels = ML_LEVELS)]
d3p[, high_lab := factor(sprintf("High-risk %d days", q_high),
                         levels = sprintf("High-risk %d days", HIGHS))]

fig3_max <- if (is.na(FIG3_MAX)) max(d3p$P) else FIG3_MAX
d3p[, dark := P > 0.55 * fig3_max]     # white text on dark squares

p3 <- ggplot(d3p, aes(ml, vacc_lab, fill = P)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_gradientn(colours = FIG3_COLORS, limits = c(0, fig3_max),
                       oob = scales::squish, labels = pct,
                       breaks = sort(unique(c(TARGETS[TARGETS < fig3_max],
                                              pretty(c(0, fig3_max), 4)))),
                       name = "Outbreak risk") +
  facet_grid(school ~ high_lab, scales = "free_x", space = "free_x") +
  labs(x = "Medium / low-risk quarantine days", y = "Students vaccinated",
       title = "3. Every quarantine combination: which ones keep risk low?",
       subtitle = paste0("Each square is one policy. Light = low risk, dark = high risk.",
                         if (!is.na(FIG3_HIGHLIGHT))
                           sprintf(" Black outline = risk %s or less.", pct(FIG3_HIGHLIGHT)) else "",
                         "\nRead left to right to find the shortest policy with a low risk")) +
  guides(fill = guide_colourbar(barwidth = unit(8, "cm"))) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 9))

if (FIG3_LABELS)
  p3 <- p3 + geom_text(aes(label = pct(P, accuracy = 1), colour = dark),
                       size = 2.7, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"))

if (!is.na(FIG3_HIGHLIGHT))
  p3 <- p3 + geom_tile(data = d3p[P <= FIG3_HIGHLIGHT], fill = NA,
                       colour = "black", linewidth = 0.8,
                       width = 0.92, height = 0.92)

fig(p3, "fig3_policy_options.png", 11, 7)
```

![](md_figures/03_quarantine_policy_analysis_dur-fig3-1.png)

## 4. How many low-risk quarantine days are worth it?

For each fixed high/medium policy: the fewest low-risk days that give
80% of the benefit low-risk quarantine can provide. Each cell shows
those days and the risk without → with them. When low-risk quarantine
cannot lower the risk by 2% or more, the cell shows **0 days**, the risk
with no low-risk quarantine, and in brackets the largest drop any number
of low-risk days gave.

``` r
d4 <- summ[q_med > 0, {
  o     <- .SD[order(q_low)]
  risk0 <- o$P[o$q_low == 0]
  gain  <- cummax(risk0 - o$P)     # more days never counted as worse (smooths noise)
  G     <- max(gain)
  best  <- if (G < MIN_GAIN - 1e-9) 0L else o$q_low[which(gain >= OPT_SHARE * G)[1]]   # 1e-9: rounding
  ## low-risk days giving the largest benefit (the longest allowed if none helps)
  max_low <- if (G > 0) o$q_low[which.max(gain)] else max(o$q_low)
  .(opt_low  = best,
    max_low  = max_low,
    risk0    = risk0,
    risk_opt = risk0 - gain[o$q_low == best],
    risk_max = risk0 - G)
}, by = .(school, vacc, vacc_lab, q_high, q_med)]

hm_levels <- unique(d4[order(-q_high, -q_med), sprintf("High %d \u00b7 Med %d", q_high, q_med)])
d4[, hm := factor(sprintf("High %d \u00b7 Med %d", q_high, q_med), levels = rev(hm_levels))]

p4 <- ggplot(d4, aes(vacc_lab, hm, fill = opt_low)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(opt_low == 0,
                               sprintf("0 days\n%s (\u2212%s max)",
                                       pct(risk0, accuracy = 1),
                                       pct(risk0 - risk_max, accuracy = 0.1)),
                               sprintf("%d days\n%s \u2192 %s", opt_low,
                                       pct(risk0, accuracy = 1), pct(risk_opt, accuracy = 1))),
                colour = opt_low > 10),
            size = 2.8, lineheight = 0.9, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, max(LOWS)),
                      name = "Optimal low-risk days") +
  facet_wrap(~ school) +
  labs(x = "Students vaccinated", y = "Fixed policy (days)",
       title = "4. How many low-risk quarantine days are worth it?",
       subtitle = sprintf(
         "Each cell: best number of low-risk days, and outbreak risk without \u2192 with them.\n\"0 days\" = low-risk quarantine lowers the risk by less than %s (largest possible drop in brackets)",
         pct(MIN_GAIN))) +
  theme(panel.grid = element_blank())

fig(p4, "fig4_optimal_low.png", 12, 5.8)
```

![](md_figures/03_quarantine_policy_analysis_dur-fig4-1.png)

### How sure is “not worth it”?

A cell is *not worth it* (shown as **0 days** in Figure 4) when the
largest benefit low-risk quarantine gives is under 2 percentage points.
The tables below list only those cells, school by school. For each, the
risk with 0 low-risk days is compared with the risk at the number of
low-risk days that gave the largest benefit (or the longest allowed, if
none helped), with a **95% confidence interval** for the benefit
(Newcombe method for a difference of two proportions; every setting has
its own seed, so the two risks are independent).

- **Benefit** = risk with 0 low-risk days − risk with those days, in
  percentage points (positive = lower risk).
- **Clearly under 2 points?** *Yes* if the whole interval is below 2
  points, so the true benefit is small. *No* if the interval reaches 2
  points or more: the benefit is probably small, but these simulations
  can’t rule out a benefit of that size.

The intervals show simulation uncertainty only, not uncertainty in the
model’s assumptions.

``` r
## Only the "Not worth it" cells
t4 <- d4[opt_low == 0]
t4 <- cnt[q_low == 0, .(school, vacc, q_high, q_med, x0 = x, n0 = n)][
  t4, on = .(school, vacc, q_high, q_med)]
t4 <- cnt[, .(school, vacc, q_high, q_med, max_low = q_low, x1 = x, n1 = n)][
  t4, on = .(school, vacc, q_high, q_med, max_low)]

## Benefit = risk with 0 low-risk days - risk with max_low days
ben <- newcombe_diff(t4$x0, t4$n0, t4$x1, t4$n1)
t4[, `:=`(ben = ben$diff, ben_lo = ben$lo, ben_hi = ben$hi)]

pts <- function(v) sprintf("%+.1f", 100 * v)
t4[, `:=`(
  Vaccinated              = as.character(vacc_lab),
  `High / medium days`    = sprintf("%d / %d", q_high, q_med),
  `Low-risk days compared` = sprintf("0 vs %d", max_low),
  `Risk, 0 days`          = pct(x0 / n0, accuracy = 0.1),
  `Risk, with days`       = pct(x1 / n1, accuracy = 0.1),
  `Benefit, points (95% CI)` = sprintf("%s (%s to %s)", pts(ben), pts(ben_lo), pts(ben_hi)),
  `Clearly under limit?`  = ifelse(ben_hi < MIN_GAIN, "Yes", "No"))]
setorder(t4, school, vacc, -q_high, -q_med)

T4_COLS <- c("Vaccinated", "High / medium days", "Low-risk days compared",
             "Risk, 0 days", "Risk, with days",
             "Benefit, points (95% CI)", "Clearly under limit?")
```

<div class="panel-tabset">

#### Elementary

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 14%" />
<col style="width: 17%" />
<col style="width: 9%" />
<col style="width: 12%" />
<col style="width: 18%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Vaccinated</th>
<th style="text-align: left;">High / medium days</th>
<th style="text-align: left;">Low-risk days compared</th>
<th style="text-align: right;">Risk, 0 days</th>
<th style="text-align: right;">Risk, with days</th>
<th style="text-align: right;">Benefit, points (95% CI)</th>
<th style="text-align: center;">Clearly under 2 points?</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21 / 21</td>
<td style="text-align: left;">0 vs 21</td>
<td style="text-align: right;">2.9%</td>
<td style="text-align: right;">1.3%</td>
<td style="text-align: right;">+1.7 (+1.1 to +2.3)</td>
<td style="text-align: center;">No</td>
</tr>
<tr class="even">
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21 / 14</td>
<td style="text-align: left;">0 vs 14</td>
<td style="text-align: right;">3.9%</td>
<td style="text-align: right;">2.3%</td>
<td style="text-align: right;">+1.6 (+0.9 to +2.4)</td>
<td style="text-align: center;">No</td>
</tr>
<tr class="odd">
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21 / 7</td>
<td style="text-align: left;">0 vs 5</td>
<td style="text-align: right;">6.6%</td>
<td style="text-align: right;">5.8%</td>
<td style="text-align: right;">+0.8 (-0.2 to +1.9)</td>
<td style="text-align: center;">Yes</td>
</tr>
</tbody>
</table>

#### Middle

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 14%" />
<col style="width: 17%" />
<col style="width: 9%" />
<col style="width: 12%" />
<col style="width: 18%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Vaccinated</th>
<th style="text-align: left;">High / medium days</th>
<th style="text-align: left;">Low-risk days compared</th>
<th style="text-align: right;">Risk, 0 days</th>
<th style="text-align: right;">Risk, with days</th>
<th style="text-align: right;">Benefit, points (95% CI)</th>
<th style="text-align: center;">Clearly under 2 points?</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21 / 21</td>
<td style="text-align: left;">0 vs 5</td>
<td style="text-align: right;">15.2%</td>
<td style="text-align: right;">14.2%</td>
<td style="text-align: right;">+1.0 (-0.6 to +2.5)</td>
<td style="text-align: center;">No</td>
</tr>
<tr class="even">
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21 / 14</td>
<td style="text-align: left;">0 vs 14</td>
<td style="text-align: right;">14.2%</td>
<td style="text-align: right;">16.2%</td>
<td style="text-align: right;">-2.0 (-3.5 to -0.4)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="odd">
<td style="text-align: left;">60%</td>
<td style="text-align: left;">21 / 7</td>
<td style="text-align: left;">0 vs 5</td>
<td style="text-align: right;">16.2%</td>
<td style="text-align: right;">16.1%</td>
<td style="text-align: right;">+0.2 (-1.4 to +1.8)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="even">
<td style="text-align: left;">60%</td>
<td style="text-align: left;">14 / 14</td>
<td style="text-align: left;">0 vs 5</td>
<td style="text-align: right;">30.0%</td>
<td style="text-align: right;">28.6%</td>
<td style="text-align: right;">+1.4 (-0.6 to +3.4)</td>
<td style="text-align: center;">No</td>
</tr>
<tr class="odd">
<td style="text-align: left;">60%</td>
<td style="text-align: left;">14 / 7</td>
<td style="text-align: left;">0 vs 7</td>
<td style="text-align: right;">31.2%</td>
<td style="text-align: right;">29.3%</td>
<td style="text-align: right;">+1.9 (-0.1 to +3.9)</td>
<td style="text-align: center;">No</td>
</tr>
<tr class="even">
<td style="text-align: left;">60%</td>
<td style="text-align: left;">7 / 7</td>
<td style="text-align: left;">0 vs 7</td>
<td style="text-align: right;">76.3%</td>
<td style="text-align: right;">74.9%</td>
<td style="text-align: right;">+1.4 (-0.5 to +3.2)</td>
<td style="text-align: center;">No</td>
</tr>
<tr class="odd">
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21 / 21</td>
<td style="text-align: left;">0 vs 3</td>
<td style="text-align: right;">9.2%</td>
<td style="text-align: right;">8.0%</td>
<td style="text-align: right;">+1.2 (-0.0 to +2.5)</td>
<td style="text-align: center;">No</td>
</tr>
<tr class="even">
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21 / 14</td>
<td style="text-align: left;">0 vs 14</td>
<td style="text-align: right;">8.8%</td>
<td style="text-align: right;">8.0%</td>
<td style="text-align: right;">+0.8 (-0.4 to +2.0)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="odd">
<td style="text-align: left;">70%</td>
<td style="text-align: left;">21 / 7</td>
<td style="text-align: left;">0 vs 7</td>
<td style="text-align: right;">9.6%</td>
<td style="text-align: right;">8.0%</td>
<td style="text-align: right;">+1.5 (+0.3 to +2.8)</td>
<td style="text-align: center;">No</td>
</tr>
<tr class="even">
<td style="text-align: left;">70%</td>
<td style="text-align: left;">14 / 7</td>
<td style="text-align: left;">0 vs 7</td>
<td style="text-align: right;">18.8%</td>
<td style="text-align: right;">18.0%</td>
<td style="text-align: right;">+0.8 (-0.9 to +2.5)</td>
<td style="text-align: center;">No</td>
</tr>
<tr class="odd">
<td style="text-align: left;">70%</td>
<td style="text-align: left;">7 / 7</td>
<td style="text-align: left;">0 vs 5</td>
<td style="text-align: right;">56.9%</td>
<td style="text-align: right;">56.7%</td>
<td style="text-align: right;">+0.2 (-2.0 to +2.3)</td>
<td style="text-align: center;">No</td>
</tr>
<tr class="even">
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21 / 21</td>
<td style="text-align: left;">0 vs 21</td>
<td style="text-align: right;">3.3%</td>
<td style="text-align: right;">2.6%</td>
<td style="text-align: right;">+0.7 (-0.0 to +1.5)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="odd">
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21 / 14</td>
<td style="text-align: left;">0 vs 7</td>
<td style="text-align: right;">3.4%</td>
<td style="text-align: right;">3.1%</td>
<td style="text-align: right;">+0.3 (-0.5 to +1.1)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="even">
<td style="text-align: left;">80%</td>
<td style="text-align: left;">21 / 7</td>
<td style="text-align: left;">0 vs 5</td>
<td style="text-align: right;">3.4%</td>
<td style="text-align: right;">3.3%</td>
<td style="text-align: right;">+0.2 (-0.6 to +0.9)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="odd">
<td style="text-align: left;">80%</td>
<td style="text-align: left;">14 / 14</td>
<td style="text-align: left;">0 vs 14</td>
<td style="text-align: right;">7.3%</td>
<td style="text-align: right;">6.3%</td>
<td style="text-align: right;">+1.0 (-0.1 to +2.1)</td>
<td style="text-align: center;">No</td>
</tr>
<tr class="even">
<td style="text-align: left;">80%</td>
<td style="text-align: left;">14 / 7</td>
<td style="text-align: left;">0 vs 7</td>
<td style="text-align: right;">6.8%</td>
<td style="text-align: right;">7.8%</td>
<td style="text-align: right;">-0.9 (-2.1 to +0.2)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="odd">
<td style="text-align: left;">80%</td>
<td style="text-align: left;">7 / 7</td>
<td style="text-align: left;">0 vs 7</td>
<td style="text-align: right;">27.6%</td>
<td style="text-align: right;">28.7%</td>
<td style="text-align: right;">-1.1 (-3.0 to +0.9)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="even">
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21 / 21</td>
<td style="text-align: left;">0 vs 7</td>
<td style="text-align: right;">0.3%</td>
<td style="text-align: right;">0.2%</td>
<td style="text-align: right;">+0.1 (-0.2 to +0.3)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="odd">
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21 / 14</td>
<td style="text-align: left;">0 vs 3</td>
<td style="text-align: right;">0.3%</td>
<td style="text-align: right;">0.2%</td>
<td style="text-align: right;">+0.1 (-0.2 to +0.3)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="even">
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21 / 7</td>
<td style="text-align: left;">0 vs 5</td>
<td style="text-align: right;">0.4%</td>
<td style="text-align: right;">0.2%</td>
<td style="text-align: right;">+0.1 (-0.1 to +0.4)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="odd">
<td style="text-align: left;">90%</td>
<td style="text-align: left;">14 / 14</td>
<td style="text-align: left;">0 vs 10</td>
<td style="text-align: right;">0.8%</td>
<td style="text-align: right;">0.6%</td>
<td style="text-align: right;">+0.1 (-0.2 to +0.5)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="even">
<td style="text-align: left;">90%</td>
<td style="text-align: left;">14 / 7</td>
<td style="text-align: left;">0 vs 5</td>
<td style="text-align: right;">0.8%</td>
<td style="text-align: right;">0.6%</td>
<td style="text-align: right;">+0.2 (-0.2 to +0.6)</td>
<td style="text-align: center;">Yes</td>
</tr>
<tr class="odd">
<td style="text-align: left;">90%</td>
<td style="text-align: left;">7 / 7</td>
<td style="text-align: left;">0 vs 7</td>
<td style="text-align: right;">4.2%</td>
<td style="text-align: right;">4.8%</td>
<td style="text-align: right;">-0.6 (-1.5 to +0.3)</td>
<td style="text-align: center;">Yes</td>
</tr>
</tbody>
</table>

#### High

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 14%" />
<col style="width: 17%" />
<col style="width: 9%" />
<col style="width: 12%" />
<col style="width: 18%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Vaccinated</th>
<th style="text-align: left;">High / medium days</th>
<th style="text-align: left;">Low-risk days compared</th>
<th style="text-align: right;">Risk, 0 days</th>
<th style="text-align: right;">Risk, with days</th>
<th style="text-align: right;">Benefit, points (95% CI)</th>
<th style="text-align: center;">Clearly under 2 points?</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">90%</td>
<td style="text-align: left;">21 / 7</td>
<td style="text-align: left;">0 vs 5</td>
<td style="text-align: right;">13.2%</td>
<td style="text-align: right;">11.6%</td>
<td style="text-align: right;">+1.6 (+0.2 to +3.0)</td>
<td style="text-align: center;">No</td>
</tr>
</tbody>
</table>

</div>

## 5. How big do outbreaks get?

``` r
SIZE_LABS <- c("10 or fewer", "11\u201335", "36\u201350", "More than 50")
SIZE_COLS <- c("#fee5d9", "#fcae91", "#fb6a4a", "#a50f15")

## Share of runs in each size band, for the policy chosen per school/vacc
size_shares <- function(policy) {
  d <- sims[policy, on = .(school, vacc, q_high, q_med, q_low)]
  d[, size_cat := cut(size, breaks = c(-Inf, 10, 35, 50, Inf), labels = SIZE_LABS)]
  d <- d[, .N, by = .(school, vacc, q_high, q_med, q_low, size_cat)][
    , share := N / sum(N), by = .(school, vacc, q_high, q_med, q_low)]
  d[, vacc_lab := factor(vacc_lab(vacc), levels = VACC_LABS)]
  d[]
}

## Optimal policy per school and vaccination (explained in 5b)
opt_pick <- function(dt) {
  dt    <- copy(dt)
  risk0 <- dt[q_high == min(HIGHS) & q_med == 0 & q_low == 0, P]
  best  <- min(dt$P)
  dt[, days := q_high + q_med + q_low]
  need <- if (risk0 - best < MIN_GAIN - 1e-9) 0 else OPT_SHARE * (risk0 - best)
  dt[risk0 - P >= need - 1e-12][order(days, P)][1, .(q_high, q_med, q_low, P,
                                                     P0 = risk0, Pmin = best)]
}

opt_pol <- summ[, opt_pick(.SD), by = .(school, vacc)]
opt_pol[, policy := sprintf("%d/%d/%d", q_high, q_med, q_low)]

## A worked example for the explanation in 5b: the first school and
## vaccination level where quarantine is worth it
ex <- opt_pol[P0 - Pmin >= MIN_GAIN][1]
if (is.na(ex$school)) ex <- opt_pol[1]      # quarantine never worth it anywhere
ex_G      <- ex$P0 - ex$Pmin
ex_target <- ex$P0 - OPT_SHARE * ex_G

size_plot <- function(d, facet = facet_wrap(~ school)) {
  ggplot(d, aes(vacc_lab, share, fill = size_cat)) +
    geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = SIZE_COLS, name = "Total cases", drop = FALSE) +
    facet +
    labs(x = "Students vaccinated", y = "Share of simulations") +
    theme(panel.grid.major.x = element_blank())
}
```

<div class="panel-tabset">

### 5a. Only high-risk contacts quarantined

``` r
## Only high-risk contacts quarantined: 7/0/0, 14/0/0 and 21/0/0
hi_pol <- unique(summ[q_med == 0 & q_low == 0, .(school, vacc, q_high, q_med, q_low)])

d5a <- size_shares(hi_pol)
d5a[, policy_lab := factor(sprintf("High-risk %d days (%d/0/0)", q_high, q_high),
                           levels = sprintf("High-risk %d days (%d/0/0)", HIGHS, HIGHS))]

p5a <- size_plot(d5a, facet_grid(policy_lab ~ school)) +
  scale_y_continuous(labels = pct) +
  labs(title = "5a. How big outbreaks get when only high-risk contacts are quarantined",
       subtitle = "Rows = high-risk quarantine days (medium and low = 0). Light = died out, dark = spread widely")

fig(p5a, "fig5a_outbreak_size_high_only.png", 10, 8.5)
```

![](md_figures/03_quarantine_policy_analysis_dur-fig5a-1.png)

### 5a-1. Only high-risk contacts quarantined, equal school sizes

Same policies as 5a, but **every school has the same number of
students** (240 students, split equally across its grades: Elementary 6
× 40, Middle and High 4 × 60). Everything else is unchanged: the same
contact matrices (contacts per student per day), transmission rates and
vaccination shares. So any difference left between schools comes from
**how students mix**, not from how many students there are. These runs
have their own seeds, independent of the main grid.

``` r
EQUAL_N <- 240L           # students per school (6 x 40, or 4 x 60)
SEED_EQ <- SEED + 51L     # own seeds, independent of the main grid
ENR_EQ  <- sapply(SCHOOLS, function(s) {
  g <- nrow(SCHOOL_CMAT[[s]]); rep(EQUAL_N %/% g, g) }, simplify = FALSE)
design_eq <- policy_grid(VACCS, HIGHS, 0L, 0L)      # 7/0/0, 14/0/0, 21/0/0

## One block of N_SIMS seeds per school x setting, all drawn at once (no repeats)
sim_equal_size <- function() {
  scen <- CJ(school = SCHOOLS, i = seq_len(nrow(design_eq)), sorted = FALSE)
  set.seed(SEED_EQ)
  all_seeds <- sample.int(.Machine$integer.max, nrow(scen) * N_SIMS)
  blocks <- split(all_seeds, rep(seq_len(nrow(scen)), each = N_SIMS))
  ncores <- if (exists("N_THREADS")) N_THREADS else 1L
  rbindlist(lapply(seq_len(nrow(scen)), function(k) {
    s <- scen$school[k]; d <- design_eq[scen$i[k]]; sizes <- ENR_EQ[[s]]
    run_chunk <- function(sd) {
      m <- ModelMeaslesMixingRiskQuarantine(
        n = sum(sizes), prevalence = 1 / sum(sizes),
        contact_matrix = SCHOOL_CMAT[[s]], transmission_rate = TRANS_RATE[[s]],
        prop_vaccinated = d$vacc, detection_rate_quarantine = 0.0,
        contact_tracing_days_window = 7,
        quarantine_period_high   = as.integer(d$q_high),
        quarantine_period_medium = as.integer(d$q_med),
        quarantine_period_low    = as.integer(d$q_low))
      for (g in sizes) m <- add_entity(m, entity("Grade", as.integer(g), FALSE))
      verbose_off(m)
      vapply(sd, function(x) { run(m, ndays = N_DAYS, seed = x); nrow(get_transmissions(m)) }, numeric(1))
    }
    sd <- blocks[[k]]
    nc <- max(1L, min(ncores, length(sd)))
    ch <- split(sd, ceiling(seq_along(sd) / ceiling(length(sd) / nc)))
    message(sprintf("Equal size | %-10s | vacc %3.0f%% | %d/0/0", s, 100 * d$vacc, d$q_high))
    data.table(school = s, d, seed = sd,
               size = unlist(parallel::mclapply(ch, run_chunk, mc.cores = nc), use.names = FALSE))
  }))
}
sims_eq <- run_or_load(sprintf("equal_size%d_high_only_v%s", EQUAL_N, VACC_TAG), sim_equal_size())

## Share of runs in each size band (same bands as 5a)
d5a1 <- copy(sims_eq)[, size_cat := cut(size, breaks = c(-Inf, 10, 35, 50, Inf), labels = SIZE_LABS)]
d5a1 <- d5a1[, .N, by = .(school, vacc, q_high, size_cat)][, share := N / sum(N), by = .(school, vacc, q_high)]
d5a1[, `:=`(school   = factor(school, SCHOOLS),
            vacc_lab = factor(vacc_lab(vacc), levels = VACC_LABS),
            policy_lab = factor(sprintf("High-risk %d days (%d/0/0)", q_high, q_high),
                                levels = sprintf("High-risk %d days (%d/0/0)", HIGHS, HIGHS)))]

p5a1 <- size_plot(d5a1, facet_grid(policy_lab ~ school)) +
  scale_y_continuous(labels = pct) +
  labs(title = sprintf("5a-1. Same as 5a, every school with %d students", EQUAL_N),
       subtitle = "Elementary 6 x 40, Middle 4 x 60, High 4 x 60 students. Same contact matrices and transmission rates as 5a")

fig(p5a1, "fig5a1_outbreak_size_high_only_equal_size.png", 10, 8.5)
```

![](md_figures/03_quarantine_policy_analysis_dur-fig5a1-1.png)

**Outbreak risk (more than 10 cases): real school sizes (5a) vs equal
sizes (5a-1).**

``` r
risk_real <- cnt[q_med == 0 & q_low == 0, .(school, vacc, q_high, real = x / n)]
risk_eq   <- sims_eq[, .(equal = mean(size > K_MAIN)), by = .(school, vacc, q_high)]
tab <- merge(risk_real, risk_eq, by = c("school", "vacc", "q_high"))
tab[, school := factor(school, SCHOOLS)]
setorder(tab, school, q_high, vacc)
tab[, cell := sprintf("%s  %s", pct(real, accuracy = 1), pct(equal, accuracy = 1))]
tab_w <- dcast(tab, school + q_high ~ vacc_lab(vacc), value.var = "cell")
setcolorder(tab_w, c("school", "q_high", VACC_LABS))
setnames(tab_w, c("school", "q_high"), c("School", "High-risk days"))
knitr::kable(tab_w, align = "lr",
             caption = sprintf("Real enrollment (%s) \u2192 equal enrollment (%d students each). Columns = students vaccinated; medium and low = 0 days.",
                               paste(sprintf("%s %d", SCHOOLS, sapply(ENROLLMENT[SCHOOLS], sum)), collapse = ", "), EQUAL_N))
```

<table>
<caption>Real enrollment (Elementary 339, Middle 466, High 473) → equal
enrollment (240 students each). Columns = students vaccinated; medium
and low = 0 days.</caption>
<thead>
<tr class="header">
<th style="text-align: left;">School</th>
<th style="text-align: right;">High-risk days</th>
<th style="text-align: left;">60%</th>
<th style="text-align: right;">70%</th>
<th style="text-align: left;">80%</th>
<th style="text-align: right;">90%</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: right;">7</td>
<td style="text-align: left;">96% 95%</td>
<td style="text-align: right;">91% 88%</td>
<td style="text-align: left;">72% 65%</td>
<td style="text-align: right;">25% 17%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Elementary</td>
<td style="text-align: right;">14</td>
<td style="text-align: left;">79% 77%</td>
<td style="text-align: right;">68% 64%</td>
<td style="text-align: left;">44% 41%</td>
<td style="text-align: right;">11% 8%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: right;">21</td>
<td style="text-align: left;">69% 65%</td>
<td style="text-align: right;">53% 50%</td>
<td style="text-align: left;">33% 29%</td>
<td style="text-align: right;">8% 5%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: right;">7</td>
<td style="text-align: left;">77% 73%</td>
<td style="text-align: right;">58% 48%</td>
<td style="text-align: left;">31% 18%</td>
<td style="text-align: right;">5% 1%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Middle</td>
<td style="text-align: right;">14</td>
<td style="text-align: left;">30% 27%</td>
<td style="text-align: right;">19% 14%</td>
<td style="text-align: left;">8% 4%</td>
<td style="text-align: right;">1% 0%</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: right;">21</td>
<td style="text-align: left;">16% 14%</td>
<td style="text-align: right;">10% 6%</td>
<td style="text-align: left;">4% 2%</td>
<td style="text-align: right;">0% 0%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: right;">7</td>
<td style="text-align: left;">98% 96%</td>
<td style="text-align: right;">94% 92%</td>
<td style="text-align: left;">82% 78%</td>
<td style="text-align: right;">41% 30%</td>
</tr>
<tr class="even">
<td style="text-align: left;">High</td>
<td style="text-align: right;">14</td>
<td style="text-align: left;">91% 88%</td>
<td style="text-align: right;">82% 79%</td>
<td style="text-align: left;">62% 57%</td>
<td style="text-align: right;">20% 16%</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: right;">21</td>
<td style="text-align: left;">84% 80%</td>
<td style="text-align: right;">72% 68%</td>
<td style="text-align: left;">50% 46%</td>
<td style="text-align: right;">14% 11%</td>
</tr>
</tbody>
</table>

Real enrollment (Elementary 339, Middle 466, High 473) → equal
enrollment (240 students each). Columns = students vaccinated; medium
and low = 0 days.

### 5b. Optimal quarantine

**How the optimal policy is computed.** For one school and one
vaccination level, let *P*(*h*,*m*,*l*) be the outbreak risk with *h*
high-, *m* medium- and *l* low-risk quarantine days: the share of
simulations with more than 10 cases. Only policies with *h* ≥ *m* ≥ *l*
are considered.

1.  **Risk with minimal quarantine:** *P*<sub>0</sub> = *P*(7,0,0).
2.  **Lowest risk any tested policy reaches:**
    *P*<sub>min</sub> = min<sub>*h*, *m*, *l*</sub>*P*(*h*,*m*,*l*).
3.  **Largest possible reduction:**
    *G* = *P*<sub>0</sub> − *P*<sub>min</sub>.
4.  **If quarantine barely helps** (*G*\< 2%), the optimal policy is
    minimal quarantine, 7 / 0 / 0.
5.  **Otherwise, a policy is good enough** if it gets at least 80% of
    the largest reduction:
    *P*<sub>0</sub> − *P*(*h*,*m*,*l*) ≥ 0.8 × *G*
6.  **The optimal policy** is the good-enough policy with the fewest
    total quarantine days *h* + *m* + *l*. If several tie, the one with
    the lowest risk is chosen.

*Example:* Elementary school, 60% vaccinated. *P*<sub>0</sub> = 95.8%
and *P*<sub>min</sub> = 30.7%, so *G* = 65.1 percentage points. A
good-enough policy must bring the risk down to 43.8% or lower (95.8% −
0.8 × 65.1 points). The shortest such policy is **21/14/3** (38 days in
total), with a risk of 41.6%.

``` r
d5b     <- size_shares(opt_pol[, .(school, vacc, q_high, q_med, q_low)])
pol_lab <- opt_pol[, .(school, vacc_lab = factor(vacc_lab(vacc), levels = VACC_LABS), policy)]

p5b <- size_plot(d5b) +
  geom_text(data = pol_lab, aes(vacc_lab, 1, label = policy), inherit.aes = FALSE,
            vjust = -0.6, size = 3.3, fontface = "bold") +
  scale_y_continuous(labels = pct, limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  labs(title = "5b. How big outbreaks get with the optimal quarantine",
       subtitle = sprintf(paste0("Above each bar: optimal quarantine days (high/medium/low), the shortest policy\n",
                                 "giving %s of the largest possible risk reduction. Light = died out, dark = spread widely"),
                          pct(OPT_SHARE)))

fig(p5b, "fig5b_outbreak_size_optimal.png", 10, 5)
```

![](md_figures/03_quarantine_policy_analysis_dur-fig5b-1.png)

</div>

## Supplement: every quarantine combination, per school

<div class="panel-tabset">

### Elementary

### Middle

![](md_figures/03_quarantine_policy_analysis_dur-supplement-1.png)

### High

![](md_figures/03_quarantine_policy_analysis_dur-supplement-2.png)![](md_figures/03_quarantine_policy_analysis_dur-supplement-3.png)

</div>

## Conclusions

``` r
out   <- character()
add   <- function(...) out <<- c(out, sprintf(...))
ascii <- function(x) gsub("\u00b7", "-", gsub("\u2192", "->", x))   # plain-text file

add("CONCLUSIONS  (outbreak = more than %d cases)", K_MAIN)
add("")

## 1. Riskiest school
base <- summ[q_high == min(HIGHS) & q_med == 0 & q_low == 0,
             .(P = mean(P)), by = school][order(-P)]
add("1. THE PROBLEM")
add("   With minimal quarantine, %s school is the riskiest (average %s)",
    as.character(base$school[1]), pct(base$P[1], 1))
add("   and %s school the least risky (average %s).",
    as.character(last(base$school)), pct(last(base$P), 1))
add("")

## 2. Levers
add("2. WHAT MATTERS MOST")
for (v in names(FIG2A)) {
  add("   %s:", ascii(panel_title(v, FIG2A[[v]])))
  for (k in seq_len(nrow(d2a[factor == v]))) {
    r <- d2a[factor == v][k]
    add("     %s: risk %s -> %s (change %s points, 95%% CI %s to %s: %s)",
        ascii(r$label), pct(r$before, 1), pct(r$after, 1),
        pts_signed(r$change), pts_signed(r$change_lo), pts_signed(r$change_hi),
        tolower(r$direction))
  }
}
add("")

## 3. Policy
add("3. WHAT A SCHOOL SHOULD DO")
for (tg in TARGETS) {
  add("   To keep risk under %s:", pct(tg))
  for (s in SCHOOLS) {
    t  <- d3[school == s & target == tg][order(vacc_lab)]
    ok <- t[!is.na(days)]
    mn <- t[label == MIN_POLICY]
    add("     %s: %s; minimal quarantine is enough %s.", s,
        if (nrow(ok)) sprintf("achievable from %s vaccination", as.character(ok$vacc_lab[1]))
        else "not achievable at any tested vaccination level",
        if (nrow(mn)) sprintf("from %s vaccination", as.character(mn$vacc_lab[1]))
        else "at no tested level")
  }
}
add("")

## 4. Optimal low-risk days
add("4. LOW-RISK QUARANTINE")
for (s in SCHOOLS) {
  o <- d4[school == s]
  add("   %s: the best low-risk quarantine is usually about %g days.", s, median(o$opt_low))
  add("     It is not worth it (risk drops less than %s) in %s of high/medium settings.",
      pct(MIN_GAIN), pct(mean(o$opt_low == 0), 1))
  add("     On average, risk goes %s -> %s with the best number of days,",
      pct(mean(o$risk0), 1), pct(mean(o$risk_opt), 1))
  add("     and only to %s even with the longest low-risk quarantine.", pct(mean(o$risk_max), 1))
}
add("")

## 5. Size
add("5. SIZE OF OUTBREAKS")
for (s in SCHOOLS) {
  b <- sims[school == s & q_high == min(HIGHS) & q_med == 0 & q_low == 0,
            .(share = mean(size > 50)), by = vacc][order(vacc)]
  add("   %s, minimal quarantine: chance of more than 50 cases falls from %s at %s to %s at %s vaccination.",
      s, pct(b$share[1], 1), pct(b$vacc[1]), pct(last(b$share), 1), pct(last(b$vacc)))

  op <- opt_pol[school == s][order(vacc)]
  bo <- sims[op, on = .(school, vacc, q_high, q_med, q_low)][
    , .(share = mean(size > 50)), by = vacc][order(vacc)]
  add("     Optimal quarantine (high/medium/low days): %s.",
      paste(sprintf("%s at %s", op$policy, pct(op$vacc)), collapse = "; "))
  add("     With it, the chance of more than 50 cases is %s at %s and %s at %s vaccination.",
      pct(bo$share[1], 1), pct(bo$vacc[1]), pct(last(bo$share), 1), pct(last(bo$vacc)))
}

writeLines(out, file.path(RESULTS_DIR, "conclusions.txt"))
```

``` r
cat("```text", out, "```", sep = "\n")
```

``` text
CONCLUSIONS  (outbreak = more than 10 cases)

1. THE PROBLEM
   With minimal quarantine, High school is the riskiest (average 79%)
   and Middle school the least risky (average 42%).

2. WHAT MATTERS MOST
   Increasing vaccination: 60% -> 70% -> 80% -> 90%  (quarantine 14/7/0 days):
     Elementary - 60% -> 70%: risk 80% -> 67% (change -12.4 points, 95% CI -14.4 to -10.5: clearly lower)
     Elementary - 70% -> 80%: risk 67% -> 45% (change -22.8 points, 95% CI -24.9 to -20.7: clearly lower)
     Elementary - 80% -> 90%: risk 45% -> 10% (change -34.9 points, 95% CI -36.7 to -33.1: clearly lower)
     Middle - 60% -> 70%: risk 31% -> 19% (change -12.5 points, 95% CI -14.3 to -10.6: clearly lower)
     Middle - 70% -> 80%: risk 19% -> 7% (change -11.9 points, 95% CI -13.4 to -10.5: clearly lower)
     Middle - 80% -> 90%: risk 7% -> 1% (change -6.0 points, 95% CI -6.9 to -5.2: clearly lower)
     High - 60% -> 70%: risk 90% -> 82% (change -8.0 points, 95% CI -9.5 to -6.5: clearly lower)
     High - 70% -> 80%: risk 82% -> 62% (change -20.5 points, 95% CI -22.4 to -18.5: clearly lower)
     High - 80% -> 90%: risk 62% -> 20% (change -41.7 points, 95% CI -43.6 to -39.7: clearly lower)
   Increasing high-risk quarantine: 14 -> 21 days  (medium 14, low 0 days):
     Elementary - 60% vaccinated: risk 62% -> 45% (change -17.0 points, 95% CI -19.2 to -14.9: clearly lower)
     Elementary - 90% vaccinated: risk 6% -> 4% (change -2.2 points, 95% CI -3.2 to -1.3: clearly lower)
     Middle - 60% vaccinated: risk 30% -> 14% (change -15.8 points, 95% CI -17.6 to -14.0: clearly lower)
     Middle - 90% vaccinated: risk 1% -> 0% (change -0.5 points, 95% CI -0.8 to -0.2: clearly lower)
     High - 60% vaccinated: risk 69% -> 54% (change -15.0 points, 95% CI -17.1 to -12.9: clearly lower)
     High - 90% vaccinated: risk 11% -> 7% (change -3.2 points, 95% CI -4.5 to -2.0: clearly lower)
   Increasing medium-risk quarantine: 0 -> 14 days  (high 14, low 0 days):
     Elementary - 60% vaccinated: risk 79% -> 62% (change -16.4 points, 95% CI -18.4 to -14.4: clearly lower)
     Elementary - 90% vaccinated: risk 11% -> 6% (change -5.3 points, 95% CI -6.6 to -4.1: clearly lower)
     Middle - 60% vaccinated: risk 30% -> 30% (change -0.1 points, 95% CI -2.1 to +1.9: no clear change)
     Middle - 90% vaccinated: risk 1% -> 1% (change -0.3 points, 95% CI -0.7 to +0.1: no clear change)
     High - 60% vaccinated: risk 91% -> 69% (change -22.5 points, 95% CI -24.1 to -20.8: clearly lower)
     High - 90% vaccinated: risk 20% -> 11% (change -9.8 points, 95% CI -11.4 to -8.2: clearly lower)
   Increasing low-risk quarantine: 0 -> 7 days  (high 14, medium 7 days):
     Elementary - 60% vaccinated: risk 80% -> 67% (change -12.4 points, 95% CI -14.3 to -10.5: clearly lower)
     Elementary - 90% vaccinated: risk 10% -> 7% (change -2.7 points, 95% CI -3.9 to -1.5: clearly lower)
     Middle - 60% vaccinated: risk 31% -> 29% (change -1.9 points, 95% CI -3.9 to +0.1: no clear change)
     Middle - 90% vaccinated: risk 1% -> 1% (change -0.1 points, 95% CI -0.5 to +0.3: no clear change)
     High - 60% vaccinated: risk 90% -> 81% (change -9.2 points, 95% CI -10.7 to -7.7: clearly lower)
     High - 90% vaccinated: risk 20% -> 16% (change -3.6 points, 95% CI -5.3 to -1.9: clearly lower)

3. WHAT A SCHOOL SHOULD DO
   To keep risk under 10%:
     Elementary: achievable from 90% vaccination; minimal quarantine is enough at no tested level.
     Middle: achievable from 70% vaccination; minimal quarantine is enough from 90% vaccination.
     High: achievable from 90% vaccination; minimal quarantine is enough at no tested level.
   To keep risk under 20%:
     Elementary: achievable from 80% vaccination; minimal quarantine is enough at no tested level.
     Middle: achievable from 60% vaccination; minimal quarantine is enough from 90% vaccination.
     High: achievable from 80% vaccination; minimal quarantine is enough at no tested level.
   To keep risk under 30%:
     Elementary: achievable from 70% vaccination; minimal quarantine is enough from 90% vaccination.
     Middle: achievable from 60% vaccination; minimal quarantine is enough from 90% vaccination.
     High: achievable from 70% vaccination; minimal quarantine is enough at no tested level.

4. LOW-RISK QUARANTINE
   Elementary: the best low-risk quarantine is usually about 7 days.
     It is not worth it (risk drops less than 2%) in 12% of high/medium settings.
     On average, risk goes 40% -> 34% with the best number of days,
     and only to 34% even with the longest low-risk quarantine.
   Middle: the best low-risk quarantine is usually about 0 days.
     It is not worth it (risk drops less than 2%) in 96% of high/medium settings.
     On average, risk goes 15% -> 15% with the best number of days,
     and only to 15% even with the longest low-risk quarantine.
   High: the best low-risk quarantine is usually about 8.5 days.
     It is not worth it (risk drops less than 2%) in 4% of high/medium settings.
     On average, risk goes 49% -> 42% with the best number of days,
     and only to 42% even with the longest low-risk quarantine.

5. SIZE OF OUTBREAKS
   Elementary, minimal quarantine: chance of more than 50 cases falls from 72% at 60% to 0% at 90% vaccination.
     Optimal quarantine (high/medium/low days): 21/14/3 at 60%; 21/14/0 at 70%; 21/14/0 at 80%; 14/14/3 at 90%.
     With it, the chance of more than 50 cases is 1% at 60% and 0% at 90% vaccination.
   Middle, minimal quarantine: chance of more than 50 cases falls from 1% at 60% to 0% at 90% vaccination.
     Optimal quarantine (high/medium/low days): 21/0/0 at 60%; 21/0/0 at 70%; 14/0/0 at 80%; 14/0/0 at 90%.
     With it, the chance of more than 50 cases is 0% at 60% and 0% at 90% vaccination.
   High, minimal quarantine: chance of more than 50 cases falls from 85% at 60% to 0% at 90% vaccination.
     Optimal quarantine (high/medium/low days): 21/21/0 at 60%; 21/21/0 at 70%; 21/14/0 at 80%; 14/14/0 at 90%.
     With it, the chance of more than 50 cases is 1% at 60% and 0% at 90% vaccination.
```
