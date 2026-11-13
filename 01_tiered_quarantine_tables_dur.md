# 1. Tiered quarantine: outbreak probability tables

Tiered-quarantine experiment using the three real per-capita SMART
contact matrices (symmetrised for reciprocity, kept at real magnitude).
R0 = 12 is reached through a per-school transmission rate, not by
scaling the matrix.

For each school we report **P(outbreak ≥ 10 / 20 / 50 cases)** across
vaccination coverage and quarantine strategies, written as **high /
medium / low-risk contact quarantine days**.

All strategies satisfy high ≥ medium ≥ low.

**Difference from the reference.** Each strategy is compared with
*21/21/21* in the same school and at the same vaccination coverage. The
“Difference” columns give the strategy’s outbreak probability minus the
reference’s, in **percentage points**, with a 95% confidence interval
(Newcombe hybrid score method for a difference of two proportions).

- A **positive** difference means more outbreaks than the reference; a
  **negative** one means fewer.
- If the interval **includes 0**, the strategy cannot be told apart from
  the reference with this number of simulations.

Every run of every school × scenario has its own random seed: all seeds
are drawn at once from one starting seed (`SEED` in `measles_model.R`)
without repeats, so the scenarios are independent samples, as the
interval assumes. The interval reflects simulation (Monte Carlo)
uncertainty only, not uncertainty in the model’s assumptions.

<div>

> **Note**
>
> Model code lives in `measles_model.R`. Simulations are cached in
> `results/`; delete the CSV to re-run. For a quick test render with
> `N_SIMS=200 quarto render 01_tiered_quarantine_tables.qmd`.

</div>

## Transmission rate check

Each school gets its own transmission rate so that R0 = 12. The table
checks this against `epiworldR::compute_reproduction_number()`.

``` r
r0_check <- rbindlist(lapply(SCHOOLS, function(s) data.table(
  School              = s,
  Students            = sum(ENROLLMENT[[s]]),
  `Transmission rate` = round(TRANS_RATE[[s]], 4),
  R0 = round(epiworldR::compute_reproduction_number(
         contact_matrix         = SCHOOL_CMAT[[s]],
         transmission_prob      = TRANS_RATE[[s]],
         infectious_period_days = INF_DAYS,
         check_reciprocity      = TRUE)$R, 2)
)))
knitr::kable(r0_check)
```

<table>
<thead>
<tr class="header">
<th style="text-align: left;">School</th>
<th style="text-align: right;">Students</th>
<th style="text-align: right;">Transmission rate</th>
<th style="text-align: right;">R0</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Elementary</td>
<td style="text-align: right;">339</td>
<td style="text-align: right;">0.0887</td>
<td style="text-align: right;">12</td>
</tr>
<tr class="even">
<td style="text-align: left;">Middle</td>
<td style="text-align: right;">466</td>
<td style="text-align: right;">0.0455</td>
<td style="text-align: right;">12</td>
</tr>
<tr class="odd">
<td style="text-align: left;">High</td>
<td style="text-align: right;">473</td>
<td style="text-align: right;">0.0743</td>
<td style="text-align: right;">12</td>
</tr>
</tbody>
</table>

## Results by school

*Unvaccinated students* is the number of students without vaccine at
that coverage (enrolment × (1 − coverage)). Only these students can be
quarantined, and an outbreak can only grow beyond this number through
rare vaccine failures.

Differences are in percentage points with 95% confidence intervals. For
example, **+3.1 (+1.2 to +5.0)** means the strategy has 3.1 percentage
points more chance of that outbreak size than the reference.

``` r
for (s in SCHOOLS) {
  cat("\n\n### ", s, "\n\n")

  cat("**Tiered strategies (high/medium/low days)**, compared with 21/21/21 at the same coverage\n\n")
  ## A separate table per coverage level: a label inside a Markdown table
  ## would end the table, so each coverage level gets its own
  for (v in VACCS_TIERED) {
    cat(sprintf("\n\n*%s vaccinated*\n\n", vacc_lab(v)))
    print(results_table(counts[school == s & set == "tiered" & vacc == v]))
  }
}
```

### Elementary

**Tiered strategies (high/medium/low days)**, compared with 21/21/21 at
the same coverage

*50% vaccinated*

<table style="width:100%;">
<colgroup>
<col style="width: 7%" />
<col style="width: 18%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 19%" />
<col style="width: 19%" />
<col style="width: 19%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">170</td>
<td style="text-align: right;">43.7%</td>
<td style="text-align: right;">20.2%</td>
<td style="text-align: right;">0.7%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">170</td>
<td style="text-align: right;">49.7%</td>
<td style="text-align: right;">23.9%</td>
<td style="text-align: right;">0.7%</td>
<td style="text-align: right;">+6.0 (+3.9 to +8.2)</td>
<td style="text-align: right;">+3.7 (+1.9 to +5.5)</td>
<td style="text-align: right;">+0.0 (-0.4 to +0.4)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">170</td>
<td style="text-align: right;">68.8%</td>
<td style="text-align: right;">47.0%</td>
<td style="text-align: right;">8.2%</td>
<td style="text-align: right;">+25.1 (+23.0 to +27.2)</td>
<td style="text-align: right;">+26.8 (+24.8 to +28.8)</td>
<td style="text-align: right;">+7.4 (+6.6 to +8.4)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">170</td>
<td style="text-align: right;">78.6%</td>
<td style="text-align: right;">64.6%</td>
<td style="text-align: right;">14.4%</td>
<td style="text-align: right;">+34.9 (+32.9 to +36.9)</td>
<td style="text-align: right;">+44.4 (+42.4 to +46.3)</td>
<td style="text-align: right;">+13.7 (+12.6 to +14.8)</td>
</tr>
</tbody>
</table>

*60% vaccinated*

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 18%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 19%" />
<col style="width: 19%" />
<col style="width: 16%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">136</td>
<td style="text-align: right;">34.9%</td>
<td style="text-align: right;">11.5%</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">136</td>
<td style="text-align: right;">40.7%</td>
<td style="text-align: right;">16.2%</td>
<td style="text-align: right;">0.2%</td>
<td style="text-align: right;">+5.8 (+3.7 to +7.9)</td>
<td style="text-align: right;">+4.7 (+3.1 to +6.2)</td>
<td style="text-align: right;">+0.2 (+0.0 to +0.4)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">136</td>
<td style="text-align: right;">60.1%</td>
<td style="text-align: right;">36.9%</td>
<td style="text-align: right;">2.0%</td>
<td style="text-align: right;">+25.2 (+23.1 to +27.3)</td>
<td style="text-align: right;">+25.4 (+23.6 to +27.2)</td>
<td style="text-align: right;">+2.0 (+1.6 to +2.5)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">136</td>
<td style="text-align: right;">70.3%</td>
<td style="text-align: right;">52.9%</td>
<td style="text-align: right;">3.8%</td>
<td style="text-align: right;">+35.4 (+33.3 to +37.4)</td>
<td style="text-align: right;">+41.4 (+39.6 to +43.2)</td>
<td style="text-align: right;">+3.8 (+3.2 to +4.4)</td>
</tr>
</tbody>
</table>

*70% vaccinated*

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 18%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 19%" />
<col style="width: 19%" />
<col style="width: 16%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">102</td>
<td style="text-align: right;">23.5%</td>
<td style="text-align: right;">4.8%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">102</td>
<td style="text-align: right;">29.9%</td>
<td style="text-align: right;">6.9%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+6.4 (+4.4 to +8.3)</td>
<td style="text-align: right;">+2.1 (+1.0 to +3.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">102</td>
<td style="text-align: right;">46.2%</td>
<td style="text-align: right;">21.8%</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">+22.7 (+20.7 to +24.7)</td>
<td style="text-align: right;">+17.0 (+15.5 to +18.4)</td>
<td style="text-align: right;">+0.1 (+0.0 to +0.3)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">102</td>
<td style="text-align: right;">57.0%</td>
<td style="text-align: right;">35.5%</td>
<td style="text-align: right;">0.2%</td>
<td style="text-align: right;">+33.4 (+31.4 to +35.4)</td>
<td style="text-align: right;">+30.7 (+29.1 to +32.3)</td>
<td style="text-align: right;">+0.2 (+0.0 to +0.4)</td>
</tr>
</tbody>
</table>

*80% vaccinated*

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 18%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 19%" />
<col style="width: 18%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">68</td>
<td style="text-align: right;">12.0%</td>
<td style="text-align: right;">0.9%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">68</td>
<td style="text-align: right;">15.5%</td>
<td style="text-align: right;">1.4%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+3.5 (+2.0 to +5.0)</td>
<td style="text-align: right;">+0.5 (+0.0 to +1.0)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">68</td>
<td style="text-align: right;">28.8%</td>
<td style="text-align: right;">6.6%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+16.8 (+15.0 to +18.5)</td>
<td style="text-align: right;">+5.7 (+4.9 to +6.5)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">68</td>
<td style="text-align: right;">36.7%</td>
<td style="text-align: right;">11.9%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+24.7 (+22.9 to +26.5)</td>
<td style="text-align: right;">+10.9 (+9.9 to +12.0)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>

*90% vaccinated*

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 19%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">34</td>
<td style="text-align: right;">2.7%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">34</td>
<td style="text-align: right;">3.2%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+0.5 (-0.2 to +1.2)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">34</td>
<td style="text-align: right;">7.0%</td>
<td style="text-align: right;">0.2%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+4.3 (+3.3 to +5.2)</td>
<td style="text-align: right;">+0.2 (+0.0 to +0.4)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">34</td>
<td style="text-align: right;">8.8%</td>
<td style="text-align: right;">0.4%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+6.1 (+5.1 to +7.1)</td>
<td style="text-align: right;">+0.4 (+0.2 to +0.6)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>

*95% vaccinated*

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 19%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">0.3%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">0.6%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+0.2 (-0.0 to +0.6)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">1.4%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+1.0 (+0.6 to +1.5)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">17</td>
<td style="text-align: right;">1.7%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+1.4 (+1.0 to +1.9)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>

### Middle

**Tiered strategies (high/medium/low days)**, compared with 21/21/21 at
the same coverage

*50% vaccinated*

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 19%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">233</td>
<td style="text-align: right;">24.2%</td>
<td style="text-align: right;">6.7%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">233</td>
<td style="text-align: right;">25.3%</td>
<td style="text-align: right;">7.1%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+1.1 (-0.7 to +3.0)</td>
<td style="text-align: right;">+0.4 (-0.7 to +1.6)</td>
<td style="text-align: right;">-0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">233</td>
<td style="text-align: right;">26.5%</td>
<td style="text-align: right;">7.3%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+2.3 (+0.4 to +4.2)</td>
<td style="text-align: right;">+0.7 (-0.4 to +1.8)</td>
<td style="text-align: right;">-0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">233</td>
<td style="text-align: right;">27.0%</td>
<td style="text-align: right;">8.8%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+2.8 (+0.9 to +4.7)</td>
<td style="text-align: right;">+2.1 (+1.0 to +3.3)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>

*60% vaccinated*

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 19%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">186</td>
<td style="text-align: right;">18.4%</td>
<td style="text-align: right;">3.3%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">186</td>
<td style="text-align: right;">17.5%</td>
<td style="text-align: right;">3.8%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">-0.8 (-2.5 to +0.9)</td>
<td style="text-align: right;">+0.5 (-0.3 to +1.3)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">186</td>
<td style="text-align: right;">18.1%</td>
<td style="text-align: right;">3.6%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">-0.3 (-2.0 to +1.4)</td>
<td style="text-align: right;">+0.3 (-0.5 to +1.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">186</td>
<td style="text-align: right;">18.6%</td>
<td style="text-align: right;">4.4%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+0.3 (-1.4 to +2.0)</td>
<td style="text-align: right;">+1.1 (+0.2 to +1.9)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>

*70% vaccinated*

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 19%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">140</td>
<td style="text-align: right;">10.7%</td>
<td style="text-align: right;">0.9%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">140</td>
<td style="text-align: right;">10.4%</td>
<td style="text-align: right;">1.2%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">-0.3 (-1.6 to +1.1)</td>
<td style="text-align: right;">+0.3 (-0.1 to +0.8)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">140</td>
<td style="text-align: right;">11.1%</td>
<td style="text-align: right;">1.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+0.4 (-1.0 to +1.8)</td>
<td style="text-align: right;">+0.1 (-0.3 to +0.6)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">140</td>
<td style="text-align: right;">10.5%</td>
<td style="text-align: right;">1.2%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">-0.2 (-1.6 to +1.2)</td>
<td style="text-align: right;">+0.4 (-0.1 to +0.9)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>

*80% vaccinated*

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 19%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">93</td>
<td style="text-align: right;">4.2%</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">93</td>
<td style="text-align: right;">4.0%</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">-0.1 (-1.0 to +0.7)</td>
<td style="text-align: right;">-0.1 (-0.2 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">93</td>
<td style="text-align: right;">4.7%</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+0.5 (-0.4 to +1.4)</td>
<td style="text-align: right;">-0.0 (-0.2 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">93</td>
<td style="text-align: right;">5.3%</td>
<td style="text-align: right;">0.2%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+1.1 (+0.2 to +2.1)</td>
<td style="text-align: right;">+0.1 (-0.1 to +0.3)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>

*90% vaccinated*

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 19%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">47</td>
<td style="text-align: right;">0.7%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">47</td>
<td style="text-align: right;">0.6%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">-0.1 (-0.4 to +0.3)</td>
<td style="text-align: right;">-0.0 (-0.1 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">47</td>
<td style="text-align: right;">0.6%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">-0.1 (-0.5 to +0.2)</td>
<td style="text-align: right;">-0.0 (-0.1 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">47</td>
<td style="text-align: right;">0.7%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+0.0 (-0.4 to +0.4)</td>
<td style="text-align: right;">-0.0 (-0.1 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>

*95% vaccinated*

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 19%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">23</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">23</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">-0.1 (-0.2 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">23</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+0.0 (-0.2 to +0.2)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">23</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">-0.0 (-0.2 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>

### High

**Tiered strategies (high/medium/low days)**, compared with 21/21/21 at
the same coverage

*50% vaccinated*

<table style="width:100%;">
<colgroup>
<col style="width: 7%" />
<col style="width: 18%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 19%" />
<col style="width: 19%" />
<col style="width: 19%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">236</td>
<td style="text-align: right;">48.3%</td>
<td style="text-align: right;">26.6%</td>
<td style="text-align: right;">2.5%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">236</td>
<td style="text-align: right;">57.4%</td>
<td style="text-align: right;">31.9%</td>
<td style="text-align: right;">2.9%</td>
<td style="text-align: right;">+9.1 (+6.9 to +11.3)</td>
<td style="text-align: right;">+5.3 (+3.3 to +7.3)</td>
<td style="text-align: right;">+0.4 (-0.3 to +1.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">236</td>
<td style="text-align: right;">82.5%</td>
<td style="text-align: right;">63.2%</td>
<td style="text-align: right;">16.9%</td>
<td style="text-align: right;">+34.2 (+32.2 to +36.1)</td>
<td style="text-align: right;">+36.6 (+34.6 to +38.6)</td>
<td style="text-align: right;">+14.4 (+13.1 to +15.6)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">236</td>
<td style="text-align: right;">91.4%</td>
<td style="text-align: right;">79.8%</td>
<td style="text-align: right;">21.8%</td>
<td style="text-align: right;">+43.1 (+41.3 to +44.9)</td>
<td style="text-align: right;">+53.2 (+51.3 to +55.0)</td>
<td style="text-align: right;">+19.2 (+17.9 to +20.6)</td>
</tr>
</tbody>
</table>

*60% vaccinated*

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 18%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 19%" />
<col style="width: 19%" />
<col style="width: 16%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">189</td>
<td style="text-align: right;">38.6%</td>
<td style="text-align: right;">18.3%</td>
<td style="text-align: right;">0.5%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">189</td>
<td style="text-align: right;">48.4%</td>
<td style="text-align: right;">22.7%</td>
<td style="text-align: right;">0.9%</td>
<td style="text-align: right;">+9.8 (+7.6 to +11.9)</td>
<td style="text-align: right;">+4.5 (+2.7 to +6.2)</td>
<td style="text-align: right;">+0.3 (-0.0 to +0.7)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">189</td>
<td style="text-align: right;">74.7%</td>
<td style="text-align: right;">53.2%</td>
<td style="text-align: right;">7.0%</td>
<td style="text-align: right;">+36.1 (+34.1 to +38.1)</td>
<td style="text-align: right;">+35.0 (+33.0 to +36.9)</td>
<td style="text-align: right;">+6.5 (+5.7 to +7.4)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">189</td>
<td style="text-align: right;">85.1%</td>
<td style="text-align: right;">68.6%</td>
<td style="text-align: right;">9.3%</td>
<td style="text-align: right;">+46.5 (+44.6 to +48.4)</td>
<td style="text-align: right;">+50.3 (+48.4 to +52.1)</td>
<td style="text-align: right;">+8.8 (+7.9 to +9.7)</td>
</tr>
</tbody>
</table>

*70% vaccinated*

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 18%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 19%" />
<col style="width: 19%" />
<col style="width: 16%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">142</td>
<td style="text-align: right;">28.9%</td>
<td style="text-align: right;">9.3%</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">142</td>
<td style="text-align: right;">35.4%</td>
<td style="text-align: right;">12.6%</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">+6.5 (+4.4 to +8.5)</td>
<td style="text-align: right;">+3.3 (+1.9 to +4.6)</td>
<td style="text-align: right;">+0.1 (-0.1 to +0.2)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">142</td>
<td style="text-align: right;">61.9%</td>
<td style="text-align: right;">37.5%</td>
<td style="text-align: right;">1.2%</td>
<td style="text-align: right;">+32.9 (+30.8 to +35.0)</td>
<td style="text-align: right;">+28.2 (+26.5 to +30.0)</td>
<td style="text-align: right;">+1.1 (+0.8 to +1.5)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">142</td>
<td style="text-align: right;">74.3%</td>
<td style="text-align: right;">50.8%</td>
<td style="text-align: right;">1.3%</td>
<td style="text-align: right;">+45.4 (+43.4 to +47.3)</td>
<td style="text-align: right;">+41.5 (+39.7 to +43.3)</td>
<td style="text-align: right;">+1.3 (+0.9 to +1.7)</td>
</tr>
</tbody>
</table>

*80% vaccinated*

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 18%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 5%" />
<col style="width: 19%" />
<col style="width: 19%" />
<col style="width: 16%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">95</td>
<td style="text-align: right;">16.5%</td>
<td style="text-align: right;">2.3%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">95</td>
<td style="text-align: right;">21.2%</td>
<td style="text-align: right;">3.5%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+4.7 (+3.0 to +6.4)</td>
<td style="text-align: right;">+1.2 (+0.5 to +2.0)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">95</td>
<td style="text-align: right;">43.1%</td>
<td style="text-align: right;">15.9%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+26.6 (+24.6 to +28.5)</td>
<td style="text-align: right;">+13.6 (+12.4 to +14.8)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">95</td>
<td style="text-align: right;">52.7%</td>
<td style="text-align: right;">22.1%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+36.2 (+34.2 to +38.1)</td>
<td style="text-align: right;">+19.8 (+18.4 to +21.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>

*90% vaccinated*

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 19%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 20%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">47</td>
<td style="text-align: right;">4.6%</td>
<td style="text-align: right;">0.2%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">47</td>
<td style="text-align: right;">6.2%</td>
<td style="text-align: right;">0.2%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+1.5 (+0.5 to +2.5)</td>
<td style="text-align: right;">+0.0 (-0.2 to +0.2)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">47</td>
<td style="text-align: right;">14.6%</td>
<td style="text-align: right;">1.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+9.9 (+8.7 to +11.2)</td>
<td style="text-align: right;">+0.8 (+0.5 to +1.2)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">47</td>
<td style="text-align: right;">17.9%</td>
<td style="text-align: right;">1.5%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+13.2 (+11.9 to +14.6)</td>
<td style="text-align: right;">+1.3 (+0.9 to +1.7)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>

*95% vaccinated*

<table>
<colgroup>
<col style="width: 8%" />
<col style="width: 19%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 6%" />
<col style="width: 17%" />
<col style="width: 17%" />
<col style="width: 17%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Strategy</th>
<th style="text-align: right;">Unvaccinated students</th>
<th style="text-align: right;">P(≥10)</th>
<th style="text-align: right;">P(≥20)</th>
<th style="text-align: right;">P(≥50)</th>
<th style="text-align: right;">Difference P(≥10)</th>
<th style="text-align: right;">Difference P(≥20)</th>
<th style="text-align: right;">Difference P(≥50)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">21/21/21</td>
<td style="text-align: right;">24</td>
<td style="text-align: right;">0.8%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
<td style="text-align: right;">Reference</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/14/14</td>
<td style="text-align: right;">24</td>
<td style="text-align: right;">1.2%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+0.4 (-0.0 to +0.9)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="odd">
<td style="text-align: left;">21/7/7</td>
<td style="text-align: right;">24</td>
<td style="text-align: right;">2.9%</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+2.1 (+1.5 to +2.7)</td>
<td style="text-align: right;">+0.1 (-0.1 to +0.2)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
<tr class="even">
<td style="text-align: left;">21/0/0</td>
<td style="text-align: right;">24</td>
<td style="text-align: right;">3.6%</td>
<td style="text-align: right;">0.1%</td>
<td style="text-align: right;">0.0%</td>
<td style="text-align: right;">+2.8 (+2.2 to +3.5)</td>
<td style="text-align: right;">+0.1 (+0.0 to +0.3)</td>
<td style="text-align: right;">+0.0 (-0.1 to +0.1)</td>
</tr>
</tbody>
</table>
