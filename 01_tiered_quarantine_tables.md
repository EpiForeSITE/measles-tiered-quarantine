# 1. Tiered quarantine: outbreak probability tables

- [Transmission rate check](#transmission-rate-check)
- [Results by school](#results-by-school)

Tiered-quarantine experiment using the three real per-capita SMART
contact matrices (symmetrised for reciprocity, kept at real magnitude).
R0 = 12 is reached through a per-school transmission rate, not by
scaling the matrix.

For each school we report **P(outbreak e 10 / 20 / 50 cases)** across
vaccination coverage and quarantine strategies, written as **high /
medium / low-risk contact quarantine days**.

All strategies satisfy high e medium e low.

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

| School     | Students | Transmission rate |  R0 |
|:-----------|---------:|------------------:|----:|
| Elementary |      156 |            0.0886 |  12 |
| Middle     |      236 |            0.0455 |  12 |
| High       |      232 |            0.0741 |  12 |

## Results by school

*Unvaccinated students* is the number of students without vaccine at
that coverage (enrolment × (1  coverage)). Only these students can be
quarantined, and an outbreak can only grow beyond this number through
rare vaccine failures.

Differences are in percentage points with 95% confidence intervals. For
example, **+3.1 (+1.2 to +5.0)** means the strategy has 3.1 percentage
points more chance of that outbreak size than the reference.

### Elementary

**Tiered strategies (high/medium/low days)**, compared with 21/21/21 at
the same coverage

*50% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|-----------------------:|-----------------------:|--------------------:|
| 21/21/21 |                    78 |  54.7% |  27.0% |   0.2% |              Reference |              Reference |           Reference |
| 21/14/14 |                    78 |  55.9% |  26.8% |   0.2% |    +1.2 (-1.0 to +3.4) |    -0.2 (-2.2 to +1.7) | +0.0 (-0.2 to +0.3) |
| 21/7/7   |                    78 |  61.2% |  37.8% |   1.2% |    +6.5 (+4.3 to +8.6) |  +10.7 (+8.7 to +12.8) | +1.0 (+0.7 to +1.4) |
| 21/0/0   |                    78 |  69.7% |  56.4% |   5.4% | +15.0 (+12.9 to +17.1) | +29.3 (+27.2 to +31.4) | +5.2 (+4.6 to +6.0) |

*60% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|-----------------------:|-----------------------:|--------------------:|
| 21/21/21 |                    62 |  43.2% |  16.3% |   0.0% |              Reference |              Reference |           Reference |
| 21/14/14 |                    62 |  46.2% |  17.0% |   0.0% |    +3.0 (+0.8 to +5.2) |    +0.7 (-1.0 to +2.3) | -0.0 (-0.1 to +0.1) |
| 21/7/7   |                    62 |  52.1% |  26.3% |   0.0% |   +8.9 (+6.7 to +11.1) |  +10.0 (+8.3 to +11.8) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    62 |  60.8% |  43.1% |   0.4% | +17.5 (+15.4 to +19.7) | +26.9 (+24.9 to +28.7) | +0.3 (+0.1 to +0.6) |

*70% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|-----------------------:|-----------------------:|--------------------:|
| 21/21/21 |                    47 |  29.6% |   5.8% |   0.0% |              Reference |              Reference |           Reference |
| 21/14/14 |                    47 |  30.9% |   6.1% |   0.0% |    +1.3 (-0.7 to +3.3) |    +0.3 (-0.7 to +1.4) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                    47 |  37.8% |  11.8% |   0.0% |   +8.2 (+6.1 to +10.2) |    +6.0 (+4.8 to +7.3) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    47 |  46.1% |  24.1% |   0.0% | +16.5 (+14.4 to +18.6) | +18.4 (+16.9 to +19.9) | +0.0 (-0.1 to +0.1) |

*80% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|-----------------------:|--------------------:|--------------------:|
| 21/21/21 |                    31 |  13.7% |   0.6% |   0.0% |              Reference |           Reference |           Reference |
| 21/14/14 |                    31 |  16.1% |   0.7% |   0.0% |    +2.4 (+0.8 to +4.0) | +0.1 (-0.2 to +0.5) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                    31 |  18.7% |   1.5% |   0.0% |    +5.0 (+3.4 to +6.6) | +0.9 (+0.5 to +1.4) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    31 |  27.9% |   4.8% |   0.0% | +14.2 (+12.4 to +15.9) | +4.2 (+3.5 to +4.9) | +0.0 (-0.1 to +0.1) |

*90% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |                    16 |   1.8% |   0.0% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |                    16 |   2.0% |   0.0% |   0.0% | +0.2 (-0.4 to +0.8) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                    16 |   2.6% |   0.0% |   0.0% | +0.9 (+0.2 to +1.5) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    16 |   5.0% |   0.0% |   0.0% | +3.2 (+2.4 to +4.0) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |

*95% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |                     8 |   0.1% |   0.0% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |                     8 |   0.2% |   0.0% |   0.0% | +0.1 (-0.1 to +0.3) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                     8 |   0.3% |   0.0% |   0.0% | +0.1 (-0.1 to +0.4) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                     8 |   0.6% |   0.0% |   0.0% | +0.4 (+0.2 to +0.7) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |

### Middle

**Tiered strategies (high/medium/low days)**, compared with 21/21/21 at
the same coverage

*50% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |                   118 |  35.6% |   9.4% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |                   118 |  36.5% |  10.5% |   0.0% | +0.9 (-1.2 to +3.0) | +1.1 (-0.2 to +2.4) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                   118 |  36.6% |  10.4% |   0.0% | +0.9 (-1.2 to +3.0) | +0.9 (-0.4 to +2.3) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                   118 |  37.2% |  12.8% |   0.1% | +1.5 (-0.6 to +3.6) | +3.4 (+2.0 to +4.8) | +0.1 (-0.1 to +0.2) |

*60% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |                    94 |  25.1% |   4.3% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |                    94 |  25.8% |   4.4% |   0.0% | +0.8 (-1.2 to +2.7) | +0.1 (-0.8 to +1.0) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                    94 |  27.1% |   3.8% |   0.0% | +2.0 (+0.1 to +3.9) | -0.5 (-1.3 to +0.4) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    94 |  27.5% |   4.5% |   0.0% | +2.4 (+0.5 to +4.4) | +0.2 (-0.7 to +1.1) | +0.0 (-0.1 to +0.1) |

*70% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |                    71 |  15.1% |   0.6% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |                    71 |  14.6% |   0.5% |   0.0% | -0.5 (-2.0 to +1.1) | -0.1 (-0.4 to +0.2) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                    71 |  13.6% |   0.5% |   0.0% | -1.5 (-3.1 to -0.0) | -0.2 (-0.5 to +0.2) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    71 |  14.9% |   0.8% |   0.0% | -0.2 (-1.7 to +1.4) | +0.1 (-0.2 to +0.5) | +0.0 (-0.1 to +0.1) |

*80% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |                    47 |   5.1% |   0.1% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |                    47 |   4.9% |   0.0% |   0.0% | -0.2 (-1.1 to +0.8) | -0.0 (-0.2 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                    47 |   5.2% |   0.0% |   0.0% | +0.1 (-0.8 to +1.1) | -0.0 (-0.2 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    47 |   5.6% |   0.0% |   0.0% | +0.5 (-0.4 to +1.5) | -0.0 (-0.2 to +0.1) | +0.0 (-0.1 to +0.1) |

*90% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |                    24 |   0.5% |   0.0% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |                    24 |   0.3% |   0.0% |   0.0% | -0.2 (-0.5 to +0.1) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                    24 |   0.4% |   0.0% |   0.0% | -0.1 (-0.4 to +0.2) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    24 |   0.4% |   0.0% |   0.0% | -0.0 (-0.3 to +0.3) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |

*95% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |                    12 |   0.1% |   0.0% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |                    12 |   0.0% |   0.0% |   0.0% | -0.1 (-0.2 to +0.1) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                    12 |   0.0% |   0.0% |   0.0% | -0.1 (-0.2 to +0.0) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    12 |   0.0% |   0.0% |   0.0% | -0.1 (-0.2 to +0.0) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |

### High

**Tiered strategies (high/medium/low days)**, compared with 21/21/21 at
the same coverage

*50% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |      Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|-----------------------:|-----------------------:|-----------------------:|
| 21/21/21 |                   116 |  62.4% |  43.5% |   6.0% |              Reference |              Reference |              Reference |
| 21/14/14 |                   116 |  63.9% |  43.8% |   4.2% |    +1.5 (-0.6 to +3.7) |    +0.3 (-1.9 to +2.5) |    -1.8 (-2.8 to -0.8) |
| 21/7/7   |                   116 |  72.6% |  58.6% |   9.8% |  +10.2 (+8.2 to +12.3) | +15.2 (+13.0 to +17.3) |    +3.7 (+2.5 to +4.9) |
| 21/0/0   |                   116 |  79.5% |  72.9% |  22.7% | +17.2 (+15.2 to +19.1) | +29.4 (+27.3 to +31.4) | +16.7 (+15.2 to +18.1) |

*60% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|-----------------------:|-----------------------:|--------------------:|
| 21/21/21 |                    93 |  52.6% |  30.4% |   0.8% |              Reference |              Reference |           Reference |
| 21/14/14 |                    93 |  54.2% |  32.4% |   0.9% |    +1.6 (-0.6 to +3.8) |    +2.0 (-0.0 to +4.0) | +0.1 (-0.3 to +0.5) |
| 21/7/7   |                    93 |  62.8% |  44.3% |   2.1% |  +10.2 (+8.0 to +12.3) | +13.8 (+11.7 to +15.9) | +1.3 (+0.8 to +1.8) |
| 21/0/0   |                    93 |  71.2% |  60.2% |   5.7% | +18.5 (+16.4 to +20.6) | +29.7 (+27.6 to +31.8) | +4.9 (+4.2 to +5.7) |

*70% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|-----------------------:|-----------------------:|--------------------:|
| 21/21/21 |                    70 |  38.2% |  16.6% |   0.1% |              Reference |              Reference |           Reference |
| 21/14/14 |                    70 |  43.8% |  19.4% |   0.1% |    +5.5 (+3.4 to +7.7) |    +2.8 (+1.1 to +4.5) | +0.0 (-0.1 to +0.2) |
| 21/7/7   |                    70 |  51.6% |  28.4% |   0.1% | +13.3 (+11.2 to +15.5) | +11.8 (+10.0 to +13.6) | +0.1 (-0.1 to +0.2) |
| 21/0/0   |                    70 |  61.1% |  43.1% |   0.4% | +22.8 (+20.7 to +25.0) | +26.6 (+24.6 to +28.5) | +0.4 (+0.1 to +0.6) |

*80% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|-----------------------:|-----------------------:|--------------------:|
| 21/21/21 |                    46 |  24.2% |   5.0% |   0.0% |              Reference |              Reference |           Reference |
| 21/14/14 |                    46 |  27.3% |   5.5% |   0.0% |    +3.0 (+1.1 to +4.9) |    +0.6 (-0.4 to +1.6) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                    46 |  32.5% |   8.6% |   0.0% |   +8.3 (+6.3 to +10.2) |    +3.6 (+2.6 to +4.8) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    46 |  42.2% |  16.4% |   0.0% | +17.9 (+15.9 to +19.9) | +11.4 (+10.1 to +12.8) | +0.0 (-0.1 to +0.1) |

*90% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |                    23 |   6.2% |   0.1% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |                    23 |   7.1% |   0.1% |   0.0% | +1.0 (-0.1 to +2.1) | +0.0 (-0.2 to +0.2) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                    23 |   9.7% |   0.3% |   0.0% | +3.5 (+2.3 to +4.7) | +0.1 (-0.1 to +0.4) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    23 |  14.1% |   0.5% |   0.0% | +7.9 (+6.6 to +9.2) | +0.4 (+0.2 to +0.7) | +0.0 (-0.1 to +0.1) |

*95% vaccinated*

| Strategy | Unvaccinated students | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|----------------------:|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |                    12 |   1.1% |   0.0% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |                    12 |   1.3% |   0.0% |   0.0% | +0.2 (-0.2 to +0.7) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |                    12 |   1.8% |   0.0% |   0.0% | +0.7 (+0.2 to +1.3) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |                    12 |   3.2% |   0.0% |   0.0% | +2.1 (+1.5 to +2.8) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
