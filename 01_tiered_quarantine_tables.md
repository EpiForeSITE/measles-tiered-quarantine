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

The diagnostic single-tier scenarios in the first table do not satisfy
high e medium e low; that constraint is applied to every other analysis.

**Difference from the reference.** Each strategy is compared with a
reference scenario in the same school and at the same vaccination
coverage: *No quarantine* for the no-vaccination table, and *21/21/21*
for the tiered tables. The “Difference” columns give the strategy’s
outbreak probability minus the reference’s, in **percentage points**,
with a 95% confidence interval (Newcombe hybrid score method for a
difference of two proportions).

- A **positive** difference means more outbreaks than the reference; a
  **negative** one means fewer.
- If the interval **includes 0**, the strategy cannot be told apart from
  the reference with this number of simulations.

Every school × scenario is simulated with its own random seed, so the
scenarios are independent samples, as the interval assumes. The interval
reflects simulation (Monte Carlo) uncertainty only, not uncertainty in
the model’s assumptions.

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

Differences are in percentage points with 95% confidence intervals. For
example, **+3.1 (+1.2 to +5.0)** means the strategy has 3.1 percentage
points more chance of that outbreak size than the reference.

### Elementary

**No vaccination: which single tier to quarantine** (compared with no
quarantine)

| Scenario             | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |      Difference P(≥20) |      Difference P(≥50) |
|:---------------------|-------:|-------:|-------:|--------------------:|-----------------------:|-----------------------:|
| No quarantine        |  97.4% |  97.4% |  97.2% |           Reference |              Reference |              Reference |
| Only high-risk (21d) |  88.5% |  85.3% |  75.3% | -8.8 (-9.9 to -7.7) | -12.1 (-13.3 to -10.9) | -21.9 (-23.3 to -20.5) |
| Only med-risk (21d)  |  96.0% |  96.0% |  94.5% | -1.3 (-2.1 to -0.5) |    -1.3 (-2.1 to -0.6) |    -2.7 (-3.6 to -1.8) |
| Only low-risk (21d)  |  96.6% |  93.5% |  57.0% | -0.8 (-1.5 to -0.0) |    -3.9 (-4.8 to -3.0) | -40.3 (-41.9 to -38.6) |

**Tiered strategies (high/medium/low days)**, compared with 21/21/21 at
the same coverage

*50% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|-----------------------:|-----------------------:|--------------------:|
| 21/21/21 |  54.2% |  32.3% |  19.5% |              Reference |              Reference |           Reference |
| 21/14/14 |  57.6% |  35.3% |  19.7% |    +3.3 (+1.2 to +5.5) |    +3.0 (+0.9 to +5.0) | +0.2 (-1.6 to +1.9) |
| 21/7/7   |  61.7% |  44.0% |  20.5% |    +7.5 (+5.3 to +9.6) |  +11.7 (+9.6 to +13.8) | +1.0 (-0.7 to +2.8) |
| 21/0/0   |  71.0% |  59.4% |  23.1% | +16.8 (+14.7 to +18.9) | +27.1 (+24.9 to +29.1) | +3.6 (+1.8 to +5.4) |

*60% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|-----------------------:|-----------------------:|--------------------:|
| 21/21/21 |  45.8% |  24.4% |  15.7% |              Reference |              Reference |           Reference |
| 21/14/14 |  48.4% |  26.2% |  15.2% |    +2.6 (+0.4 to +4.8) |    +1.8 (-0.1 to +3.7) | -0.5 (-2.1 to +1.1) |
| 21/7/7   |  53.2% |  33.0% |  15.1% |    +7.4 (+5.2 to +9.6) |   +8.6 (+6.7 to +10.6) | -0.6 (-2.1 to +1.0) |
| 21/0/0   |  59.8% |  46.0% |   9.7% | +14.0 (+11.8 to +16.2) | +21.6 (+19.5 to +23.6) | -6.0 (-7.5 to -4.6) |

*70% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|-----------------------:|-----------------------:|--------------------:|
| 21/21/21 |  31.1% |  13.9% |   2.7% |              Reference |              Reference |           Reference |
| 21/14/14 |  33.0% |  13.2% |   2.0% |    +1.9 (-0.1 to +3.9) |    -0.8 (-2.3 to +0.8) | -0.7 (-1.4 to -0.1) |
| 21/7/7   |  39.8% |  18.7% |   2.5% |   +8.8 (+6.7 to +10.9) |    +4.7 (+3.1 to +6.4) | -0.2 (-1.0 to +0.4) |
| 21/0/0   |  47.5% |  28.6% |   1.2% | +16.4 (+14.3 to +18.5) | +14.6 (+12.9 to +16.4) | -1.5 (-2.1 to -0.9) |

*80% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|-----------------------:|--------------------:|--------------------:|
| 21/21/21 |  15.3% |   6.2% |   0.0% |              Reference |           Reference |           Reference |
| 21/14/14 |  17.2% |   6.2% |   0.0% |    +1.9 (+0.3 to +3.5) | -0.1 (-1.2 to +1.0) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |  21.7% |   6.4% |   0.0% |    +6.4 (+4.7 to +8.1) | +0.2 (-0.9 to +1.2) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |  28.5% |   8.6% |   0.0% | +13.2 (+11.4 to +15.0) | +2.4 (+1.2 to +3.5) | +0.0 (-0.1 to +0.1) |

*90% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |   3.8% |   0.4% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |   3.9% |   0.4% |   0.0% | +0.1 (-0.7 to +1.0) | +0.0 (-0.3 to +0.3) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |   3.6% |   0.2% |   0.0% | -0.3 (-1.1 to +0.6) | -0.2 (-0.5 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |   5.9% |   0.2% |   0.0% | +2.0 (+1.1 to +3.0) | -0.2 (-0.5 to +0.0) | +0.0 (-0.1 to +0.1) |

*95% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |   1.2% |   0.0% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |   0.9% |   0.0% |   0.0% | -0.4 (-0.8 to +0.1) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |   1.0% |   0.0% |   0.0% | -0.2 (-0.7 to +0.3) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |   0.8% |   0.0% |   0.0% | -0.5 (-0.9 to -0.0) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |

### Middle

**No vaccination: which single tier to quarantine** (compared with no
quarantine)

| Scenario             | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |      Difference P(≥50) |
|:---------------------|-------:|-------:|-------:|-----------------------:|-----------------------:|-----------------------:|
| No quarantine        |  92.5% |  92.3% |  91.3% |              Reference |              Reference |              Reference |
| Only high-risk (21d) |  68.5% |  54.0% |  32.7% | -23.9 (-25.6 to -22.3) | -38.2 (-40.0 to -36.5) | -58.6 (-60.2 to -56.8) |
| Only med-risk (21d)  |  92.8% |  92.8% |  91.8% |    +0.4 (-0.7 to +1.5) |    +0.6 (-0.6 to +1.7) |    +0.5 (-0.7 to +1.7) |
| Only low-risk (21d)  |  92.1% |  91.8% |  90.5% |    -0.4 (-1.5 to +0.8) |    -0.5 (-1.7 to +0.7) |    -0.8 (-2.1 to +0.4) |

**Tiered strategies (high/medium/low days)**, compared with 21/21/21 at
the same coverage

*50% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |    Difference P(≥50) |
|:---------|-------:|-------:|-------:|--------------------:|--------------------:|---------------------:|
| 21/21/21 |  37.1% |  18.1% |  12.7% |           Reference |           Reference |            Reference |
| 21/14/14 |  37.8% |  20.2% |  11.2% | +0.7 (-1.4 to +2.8) | +2.1 (+0.4 to +3.8) |  -1.4 (-2.9 to -0.0) |
| 21/7/7   |  38.7% |  20.1% |   8.6% | +1.6 (-0.6 to +3.7) | +2.0 (+0.3 to +3.7) |  -4.1 (-5.5 to -2.8) |
| 21/0/0   |  40.2% |  20.5% |   3.0% | +3.1 (+1.0 to +5.3) | +2.4 (+0.7 to +4.2) | -9.7 (-10.9 to -8.6) |

*60% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |    Difference P(≥50) |
|:---------|-------:|-------:|-------:|--------------------:|--------------------:|---------------------:|
| 21/21/21 |  28.0% |  12.8% |  10.5% |           Reference |           Reference |            Reference |
| 21/14/14 |  28.7% |  12.9% |   8.1% | +0.8 (-1.2 to +2.7) | +0.1 (-1.4 to +1.5) |  -2.4 (-3.7 to -1.2) |
| 21/7/7   |  28.2% |  13.4% |   5.7% | +0.2 (-1.8 to +2.2) | +0.5 (-1.0 to +2.0) |  -4.9 (-6.1 to -3.7) |
| 21/0/0   |  29.3% |  12.3% |   1.5% | +1.4 (-0.6 to +3.3) | -0.5 (-2.0 to +0.9) | -9.0 (-10.1 to -8.0) |

*70% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |  16.0% |   6.4% |   6.0% |           Reference |           Reference |           Reference |
| 21/14/14 |  18.1% |   7.1% |   5.8% | +2.1 (+0.5 to +3.8) | +0.7 (-0.4 to +1.8) | -0.2 (-1.3 to +0.8) |
| 21/7/7   |  17.8% |   6.2% |   3.5% | +1.9 (+0.2 to +3.5) | -0.2 (-1.3 to +0.9) | -2.5 (-3.4 to -1.6) |
| 21/0/0   |  19.1% |   5.6% |   0.8% | +3.2 (+1.5 to +4.9) | -0.9 (-1.9 to +0.2) | -5.2 (-6.0 to -4.5) |

*80% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |   7.5% |   4.0% |   2.1% |           Reference |           Reference |           Reference |
| 21/14/14 |   7.7% |   2.9% |   1.2% | +0.2 (-1.0 to +1.4) | -1.1 (-1.9 to -0.3) | -0.8 (-1.4 to -0.2) |
| 21/7/7   |   7.3% |   1.9% |   0.9% | -0.2 (-1.3 to +1.0) | -2.1 (-2.8 to -1.3) | -1.1 (-1.7 to -0.6) |
| 21/0/0   |   8.5% |   0.8% |   0.1% | +1.0 (-0.2 to +2.2) | -3.2 (-3.9 to -2.6) | -1.9 (-2.4 to -1.5) |

*90% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |   2.0% |   1.5% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |   1.6% |   1.1% |   0.0% | -0.4 (-1.0 to +0.2) | -0.4 (-0.9 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |   1.3% |   0.5% |   0.0% | -0.7 (-1.3 to -0.1) | -1.0 (-1.5 to -0.6) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |   0.9% |   0.1% |   0.0% | -1.1 (-1.6 to -0.6) | -1.4 (-1.8 to -1.0) | +0.0 (-0.1 to +0.1) |

*95% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |   1.3% |   0.0% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |   0.6% |   0.0% |   0.0% | -0.7 (-1.2 to -0.3) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |   0.3% |   0.0% |   0.0% | -1.0 (-1.4 to -0.6) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |   0.1% |   0.0% |   0.0% | -1.2 (-1.6 to -0.9) | +0.0 (-0.1 to +0.1) | +0.0 (-0.1 to +0.1) |

### High

**No vaccination: which single tier to quarantine** (compared with no
quarantine)

| Scenario             | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |      Difference P(≥50) |
|:---------------------|-------:|-------:|-------:|--------------------:|--------------------:|-----------------------:|
| No quarantine        |  98.1% |  98.1% |  98.0% |           Reference |           Reference |              Reference |
| Only high-risk (21d) |  94.1% |  92.9% |  86.7% | -4.0 (-4.9 to -3.2) | -5.2 (-6.1 to -4.3) | -11.3 (-12.5 to -10.2) |
| Only med-risk (21d)  |  94.4% |  94.3% |  94.2% | -3.7 (-4.6 to -2.9) | -3.7 (-4.6 to -2.9) |    -3.8 (-4.7 to -3.0) |
| Only low-risk (21d)  |  96.1% |  93.2% |  80.8% | -2.0 (-2.7 to -1.3) | -4.9 (-5.8 to -4.0) | -17.2 (-18.5 to -15.9) |

**Tiered strategies (high/medium/low days)**, compared with 21/21/21 at
the same coverage

*50% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |      Difference P(≥50) |
|:---------|-------:|-------:|-------:|-----------------------:|-----------------------:|-----------------------:|
| 21/21/21 |  62.3% |  48.1% |  30.4% |              Reference |              Reference |              Reference |
| 21/14/14 |  65.1% |  50.8% |  31.4% |    +2.8 (+0.7 to +4.9) |    +2.7 (+0.6 to +4.9) |    +1.1 (-1.0 to +3.1) |
| 21/7/7   |  71.5% |  60.2% |  41.0% |   +9.2 (+7.1 to +11.2) |  +12.1 (+9.9 to +14.3) |  +10.7 (+8.6 to +12.8) |
| 21/0/0   |  81.0% |  76.0% |  55.2% | +18.7 (+16.8 to +20.6) | +27.9 (+25.8 to +29.9) | +24.8 (+22.7 to +26.9) |

*60% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |      Difference P(≥50) |
|:---------|-------:|-------:|-------:|-----------------------:|-----------------------:|-----------------------:|
| 21/21/21 |  52.7% |  36.6% |  21.6% |              Reference |              Reference |              Reference |
| 21/14/14 |  57.0% |  41.0% |  25.1% |    +4.3 (+2.1 to +6.5) |    +4.4 (+2.2 to +6.5) |    +3.4 (+1.6 to +5.3) |
| 21/7/7   |  63.8% |  51.6% |  28.2% |  +11.1 (+8.9 to +13.2) | +15.0 (+12.8 to +17.1) |    +6.6 (+4.7 to +8.5) |
| 21/0/0   |  72.5% |  64.7% |  38.8% | +19.8 (+17.7 to +21.9) | +28.1 (+26.0 to +30.2) | +17.2 (+15.2 to +19.1) |

*70% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|-----------------------:|-----------------------:|--------------------:|
| 21/21/21 |  41.7% |  24.5% |  16.6% |              Reference |              Reference |           Reference |
| 21/14/14 |  44.7% |  28.4% |  16.1% |    +3.0 (+0.8 to +5.2) |    +3.9 (+2.0 to +5.9) | -0.5 (-2.1 to +1.1) |
| 21/7/7   |  52.1% |  36.8% |  17.4% |  +10.4 (+8.2 to +12.5) | +12.3 (+10.3 to +14.3) | +0.9 (-0.8 to +2.5) |
| 21/0/0   |  62.4% |  49.5% |  16.8% | +20.7 (+18.5 to +22.8) | +25.0 (+23.0 to +27.1) | +0.2 (-1.4 to +1.8) |

*80% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |      Difference P(≥10) |      Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|-----------------------:|-----------------------:|--------------------:|
| 21/21/21 |  25.4% |  11.6% |   4.3% |              Reference |              Reference |           Reference |
| 21/14/14 |  29.2% |  14.9% |   3.3% |    +3.9 (+1.9 to +5.8) |    +3.2 (+1.8 to +4.7) | -1.0 (-1.8 to -0.1) |
| 21/7/7   |  34.7% |  18.5% |   4.2% |   +9.4 (+7.4 to +11.4) |    +6.9 (+5.3 to +8.4) | -0.1 (-1.0 to +0.7) |
| 21/0/0   |  44.4% |  26.6% |   2.5% | +19.1 (+17.0 to +21.1) | +15.0 (+13.3 to +16.7) | -1.8 (-2.6 to -1.0) |

*90% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |   8.9% |   4.7% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |  10.2% |   4.5% |   0.0% | +1.3 (-0.0 to +2.6) | -0.2 (-1.1 to +0.7) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |  13.4% |   3.5% |   0.0% | +4.5 (+3.1 to +5.9) | -1.2 (-2.1 to -0.4) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |  16.4% |   3.5% |   0.0% | +7.5 (+6.0 to +8.9) | -1.2 (-2.1 to -0.3) | +0.0 (-0.1 to +0.1) |

*95% vaccinated*

| Strategy | P(≥10) | P(≥20) | P(≥50) |   Difference P(≥10) |   Difference P(≥20) |   Difference P(≥50) |
|:---------|-------:|-------:|-------:|--------------------:|--------------------:|--------------------:|
| 21/21/21 |   2.7% |   0.1% |   0.0% |           Reference |           Reference |           Reference |
| 21/14/14 |   3.0% |   0.2% |   0.0% | +0.3 (-0.4 to +1.1) | +0.1 (-0.1 to +0.3) | +0.0 (-0.1 to +0.1) |
| 21/7/7   |   2.6% |   0.0% |   0.0% | -0.1 (-0.8 to +0.7) | -0.1 (-0.2 to +0.1) | +0.0 (-0.1 to +0.1) |
| 21/0/0   |   3.8% |   0.1% |   0.0% | +1.1 (+0.4 to +1.9) | -0.0 (-0.2 to +0.1) | +0.0 (-0.1 to +0.1) |
