# 3. Which quarantine policy keeps outbreaks small?

- [Simulations](#simulations)
- [1. The problem](#the-problem)
- [2. What drives outbreak risk?](#what-drives-outbreak-risk)
- [3. Which quarantine combinations keep risk
  low?](#which-quarantine-combinations-keep-risk-low)
- [4. How many low-risk quarantine days are worth
  it?](#how-many-low-risk-quarantine-days-are-worth-it)
- [5. How big do outbreaks get?](#how-big-do-outbreaks-get)
- [Supplement: every quarantine combination, per
  school](#supplement-every-quarantine-combination-per-school)
- [Conclusions](#conclusions)

This is the main analysis. Every combination of school, vaccination
coverage and **high / medium / low-risk quarantine days** (with high ≥
medium ≥ low) is simulated once, and all figures are computed from that
single set of runs.

| Figure | Question                                                                        |
|:-------|:--------------------------------------------------------------------------------|
| 1      | How risky is each school with minimal quarantine?                               |
| 2      | What happens to outbreak risk when one factor changes, and is the change clear? |
| 3      | Which quarantine combinations keep risk low?                                    |
| 4      | Given high- and medium-risk days, how many low-risk days are worth it?          |
| 5      | How big do outbreaks get, with high-risk-only and with optimal quarantine?      |
| 6      | How does risk change with low-risk days, across thresholds?                     |

Plain-language conclusions, computed from the results, are at the end
and are also written to `results/conclusions.txt`.

<div>

> **Note**
>
> Model code lives in `measles_model.R`. Simulations are cached in
> `results/`, and PNG copies of every figure are saved to `figures/`.

</div>

## Simulations

136 quarantine/vaccination settings per school × 3 schools × 4,000
simulations = 1,632,000 runs.

<div>

> **What “outbreak risk” means in this report**
>
> Each simulation starts with one infected student. An **outbreak** is a
> simulation that ends with **more than 25 total cases**. The **outbreak
> risk** of a setting is the share of its 4,000 simulations that ended
> this way, i.e. P(final outbreak size \> 25). Figures 1 to 5 and the
> conclusions use this definition; Figure 6 repeats it for several other
> thresholds.

</div>

## 1. The problem

With minimal quarantine (only high-risk contacts quarantined), how
likely is an outbreak?

![](md_figures/03_quarantine_policy_analysis-fig1-1.png)

## 2. What drives outbreak risk?

Each panel changes **one** factor while the others stay fixed, and shows
the outbreak risk before (grey dot) and after (arrow). The fixed values
are named in each panel’s title. Vaccination is compared at one fixed
quarantine policy; each quarantine tier is compared at the lowest and
highest vaccination level.

Each risk comes from a finite number of simulations, so two settings can
differ a little just by chance. Arrows are therefore coloured by whether
the change is **clear**: blue if the risk is clearly lower, red if
clearly higher, and grey if the 95% confidence interval of the change
includes 0 (the change can’t be told apart from simulation noise). The
table below the figure gives the numbers.

![](md_figures/03_quarantine_policy_analysis-fig2a-1.png)

**Small benefits, with 95% confidence intervals.** Large changes are
clear from the figure, so this table only lists the comparisons where
the **benefit is under 5 percentage points**, where simulation noise
matters most. Benefit = risk before − risk after, in percentage points:
positive = lower risk, negative = higher risk. *Clearly lower / higher*
means the whole interval is above / below 0; *No clear change* means the
interval includes 0. Intervals reflect simulation noise only, not
uncertainty in the model’s assumptions.

*Increasing vaccination: 60% → 85% (quarantine 14/7/0 days)*

| Setting | Risk before | Risk after | Benefit, points (95% CI) | Verdict       |
|:--------|------------:|-----------:|-------------------------:|:--------------|
| Middle  |        1.4% |       0.0% |      +1.4 (+1.1 to +1.8) | Clearly lower |

*Increasing high-risk quarantine: 14 → 21 days (medium 14, low 0 days)*

| Setting                     | Risk before | Risk after | Benefit, points (95% CI) | Verdict         |
|:----------------------------|------------:|-----------:|-------------------------:|:----------------|
| Elementary · 85% vaccinated |        0.0% |       0.0% |      +0.0 (-0.1 to +0.1) | No clear change |
| Middle · 60% vaccinated     |        1.5% |       0.4% |      +1.1 (+0.7 to +1.6) | Clearly lower   |
| Middle · 85% vaccinated     |        0.0% |       0.0% |      -0.0 (-0.1 to +0.1) | No clear change |
| High · 85% vaccinated       |        1.0% |       0.2% |      +0.8 (+0.5 to +1.2) | Clearly lower   |

*Increasing medium-risk quarantine: 0 → 14 days (high 14, low 0 days)*

| Setting                     | Risk before | Risk after | Benefit, points (95% CI) | Verdict         |
|:----------------------------|------------:|-----------:|-------------------------:|:----------------|
| Elementary · 85% vaccinated |        0.0% |       0.0% |      +0.0 (-0.1 to +0.2) | No clear change |
| Middle · 60% vaccinated     |        1.5% |       1.5% |      -0.0 (-0.5 to +0.5) | No clear change |
| Middle · 85% vaccinated     |        0.0% |       0.0% |      -0.0 (-0.1 to +0.1) | No clear change |
| High · 85% vaccinated       |        1.8% |       1.0% |      +0.8 (+0.3 to +1.3) | Clearly lower   |

*Increasing low-risk quarantine: 0 → 7 days (high 14, medium 7 days)*

| Setting                     | Risk before | Risk after | Benefit, points (95% CI) | Verdict         |
|:----------------------------|------------:|-----------:|-------------------------:|:----------------|
| Elementary · 85% vaccinated |        0.0% |       0.0% |      +0.0 (-0.1 to +0.1) | No clear change |
| Middle · 60% vaccinated     |        1.4% |       1.5% |      -0.1 (-0.6 to +0.5) | No clear change |
| Middle · 85% vaccinated     |        0.0% |       0.0% |      -0.0 (-0.1 to +0.1) | No clear change |
| High · 85% vaccinated       |        1.0% |       0.6% |      +0.4 (+0.0 to +0.8) | Clearly lower   |

## 3. Which quarantine combinations keep risk low?

Every quarantine combination, coloured by its outbreak risk on a
continuous scale: **light = low risk, dark = high risk**. The number in
each square is its risk. Squares with a **black outline** have a risk of
10% or less. Columns are medium / low-risk days, grouped by high-risk
days; low-risk days are shown in weekly steps to keep the figure
readable.

![](md_figures/03_quarantine_policy_analysis-fig3-1.png)

## 4. How many low-risk quarantine days are worth it?

For each fixed high/medium policy: the fewest low-risk days that give
80% of the benefit low-risk quarantine can provide. Each cell shows
those days and the risk without → with them. When low-risk quarantine
cannot lower the risk by 2% or more, the cell shows **0 days**, the risk
with no low-risk quarantine, and in brackets the largest drop any number
of low-risk days gave.

![](md_figures/03_quarantine_policy_analysis-fig4-1.png)

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

<div class="panel-tabset">

#### Elementary

| Vaccinated | High / medium days | Low-risk days compared | Risk, 0 days | Risk, with days | Benefit, points (95% CI) | Clearly under 2 points? |
|:-----------|:-------------------|:-----------------------|-------------:|----------------:|-------------------------:|:-----------------------:|
| 80%        | 21 / 21            | 0 vs 7                 |         0.1% |            0.0% |      +0.1 (-0.0 to +0.2) |           Yes           |
| 80%        | 21 / 14            | 0 vs 14                |         0.2% |            0.0% |      +0.1 (+0.0 to +0.3) |           Yes           |
| 80%        | 21 / 7             | 0 vs 5                 |         0.4% |            0.1% |      +0.3 (+0.1 to +0.6) |           Yes           |
| 80%        | 14 / 14            | 0 vs 10                |         0.4% |            0.0% |      +0.3 (+0.1 to +0.6) |           Yes           |
| 80%        | 14 / 7             | 0 vs 7                 |         0.7% |            0.2% |      +0.5 (+0.3 to +0.9) |           Yes           |
| 80%        | 7 / 7              | 0 vs 7                 |         1.8% |            0.5% |      +1.3 (+0.8 to +1.8) |           Yes           |
| 85%        | 21 / 21            | 0 vs 21                |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 85%        | 21 / 14            | 0 vs 14                |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 85%        | 21 / 7             | 0 vs 3                 |         0.0% |            0.0% |      +0.1 (-0.1 to +0.2) |           Yes           |
| 85%        | 14 / 14            | 0 vs 3                 |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 85%        | 14 / 7             | 0 vs 5                 |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 85%        | 7 / 7              | 0 vs 5                 |         0.1% |            0.0% |      +0.1 (-0.0 to +0.2) |           Yes           |

#### Middle

| Vaccinated | High / medium days | Low-risk days compared | Risk, 0 days | Risk, with days | Benefit, points (95% CI) | Clearly under 2 points? |
|:-----------|:-------------------|:-----------------------|-------------:|----------------:|-------------------------:|:-----------------------:|
| 60%        | 21 / 21            | 0 vs 10                |         0.5% |            0.4% |      +0.1 (-0.2 to +0.4) |           Yes           |
| 60%        | 21 / 14            | 0 vs 14                |         0.4% |            0.4% |      -0.1 (-0.4 to +0.2) |           Yes           |
| 60%        | 21 / 7             | 0 vs 7                 |         0.7% |            0.3% |      +0.4 (+0.1 to +0.8) |           Yes           |
| 60%        | 14 / 14            | 0 vs 14                |         1.5% |            1.0% |      +0.5 (+0.1 to +1.0) |           Yes           |
| 60%        | 14 / 7             | 0 vs 5                 |         1.4% |            1.1% |      +0.3 (-0.2 to +0.8) |           Yes           |
| 60%        | 7 / 7              | 0 vs 7                 |         4.9% |            3.7% |      +1.2 (+0.3 to +2.1) |           No            |
| 70%        | 21 / 21            | 0 vs 3                 |         0.0% |            0.0% |      +0.1 (-0.1 to +0.2) |           Yes           |
| 70%        | 21 / 14            | 0 vs 5                 |         0.1% |            0.0% |      +0.1 (-0.0 to +0.2) |           Yes           |
| 70%        | 21 / 7             | 0 vs 3                 |         0.0% |            0.0% |      +0.1 (-0.1 to +0.2) |           Yes           |
| 70%        | 14 / 14            | 0 vs 7                 |         0.1% |            0.0% |      +0.1 (-0.0 to +0.2) |           Yes           |
| 70%        | 14 / 7             | 0 vs 7                 |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 70%        | 7 / 7              | 0 vs 7                 |         0.3% |            0.2% |      +0.1 (-0.1 to +0.4) |           Yes           |
| 80%        | 21 / 21            | 0 vs 21                |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 80%        | 21 / 14            | 0 vs 14                |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 80%        | 21 / 7             | 0 vs 7                 |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 80%        | 14 / 14            | 0 vs 14                |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 80%        | 14 / 7             | 0 vs 7                 |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 80%        | 7 / 7              | 0 vs 7                 |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 85%        | 21 / 21            | 0 vs 21                |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 85%        | 21 / 14            | 0 vs 14                |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 85%        | 21 / 7             | 0 vs 7                 |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 85%        | 14 / 14            | 0 vs 14                |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 85%        | 14 / 7             | 0 vs 7                 |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |
| 85%        | 7 / 7              | 0 vs 7                 |         0.0% |            0.0% |      +0.0 (-0.1 to +0.1) |           Yes           |

#### High

| Vaccinated | High / medium days | Low-risk days compared | Risk, 0 days | Risk, with days | Benefit, points (95% CI) | Clearly under 2 points? |
|:-----------|:-------------------|:-----------------------|-------------:|----------------:|-------------------------:|:-----------------------:|
| 80%        | 21 / 21            | 0 vs 14                |         1.8% |            1.3% |      +0.5 (+0.0 to +1.1) |           Yes           |
| 80%        | 21 / 14            | 0 vs 10                |         2.7% |            1.0% |      +1.7 (+1.1 to +2.3) |           No            |
| 80%        | 21 / 7             | 0 vs 7                 |         3.9% |            2.2% |      +1.7 (+1.0 to +2.5) |           No            |
| 85%        | 21 / 21            | 0 vs 14                |         0.3% |            0.1% |      +0.2 (+0.1 to +0.5) |           Yes           |
| 85%        | 21 / 14            | 0 vs 10                |         0.2% |            0.0% |      +0.2 (+0.0 to +0.4) |           Yes           |
| 85%        | 21 / 7             | 0 vs 3                 |         0.5% |            0.3% |      +0.2 (-0.1 to +0.5) |           Yes           |
| 85%        | 14 / 14            | 0 vs 7                 |         1.0% |            0.2% |      +0.8 (+0.5 to +1.2) |           Yes           |
| 85%        | 14 / 7             | 0 vs 5                 |         1.0% |            0.4% |      +0.6 (+0.2 to +1.0) |           Yes           |

</div>

## 5. How big do outbreaks get?

<div class="panel-tabset">

### 5a. Only high-risk contacts quarantined

![](md_figures/03_quarantine_policy_analysis-fig5a-1.png)

### 5b. Optimal quarantine

**How the optimal policy is computed.** For one school and one
vaccination level, let $P(h, m, l)$ be the outbreak risk with $h$ high-,
$m$ medium- and $l$ low-risk quarantine days: the share of simulations
with more than 25 cases. Only policies with $h \ge m \ge l$ are
considered.

1.  **Risk with minimal quarantine:** $P_0 = P(7, 0, 0)$.
2.  **Lowest risk any tested policy reaches:**
    $P_{\min} = \min_{h, m, l} P(h, m, l)$.
3.  **Largest possible reduction:** $G = P_0 - P_{\min}$.
4.  **If quarantine barely helps** ($G <$ 2%), the optimal policy is
    minimal quarantine, 7 / 0 / 0.
5.  **Otherwise, a policy is good enough** if it gets at least 80% of
    the largest reduction: $$P_0 - P(h, m, l) \ge 0.8 \times G$$
6.  **The optimal policy** is the good-enough policy with the fewest
    total quarantine days $h + m + l$. If several tie, the one with the
    lowest risk is chosen.

*Example:* Elementary school, 60% vaccinated. $P_0$ = 60.5% and
$P_{\min}$ = 6.4%, so $G$ = 54.2 percentage points. A good-enough policy
must bring the risk down to 17.2% or lower (60.5% − 0.8 × 54.2 points).
The shortest such policy is **21/7/5** (33 days in total), with a risk
of 15.2%.

![](md_figures/03_quarantine_policy_analysis-fig5b-1.png)

</div>

## Supplement: every quarantine combination, per school

<div class="panel-tabset">

### Elementary

### Middle

![](md_figures/03_quarantine_policy_analysis-supplement-1.png)

### High

![](md_figures/03_quarantine_policy_analysis-supplement-2.png)![](md_figures/03_quarantine_policy_analysis-supplement-3.png)

</div>

## Conclusions

``` text
CONCLUSIONS  (outbreak = more than 25 cases)

1. THE PROBLEM
   With minimal quarantine, High school is the riskiest (average 43%)
   and Middle school the least risky (average 2%).

2. WHAT MATTERS MOST
   Increasing vaccination: 60% -> 85%  (quarantine 14/7/0 days):
     Elementary: risk 34% -> 0% (change -33.7 points, 95% CI -35.1 to -32.2: clearly lower)
     Middle: risk 1% -> 0% (change -1.4 points, 95% CI -1.8 to -1.1: clearly lower)
     High: risk 54% -> 1% (change -52.5 points, 95% CI -54.1 to -50.9: clearly lower)
   Increasing high-risk quarantine: 14 -> 21 days  (medium 14, low 0 days):
     Elementary - 60% vaccinated: risk 25% -> 17% (change -7.8 points, 95% CI -9.5 to -6.0: clearly lower)
     Elementary - 85% vaccinated: risk 0% -> 0% (change -0.0 points, 95% CI -0.1 to +0.1: no clear change)
     Middle - 60% vaccinated: risk 2% -> 0% (change -1.1 points, 95% CI -1.6 to -0.7: clearly lower)
     Middle - 85% vaccinated: risk 0% -> 0% (change +0.0 points, 95% CI -0.1 to +0.1: no clear change)
     High - 60% vaccinated: risk 46% -> 33% (change -13.2 points, 95% CI -15.3 to -11.0: clearly lower)
     High - 85% vaccinated: risk 1% -> 0% (change -0.8 points, 95% CI -1.2 to -0.5: clearly lower)
   Increasing medium-risk quarantine: 0 -> 14 days  (high 14, low 0 days):
     Elementary - 60% vaccinated: risk 43% -> 25% (change -18.4 points, 95% CI -20.5 to -16.4: clearly lower)
     Elementary - 85% vaccinated: risk 0% -> 0% (change -0.0 points, 95% CI -0.2 to +0.1: no clear change)
     Middle - 60% vaccinated: risk 2% -> 2% (change +0.0 points, 95% CI -0.5 to +0.5: no clear change)
     Middle - 85% vaccinated: risk 0% -> 0% (change +0.0 points, 95% CI -0.1 to +0.1: no clear change)
     High - 60% vaccinated: risk 65% -> 46% (change -19.3 points, 95% CI -21.4 to -17.1: clearly lower)
     High - 85% vaccinated: risk 2% -> 1% (change -0.8 points, 95% CI -1.3 to -0.3: clearly lower)
   Increasing low-risk quarantine: 0 -> 7 days  (high 14, medium 7 days):
     Elementary - 60% vaccinated: risk 34% -> 18% (change -15.7 points, 95% CI -17.6 to -13.8: clearly lower)
     Elementary - 85% vaccinated: risk 0% -> 0% (change -0.0 points, 95% CI -0.1 to +0.1: no clear change)
     Middle - 60% vaccinated: risk 1% -> 2% (change +0.1 points, 95% CI -0.5 to +0.6: no clear change)
     Middle - 85% vaccinated: risk 0% -> 0% (change +0.0 points, 95% CI -0.1 to +0.1: no clear change)
     High - 60% vaccinated: risk 54% -> 41% (change -12.1 points, 95% CI -14.3 to -9.9: clearly lower)
     High - 85% vaccinated: risk 1% -> 1% (change -0.4 points, 95% CI -0.8 to -0.0: clearly lower)

3. WHAT A SCHOOL SHOULD DO
   To keep risk under 10%:
     Elementary: achievable from 60% vaccination; minimal quarantine is enough from 80% vaccination.
     Middle: achievable from 60% vaccination; minimal quarantine is enough from 60% vaccination.
     High: achievable from 70% vaccination; minimal quarantine is enough from 85% vaccination.
   To keep risk under 20%:
     Elementary: achievable from 60% vaccination; minimal quarantine is enough from 80% vaccination.
     Middle: achievable from 60% vaccination; minimal quarantine is enough from 60% vaccination.
     High: achievable from 60% vaccination; minimal quarantine is enough from 85% vaccination.
   To keep risk under 30%:
     Elementary: achievable from 60% vaccination; minimal quarantine is enough from 80% vaccination.
     Middle: achievable from 60% vaccination; minimal quarantine is enough from 60% vaccination.
     High: achievable from 60% vaccination; minimal quarantine is enough from 80% vaccination.

4. LOW-RISK QUARANTINE
   Elementary: the best low-risk quarantine is usually about 2.5 days.
     It is not worth it (risk drops less than 2%) in 50% of high/medium settings.
     On average, risk goes 9% -> 5% with the best number of days,
     and only to 5% even with the longest low-risk quarantine.
   Middle: the best low-risk quarantine is usually about 0 days.
     It is not worth it (risk drops less than 2%) in 100% of high/medium settings.
     On average, risk goes 0% -> 0% with the best number of days,
     and only to 0% even with the longest low-risk quarantine.
   High: the best low-risk quarantine is usually about 5 days.
     It is not worth it (risk drops less than 2%) in 33% of high/medium settings.
     On average, risk goes 20% -> 15% with the best number of days,
     and only to 14% even with the longest low-risk quarantine.

5. SIZE OF OUTBREAKS
   Elementary, minimal quarantine: chance of more than 50 cases falls from 6% at 60% to 0% at 85% vaccination.
     Optimal quarantine (high/medium/low days): 21/7/5 at 60%; 14/7/5 at 70%; 7/7/5 at 80%; 7/0/0 at 85%.
     With it, the chance of more than 50 cases is 0% at 60% and 0% at 85% vaccination.
   Middle, minimal quarantine: chance of more than 50 cases falls from 0% at 60% to 0% at 85% vaccination.
     Optimal quarantine (high/medium/low days): 21/0/0 at 60%; 7/0/0 at 70%; 7/0/0 at 80%; 7/0/0 at 85%.
     With it, the chance of more than 50 cases is 0% at 60% and 0% at 85% vaccination.
   High, minimal quarantine: chance of more than 50 cases falls from 45% at 60% to 0% at 85% vaccination.
     Optimal quarantine (high/medium/low days): 21/14/5 at 60%; 21/7/3 at 70%; 14/7/5 at 80%; 21/0/0 at 85%.
     With it, the chance of more than 50 cases is 1% at 60% and 0% at 85% vaccination.
```
