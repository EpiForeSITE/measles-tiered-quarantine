# Measles tiered quarantine in schools

Agent-based simulations of measles spread in elementary, middle and high
schools, testing **tiered quarantine**: contacts of a case are quarantined
for a different number of days depending on whether they are high-, medium-
or low-risk contacts. Strategies are written as **high / medium / low days**,
e.g. `21/7/0`.

The model uses [`epiworldR`](https://github.com/UofUEpiBio/epiworldR) and the
[`measles`](https://github.com/UofUEpiBio/measles) package
(`ModelMeaslesMixingRiskQuarantine`).

## Files

| File | What it does |
|:--|:--|
| `measles_model.R` | Shared code: settings, contact matrices, transmission rates, `simulator()`, design helpers and plot theme. Sourced by every `.qmd`. |
| `01_tiered_quarantine_tables.qmd` | R0 check and tables of P(outbreak ≥ 10 / 20 / 50) for single-tier and tiered strategies across vaccination coverage, with Wilcoxon rank-sum p-values against a reference strategy. |
| `02_outbreak_size_distributions.qmd` | Full outbreak-size distributions for 11 strategies: density curves comparing schools, histograms per school, and summary tables. |
| `03_quarantine_policy_analysis.qmd` | Main analysis over every allowed high/medium/low combination: baseline risk, what drives risk (decision tree), policy options, optimal low-risk days, outbreak sizes, threshold curves, and plain-language conclusions. |

Each document runs on its own; there is no required order.

Rendering a `.qmd` produces two versions of the report:

- a **`.md`** file (e.g. `03_quarantine_policy_analysis.md`) with its
  figures in `md_figures/`, which GitHub displays directly. It shows the
  results only; the code is in the `.qmd`.
- a self-contained **`.html`** file with folded code, for viewing locally
  or sharing (git-ignored).

Commit the `.md` files and `md_figures/` so the results can be read on
GitHub without running anything.

## Model setup

- **Contact matrices:** real per-capita SMART contact matrices by grade,
  symmetrised for reciprocity and kept at real magnitude.
- **Enrollment:** 26 students per grade (elementary, K–5), 59 (middle, 5–8),
  58 (high, 9–12).
- **Transmission:** R0 = 12 is reached with a per-school transmission rate,
  `p = R0 / (infectious days × largest eigenvalue)`, with 4 infectious days.
- **Runs:** one index case, 100 days, 2,000 simulations per setting,
  seed 221. Outbreak size = everyone no longer susceptible at day 100.
- **Quarantine:** contact tracing window of 7 days, detection rate 0.
  Except for the diagnostic scenarios in document 1, all strategies satisfy
  high ≥ medium ≥ low.

Note on thresholds: documents 1 and 2 report P(size **≥** K), while
document 3 defines an outbreak as **more than** K cases, following the
original scripts.

## Requirements

- R ≥ 4.1 and [Quarto](https://quarto.org)
- R packages: `epiworldR`, `measles`, `data.table`, `ggplot2`, `scales`,
  `rpart`, `knitr`, `rmarkdown`

```r
install.packages(c("epiworldR", "data.table", "ggplot2", "scales",
                   "rpart", "knitr", "rmarkdown", "remotes"))
remotes::install_github("UofUEpiBio/measles")
```

## Running

```bash
quarto render 01_tiered_quarantine_tables.qmd
quarto render 02_outbreak_size_distributions.qmd
quarto render 03_quarantine_policy_analysis.qmd
```

Simulations are slow at full size. For a quick test, lower the number of
simulations with an environment variable:

```bash
N_SIMS=200 quarto render 03_quarantine_policy_analysis.qmd
```

On a SLURM cluster the number of threads is taken from
`SLURM_CPUS_PER_TASK`; otherwise all detected cores are used. If
`~/R_libs` exists it is used as the package library.

## Outputs

Besides the `.md`/`.html` reports and `md_figures/`, running a document
creates these folders (both are git-ignored):

- `results/` holds the raw simulation output, cached as
  `<analysis>_n<N_SIMS>.csv`. If the file exists, the document loads it
  instead of re-simulating; delete it to re-run. Document 3 also writes
  the probability tables behind Figure 6 and `conclusions.txt`.
- `figures/` holds a PNG copy of every figure.

To reuse raw output from the earlier stand-alone scripts, move
`all_school_simulation_sizes.csv` to `results/policy_grid_n2000.csv`
(same columns).
