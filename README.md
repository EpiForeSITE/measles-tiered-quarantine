# Tiered quarantine for measles exposures in schools

Agent-based simulations of measles outbreaks in an elementary, a middle and a
high school, comparing **tiered quarantine strategies**: different quarantine
lengths for high-, medium- and low-risk contacts of a detected case.

A strategy is written **H/M/L** = days at home for

| Tier | Who |
|---|---|
| **H**igh risk | unvaccinated students in the same grade as the detected case |
| **M**edium risk | unvaccinated recent contacts of the detected case (contact tracing) |
| **L**ow risk | all other unvaccinated students |

For example, `21/14/0` = 21 days for same-grade contacts, 14 days for traced
contacts, none for the rest. `21/21/21` (everyone unvaccinated home for 21
days) is the reference. Vaccinated students are never quarantined.

The model is `ModelMeaslesMixingRiskQuarantine` from the
[measles](https://github.com/UofUEpiBio/measles) R package (built on
[epiworldR](https://github.com/UofUEpiBio/epiworldR)).

---

## Model set-up (`measles_model.R`)

Every analysis starts with `source("measles_model.R")` (or
`measles_model_dur.R`, below). It holds the shared settings:

| Setting | Value |
|---|---|
| Schools | Elementary (grades K–5), Middle (5–8), High (9–12) |
| Students | Participants per grade in the SMART sensor study (`ENROLLMENT`): 339 / 466 / 473 |
| Contacts | Per-capita SMART contact matrices between grades, reciprocity-symmetrised (`SCHOOL_CMAT`) |
| R0 | 12, reached with one transmission rate per school: β = R0 / (contagious days × largest eigenvalue) |
| Contagious period | 4 days, before the rash only |
| Each run | 1 index case, 100 days |
| Outcome | Final outbreak size; P(≥10 cases), P(≥20 cases) |

Contact data: Grantz et al., *Sci Rep* 2021 and Guclu et al., *PLoS ONE* 2016
(SMART study; a contact = another sensor within about 3 m).

---

## Two versions of the disease model

| | Original package | Patched package (`_dur` files) |
|---|---|---|
| Latent period | Geometric (memoryless), mean 12 days | **Log-normal**, median 12.5 days, dispersion 1.23 (Lessler et al. 2009), drawn once at infection |
| Contagious period | Geometric, mean 4 days | **Gamma**, mean 4, SD 1, drawn once |
| Contact probability | C/N (can exceed 1) | min(1, C/N) |
| Extra options | – | contagious rash (`Rash contagious days`), per-tier quarantine counters (`Tracked ...`) |

The patched model is `measlesmixingriskquarantine.hpp`. All new settings are
off by default, so with them off it gives the same results as the original.
`measles_model_dur.R` loads the patched package, sources `measles_model.R`
unchanged, and switches every model to the new durations.

---

## Files

### Shared code

| File | What it is |
|---|---|
| `measles_model.R` | Settings, contact matrices, `simulator()`, design helpers, Wilson / Newcombe CIs, plot theme |
| `measles_model_dur.R` | Same, with the log-normal / gamma durations. Options: `R0`, `RASH_DAYS`, `SCHOOL_SIZE` (below) |
| `measlesmixingriskquarantine.hpp` | Patched model (copy into the measles source before installing) |
| `run_dur.sh` | Slurm script: installs the patched package once, renders the `_dur` qmd files |
| `install_crn.R` | Installs the "fair pairs" (common random numbers) version into `~/R_libs_crn` |

### Analyses (Quarto)

| File | Question |
|---|---|
| `01_tiered_quarantine_tables.qmd` | Outbreak probability for each strategy, school and vaccination level |
| `02_outbreak_size_distributions.qmd` | Distribution of final outbreak sizes (histograms) |
| `03_quarantine_policy_analysis.qmd` | Which policy keeps outbreaks small? Vaccination vs quarantine, one change at a time |
| `04_model_mechanics.qmd` | How the model works: anatomy of an outbreak |
| `07_fair_pairs_tables.qmd` | Strategies compared on the same outbreak (common random numbers) |
| `08_independent_tables.qmd` | Strategies compared with independent runs |
| `09_model_behaviour.qmd` | Disease progression, quarantine, where infections come from |
| `10_quarantine_burden_dur.qmd` | What quarantine costs a school: rounds, students sent home, school days lost |
| `10b_outbreak_and_burden_table_dur.qmd` | One slide table: outbreak risk, benefit vs 7/0/0, and burden |
| `sensitivity.qmd` | Sensitivity analysis: tornado bars at 60/70/80% coverage and cost–risk arrows |

Files ending in `_dur.qmd` are the same analyses with the patched model
(they source `measles_model_dur.R`).

### Helper scripts (`Rscript <file>`, next to `measles_model.R`)

| File | What it does |
|---|---|
| `student_timeline.R` | Day-by-day state of every student in one outbreak |
| `outbreak_tree.R` | Who infected whom (transmission tree) |
| `quarantine_spread.R` | Where infections come from and whether quarantine helped |
| `concentration_test.R` | Tests whether students left at school mix more ("concentration") |
| `find_paired_runs.R` | Seeds where a stricter strategy did worse / better / the same |
| `fair_pair.R`, `trace_pair.R` | Fair comparison of two strategies on one outbreak; follow one infection |

---

## Installing (CHPC)

Three separate R libraries keep the versions apart:

| Folder | Contains |
|---|---|
| `~/R_libs` | Normal packages, including the original measles |
| `~/R_libs_crn` | Fair-pairs measles (`Rscript install_crn.R`) |
| `~/R_libs_dur` | Patched measles with the new durations (installed by `run_dur.sh`) |

Patched package: put `measlesmixingriskquarantine.hpp` in the folder you submit
from. `run_dur.sh` copies it into `~/crn_build/measles_crn/inst/include/measles/`,
compiles it into `~/R_libs_dur`, and checks that the compiled package contains
the new code. It reinstalls only when the `.hpp` changes, and two jobs never
install at the same time.

---

## Running

### Original model
```bash
quarto render 01_tiered_quarantine_tables.qmd
N_SIMS=200 quarto render 01_tiered_quarantine_tables.qmd    # quick test
```
(Load the same GCC module that built `~/R_libs`; see Troubleshooting.)

### Patched model, on Slurm
```bash
unset SCHOOL_SIZE RASH_DAYS N_RUNS N_SIMS R0
QMD="01_tiered_quarantine_tables_dur.qmd 02_outbreak_size_distributions_dur.qmd 03_quarantine_policy_analysis_dur.qmd" sbatch run_dur.sh
```

Options (put them in front of `sbatch`):

| Variable | Effect | Results folder |
|---|---|---|
| `QMD="a.qmd b.qmd"` | Files to render, in order (default `09_model_behaviour.qmd`) | |
| `N_SIMS=200` | Runs per scenario (default 4000) | |
| `R0=9` | Different R0; transmission rates recomputed | `results_dur_R0_9/` |
| `RASH_DAYS=4` | Unnoticed rash contagious for 4 days (β recalibrated) | `results_dur_rash/` |
| `SCHOOL_SIZE=400` | Every school 400 students, grade shares as in `ENROLLMENT` | `results_dur_size400/` |
| `MD=0` | HTML only (default: HTML + Markdown) | |

The HTML file names do not change with these options: rename them after each
job (e.g. add `_R0_9`) before running another setting.

### Sensitivity analysis
```bash
N_RUNS=200 QMD="sensitivity.qmd" sbatch run_dur.sh     # quick test
QMD="sensitivity.qmd" sbatch run_dur.sh                # 2,000 runs, about 1.5–2.5 h on 4 CPUs
```
Each parameter is set to a low and a high value with everything else at base,
at 60%, 70% and 80% coverage, with the transmission rate recalibrated whenever
contacts, school size or the contagious period change. Ranges are ±50% of base,
except vaccine efficacy 93–100%, quarantine compliance 90% and contact tracing
success 90%. Change `focus:` (strategy) and `outcome:` (`own` risk or `diff`
vs 21/21/21) at the top of the file. Output in `results_sens/`.

---

## Output and caching

| Where | What |
|---|---|
| `<file>.html` | Report with tables and figures |
| `<file>.md` + `md_figures/` | Markdown version (keep the two together) |
| `results*/` | Cached simulations, `<name>_n<N_SIMS>.csv` |
| `figures*/` | PNG copy of every figure |

Simulations are slow, so `run_or_load()` saves every set of runs. Rendering
again reuses them; **delete the cached file to re-run** after changing
settings that are not part of its name.

---

## Troubleshooting

| Message | Cause and fix |
|---|---|
| `GLIBCXX_3.4.32 not found` | A package from `~/R_libs` was built with a newer GCC than the job has. Render the `_dur` file through `run_dur.sh`, or `module load gcc/13.x` before rendering the original files. |
| `No space left on device` | Home folder full. `df -h ~`; delete old `dur08-*.out/.err`, `*_files/`, test caches (`results*/*_n200_*`). |
| `there is no package called 'measles'` during install | Two jobs installed at once, or the disk is full. `scancel` the extra jobs, `rm -rf ~/R_libs_dur/00LOCK-* ~/R_libs_dur/measles`, submit one job. |
| A change has no effect | Cached results were reused. Delete the matching file in `results*/`. |
| Wrong number of runs | A variable is still exported from an earlier test: `unset N_RUNS N_SIMS SCHOOL_SIZE RASH_DAYS R0`. |
| `make: No rule to make target` | Old instructions. Nothing here uses `make`; submit `run_dur.sh` with `QMD=...`. |

---

## References

- Grantz KH et al. Age-specific social mixing of school-aged children in a US setting using proximity detecting sensors and contact surveys. *Sci Rep* 2021.
- Guclu H et al. Social contact networks and mixing among students in K-12 schools in Pittsburgh, PA. *PLoS ONE* 2016.
- Lessler J et al. Incubation periods of acute respiratory viral infections: a systematic review. *Lancet Infect Dis* 2009.
- Kremer C et al. Informing quarantine policy for measles control in primary school and daycare settings, Flanders, Belgium. *Eurosurveillance* 2026;31(22).
- Vega Yon G et al. measles: Measles Epidemiological Models (R package). github.com/UofUEpiBio/measles.
