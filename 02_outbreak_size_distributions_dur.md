# 2. Outbreak-size distributions by quarantine strategy

How big do outbreaks get under each quarantine strategy? For every
vaccination level and one of 11 strategies (**high / medium / low-risk
quarantine days**), we look at the full distribution of final outbreak
sizes rather than a single probability.

- **Histograms** are drawn as outlines, one coloured line per school, so
  the three schools can be compared in the same panel. There is one
  panel per strategy and one graph per vaccination level.
- **Tables** give the probability of an outbreak of at least 10, 20 and
  50 cases, for each school and strategy.

<div>

> **Note**
>
> Model code lives in `measles_model.R`. Simulations are cached in
> `results/`, and PNG copies of every figure are saved to `figures/`.

</div>

## Histograms, all schools

::: {.panel-tabset}

``` r
for (v in VACCS) {
  cat("\n\n## ", vacc_lab(v), " vaccinated\n\n", sep = "")

  p <- plot_hist(v)

  print(save_fig(
    p,
    sprintf("fig_hist_all_schools_%02d.png", round(100 * v)),
    11, 10
  ))

  cat("\n\n")

}
```

## 50% vaccinated

## 60% vaccinated

![](md_figures/02_outbreak_size_distributions_dur-histograms-1.png)

## 70% vaccinated

![](md_figures/02_outbreak_size_distributions_dur-histograms-2.png)

## 80% vaccinated

![](md_figures/02_outbreak_size_distributions_dur-histograms-3.png)

## 85% vaccinated

![](md_figures/02_outbreak_size_distributions_dur-histograms-4.png)![](md_figures/02_outbreak_size_distributions_dur-histograms-5.png)
