#' Runs the simulation with specified parameters.
#' @param duration A vector of length 3 indicating the number of quarantine days for each risk level.
#' @param vaccinated Proportion of vaccinated individuals in the population.
#' @return
#' A data.table containing the total history of the simulation.
simulator <- function(
  duration = c(21L, 21L, 21L),
  vaccinated = 0.9
) {
  # Building the model
  model_baseline <- ModelMeaslesMixingRiskQuarantine(
    n = n_agents,
    prevalence = 1 / 600,
    contact_matrix = contact_matrix,
    transmission_rate = p_infect,
    prop_vaccinated = vaccinated,
    detection_rate_quarantine = 0.0,
    contact_tracing_days_prior = 7,
    quarantine_period_high = duration[1],
    quarantine_period_medium = duration[2],
    quarantine_period_low = duration[3]
  )

  # Creating entities
  for (i in 1:n_classes) {
    model_baseline |>
      add_entity(
        entity("Class", n_agents_per_class, as_proportion = FALSE)
      )
  }

  # Running the simulation multiple times
  model_baseline |>
    run_multiple(
      ndays = n_days,
      nsims = n_sims,
      seed = 221,
      saver = make_saver("total_hist"),
      nthreads = n_threads,
      verbose = interactive()
    )

  # Getting the results
  ans_baseline <- model_baseline |>
    run_multiple_get_results(freader = data.table::fread, nthreads = 1L)

  # Returning the total result
  ans <- ans_baseline$total_hist

  # Extracting the final counts
  ans[date == max(date),][
    (state != "Susceptible") &
          (state != "Susceptible Quarantine"),
    .(
      total_infected = sum(counts)
    ),
    by = .(sim_num)
  ]
}

#' Summarizes the simulation results into a table.
#' @param ans A list of simulation results.
#' @return
#' A formatted table summarizing the outbreak sizes and confidence intervals.
tabulator <- function(ans) {
  scenario_names <- names(ans)
    data.table(
      Scenario = scenario_names,
      `Mean` = sapply(
        scenario_names,
        function(x) sprintf(
          "%.2f",
          mean(ans[[x]]$total_infected)
        )
      ) ,
      `Median` = sapply(
        scenario_names,
        function(x) sprintf(
          "%.2f",
          median(ans[[x]]$total_infected)
        )
      ),
      `95% CI` = sapply(
        scenario_names,
        function(x) {
          ci <- quantile(
            ans[[x]]$total_infected,
            probs = c(0.025, 0.975)
          )
          sprintf("(% 5.2f, % 5.2f)", ci[1], ci[2])
        }
      )
    ) |> knitr::kable(caption = "Outbreak sizes across different quarantine scenarios.")
}