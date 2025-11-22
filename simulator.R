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
#' @param sizes A numeric vector of outbreak size thresholds to compute probabilities for. Default is c(10, 20, 50).
#' @return
#' A formatted table summarizing the probability of outbreaks of various sizes.
tabulator <- function(ans, sizes = c(10, 20, 50)) {
  if (length(sizes) == 0 || !is.numeric(sizes) || any(sizes <= 0)) {
    stop("'sizes' must be a non-empty numeric vector with positive values")
  }
  if (is.null(ans) || length(ans) == 0 || is.null(names(ans))) {
    stop("'ans' must be a non-empty named list of simulation results")
  }
  scenario_names <- names(ans)
  
  # Create the base data.table with scenario names
  result <- data.table(Scenario = scenario_names)
  
  # Add a column for each size threshold
  for (size in sizes) {
    col_name <- paste0("P(≥", size, ")")
    result[[col_name]] <- sapply(
      scenario_names,
      function(x) sprintf(
        "%.3f",
        mean(ans[[x]]$total_infected >= size)
      )
    )
  }
  
  result |> knitr::kable(caption = "Probability of outbreak sizes across different quarantine scenarios.")
}