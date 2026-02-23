vax_rates <- c(0.95, 0.9, .8, .5)

for (vax_rate in vax_rates) {

  scenario_name <- sprintf("scenario-vax-%.2f", vax_rate)
  fn_md <- sprintf("scenarios/%s.md", scenario_name)
  fn_qmd <- sprintf("scenarios/%s.qmd", scenario_name)

  message("Processing scenario ", vax_rate)

  file.copy("template.qmd", fn_qmd, overwrite = TRUE)
  
  quarto::quarto_render(
    input = fn_qmd,
    output_file = paste0(scenario_name, ".md"),
    execute_params = list(
      vax_rate = vax_rate,
      n_sims   = 10000
      )
  )

}


