README.md: README.qmd
	quarto render README.qmd

scenarios: template.qmd
	R CMD BATCH --vanilla scenarios.R scenarios.Rout &

clean:
	rm scenarios/*.rda

.PHONY: scenarios
