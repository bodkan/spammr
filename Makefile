website:
	rm -rf docs/ man/; mkdir -p man/figures
	cp vignettes/figures/shiny_graph.png man/figures/
	cp vignettes/figures/shiny_maps.png man/figures/
	R -e 'devtools::document(); source("logo.R"); knitr::knit("README.Rmd", output = "README.md"); options(pkgdown.internet = FALSE); pkgdown::build_site();'
