logo:
	R -e 'source("logo.R")'

website:
	rm -rf docs/ man/; mkdir -p man/figures
	cp vignettes/figures/shiny_graph.png man/figures/
	cp vignettes/figures/shiny_maps.png man/figures/
	R -e 'devtools::document(); knitr::knit("README.Rmd", output = "README.md"); options(pkgdown.internet = FALSE); pkgdown::build_site();'
