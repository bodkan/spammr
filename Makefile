web:
	R -e 'source("logo.R"); knitr::knit("README.Rmd", output = "README.md"); devtools::document(); options(pkgdown.internet = FALSE); pkgdown::build_site();'
