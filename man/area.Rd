% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{area}
\alias{area}
\title{Calculate the area covered by the given slendr object}
\usage{
area(x)
}
\arguments{
\item{x}{Object of the class \code{slendr}}
}
\value{
Area covered by the input object. If the coordinate reference system
was specified, the area in projected units (i.e. m^2) is returned.
Otherwise the function returns area without units.
}
\description{
Calculate the area covered by the given slendr object
}
