% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{boundary}
\alias{boundary}
\title{Update the population range}
\usage{
boundary(
  pop,
  time,
  center = NULL,
  radius = NULL,
  polygon = NULL,
  overlap = 0.8
)
}
\arguments{
\item{pop}{Object of the class \code{slendr_pop}}

\item{time}{Time of the change}

\item{center}{Two-dimensional vector specifying the center of the circular range}

\item{radius}{Radius of the circular range}

\item{polygon}{List of vector pairs, defining corners of the polygon range (see also the
\code{region} argument) or a geographic region of the class \code{slendr_region} from which the
polygon coordinates will be extracted}

\item{overlap}{Minimum required overlap with the previous active population boundary (set to 0 to
disable the check)}
}
\value{
Object of the class \code{slendr_pop}
}
\description{
This function allows a more manual control of spatial map changes in addition to the
\code{expand} and \code{move} functions
}
