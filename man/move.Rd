% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{move}
\alias{move}
\title{Move the population to a new location in a given amount of time}
\usage{
move(
  pop,
  trajectory,
  end,
  start,
  overlap = 0.8,
  snapshots = NULL,
  verbose = TRUE
)
}
\arguments{
\item{pop}{Object of the class \code{slendr_pop}}

\item{trajectory}{List of two-dimensional vectors (longitude, latitude)
specifying the migration trajectory}

\item{start, end}{Start/end points of the population migration}

\item{overlap}{Minimum overlap between subsequent spatial boundaries}

\item{snapshots}{The number of intermediate snapshots (overrides the
\code{overlap} parameter)}

\item{verbose}{Show the progress of searching through the number of
sufficient snapshots?}
}
\value{
Object of the class \code{slendr_pop}
}
\description{
This function defines a displacement of a population along a given trajectory
in a given time frame
}
