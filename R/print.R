#' Print a short summary of a \code{slendr_pop} object
#'
#' All spatial objects in the slendr package are internally
#' represented by a Simple Features (sf) object. This fact is hidden
#' in most circumstances this, as the goal of the slendr package is to
#' provide functionality at a much higher level (population
#' boundaries, geographic regions, instead of individual polygons and
#' other "low-level" geometric objects). The logical argument
#' \code{sf} allows the user to inspect the underlying object
#' structure. Similarly, by default, \code{print} does not print all
#' spatial map information to reduce clutter. If many individual
#' spatial maps for a population are present, a full table of spatial
#' snapshots can be printed by typing \code{x[]}.
#'
#' @param x Object of a class \code{slendr_pop}
#'
#' @export
print.slendr_pop <- function(x, ...) {
  print_header_info(x)

  cat("name:", unique(x$pop), "\n")

  # removal time of the population
  if (attr(x, "remove") == -1)
    cat("stays until the end of the simulation\n")
  else
    cat("scheduled removal at time ", attr(x, "remove"), "\n")

  # aquatic or terrestrial?
  if (attr(x, "aquatic") == FALSE)
    cat("habitat: terrestrial\n")
  else
    cat("habitat: aquatic\n")

  print_pop_history(x)

  cat("\nnumber of spatial maps:", nrow(x), "\n")

  print_map_info(x)
}


#' Print a short summary of a \code{slendr_region} object
#'
#' All spatial objects in the slendr package are internally
#' represented by a Simple Features (sf) object. This fact is hidden
#' in most circumstances this, as the goal of the slendr package is to
#' provide functionality at a much higher level (population
#' boundaries, geographic regions, instead of individual polygons and
#' other "low-level" geometric objects). The logical argument
#' \code{sf} allows the user to inspect the underlying object
#' structure. Similarly, by default, \code{print} does not print all
#' spatial map information to reduce clutter. If many individual
#' spatial maps for a population are present, a full table of spatial
#' snapshots can be printed by typing \code{x[]}.
#'
#' @param x Object of a class \code{slendr_region}
#'
#' @export
print.slendr_region <- function(x, ...) {
  print_header_info(x)

  if (nrow(x))
    cat("name:", x$region, "\n\n")
  else
    cat("[this region object is empty]\n\n")

  print_map_info(x)
}


#' Print a short summary of a \code{slendr_map} object
#'
#' All spatial objects in the slendr package are internally
#' represented by a Simple Features (sf) object. This fact is hidden
#' in most circumstances this, as the goal of the slendr package is to
#' provide functionality at a much higher level (population
#' boundaries, geographic regions, instead of individual polygons and
#' other "low-level" geometric objects). The logical argument
#' \code{sf} allows the user to inspect the underlying object
#' structure. Similarly, by default, \code{print} does not print all
#' spatial map information to reduce clutter. If many individual
#' spatial maps for a population are present, a full table of spatial
#' snapshots can be printed by typing \code{x[]}.
#'
#' @param x Object of a class \code{slendr_map}
#'
#' @export
print.slendr_map <- function(x, ...) {
  print_header_info(x)
  print_map_info(x)
}


#' Print a short summary of a \code{slendr_model} object
#'
#' All spatial objects in the slendr package are internally
#' represented by a Simple Features (sf) object. This fact is hidden
#' in most circumstances this, as the goal of the slendr package is to
#' provide functionality at a much higher level (population
#' boundaries, geographic regions, instead of individual polygons and
#' other "low-level" geometric objects). The logical argument
#' \code{sf} allows the user to inspect the underlying object
#' structure. Similarly, by default, \code{print} does not print all
#' spatial map information to reduce clutter. If many individual
#' spatial maps for a population are present, a full table of spatial
#' snapshots can be printed by typing \code{x[]}.
#'
#' @param x Object of a class \code{slendr_model}
#'
#' @export
print.slendr_model <- function(x, ...) {
  print_header_info(x)

  cat("populations:", paste0(x$splits$pop, collapse = ", "), "\n")
  cat("geneflow events: ")
  if (!is.null(x$geneflows))
    cat(nrow(x$geneflows), "\n")
  else
    cat("[no geneflow]\n")
  cat("generation time:", x$generation_time, "\n")
  cat("time direction:", x$direction, "\n")
  cat("number of spatial maps:", nrow(x$maps), "\n")
  cat("resolution:", x$resolution, "distance unit per pixel\n\n")
  cat("configuration files in:", normalizePath(x$directory), "\n\n")
  cat(
    "A detailed model specification can be found in `$splits`, `$geneflows`,
`$maps`, `$populations`, and other components of the model object (for
a complete list see `names(<model object>)`). You can also examine
the serialized configuration files in the model directory.\n")
}


get_slendr_type <- function(x) {
  grep("slendr_", class(x), value = TRUE) %>%
    gsub("slendr_", "", .) %>%
    gsub("pop", "population", .)
}


print_header_info <- function(x) {
  type <- get_slendr_type(x)

  header <- sprintf("slendr '%s' object", type)
  sep <- paste(rep("-", nchar(header)), collapse = "")

  cat(header, "\n")
  cat(sep, "\n")
}


print_map_info <- function(x) {
  type <- get_slendr_type(x)

  cat("map: ")
  if (type == "map" | !is.null(attr(x, "map"))) {
    crs <- sf::st_crs(x)$epsg
    if (is.na(crs)) {
      cat("abstract spatial landscape ")
      if (nrow(x))
        cat("with custom features\n")
      else
        cat("with no features\n")
      units <- ""
    } else {
      crs <- paste("EPSG", crs)
      cat("internal coordinate reference system:", crs, "\n")
      units <- " (in degrees longitude and latitude)"
    }

    xrange <- attr(x, "xrange")
    yrange <- attr(x, "yrange")
    cat(sprintf("spatial limits%s:\n  - vertical %d ... %d\n  - horizontal %d ... %d\n",
                units, xrange[1], xrange[2], yrange[1], yrange[2]))
  } else
    cat("[no map defined]\n")
}

print_pop_history <- function(x) {
  cat("\npopulation history overview:\n")
  history <- attr(x, "history")
  for (event in history) {
    cat("  - time ")

    # population split
    if (event$event == "split") {
      cat(event$time, ": ", sep = "")
      parent <- attr(x, "parent")
      if (is.character(parent) && parent == "ancestor")
        cat("created as an ancestral population")
      else {
        cat("split from", parent$pop[1])
      }
    }

    # spatial dynamics events
    else if (event$event == "move") {
      cat(sprintf("%d-%d: movement across a landscape", event$tstart, event$tend))
    } else if (event$event == "expand") {
      cat(sprintf("%d-%d: range expansion", event$tstart, event$tend))
    } else if (event$event == "range") {
      cat(sprintf("%d: change of the spatial boundary", event$time))
    }

    # population size change
    else if (event$event == "resize" && event$how == "step") {
      cat(sprintf("%d: resize from %d to %d individuals",
                  event$tresize, event$prev_N, event$N))
    } else if (event$event == "resize" && event$how == "exponential") {
      cat(sprintf("%d-%d: exponential resize from %d to %d individuals",
                  event$tresize, event$tend, event$prev_N, event$N))
    }

    # change of dispersal parameters
    else if (event$event == "dispersal") {
      cat(sprintf("%d: change in spatial interaction", event$tdispersal))
      if (!is.na(event$competition_dist)) cat("\n        - competition distance:", event$competition_dist)
      if (!is.na(event$mate_dist)) cat("\n        - mate choice distance:", event$competition_dist)
      if (!is.na(event$dispersal_dist)) cat("\n        - dispersal from parent:", event$competition_dist)
    } else
      stop("Unknown event type", call. = FALSE)

    cat("\n")
  }
}
