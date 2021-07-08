#' Define a population and its original spatial range
#'
#' Defines a spatial range of a population and its most important parameters.
#'
#' There are three ways to specify a spatial boundary: i) circular range
#' specified using a center coordinate and a radisu, ii) polygon specified as a
#' list of two-dimensional vector coordinates, iii) polygon as in ii), but
#' defined (and named) using the \code{region} function.
#'
#' @param name Name of the population
#' @param time Time of the population's first appearance
#' @param N Number of individuals at the time of first appearance
#' @param parent Parent population object or "ancestor" (indicating that the
#'   population does not have an ancestor, and that it is the first population
#'   in its "lineage")
#' @param map Object of the type \code{slendr_map} which defines the world
#'   context (created using the \code{world} function)
#' @param center Two-dimensional vector specifying the center of the circular
#'   range
#' @param radius Radius of the circular range
#' @param polygon List of vector pairs, defining corners of the polygon range or
#'   a geographic region of the class \code{slendr_region} from which the
#'   polygon coordinates will be extracted (see the \code{region() function})
#' @param remove Time at which the population should be removed
#' @param intersect Intersect the population's boundaries with landscape
#'   features?
#' @param competition_dist,mate_dist Maximum spatial competition and mating
#'   choice distance
#' @param dispersal_dist Standard deviation of the normal distribution of the
#'   distance that offspring disperses from its parent
#' @param aquatic Is the species aquatic (\code{FALSE} by default, i.e.
#'   terrestrial species)?
#'
#' @return Object of the class \code{slendr_pop}
#'
#' @export
population <- function(name, time, N, parent = "ancestor", map = NULL,
                       center = NULL, radius = NULL, polygon = NULL,
                       remove = NULL, intersect = TRUE,
                       competition_dist = NA, mate_dist = NA, dispersal_dist = NA,
                       aquatic = FALSE) {
  # is this the first population defined in the model?
  if (is.character(parent) && parent == "ancestor") {
    if (is.null(map))
      stop("Ancestral population must specify its 'map'", call. = FALSE)
    else
      map <- map
  } else {
    check_split_time(time, parent)
    map <- attr(parent, "map")
  }

  # define the population range as a simple geometry object
  # and bind it with the annotation info into an sf object
  if (!is.null(polygon) & inherits(polygon, "slendr_region"))
    geometry <- sf::st_geometry(polygon)
  else
    geometry <- define_boundary(map, center, radius, polygon)

  pop <- sf::st_sf(
    data.frame(pop = name, time = time, stringsAsFactors = FALSE),
    geometry = geometry
  )
  sf::st_agr(pop) <- "constant"

  # when to clean up the population?
  attr(pop, "remove") <- if (!is.null(remove)) remove else -1

  # keep a record of the parent population
  if (inherits(parent, "slendr_pop")) {
    attr(pop, "parent") <- parent
    # keep the map as an internal attribute
    attr(pop, "map") <- map
  } else if (is.character(parent) & parent == "ancestor") {
    attr(pop, "parent") <- "ancestor"
    attr(pop, "map") <- map
  } else
    stop("Suspicious parent population", call. = FALSE)

  attr(pop, "intersect") <- intersect
  attr(pop, "aquatic") <- aquatic

  # create the first population history event - population split
  attr(pop, "history") <- list(data.frame(
    pop =  name,
    event = "split",
    time = time,
    N = N,
    competition_dist = competition_dist,
    mate_dist = mate_dist,
    dispersal_dist = dispersal_dist
  ))

  class(pop) <- set_class(pop, "pop")

  pop
}


#' Move the population to a new location in a given amount of time
#'
#' This function defines a displacement of a population along a given trajectory
#' in a given time frame
#'
#' @param pop Object of the class \code{slendr_pop}
#' @param trajectory List of two-dimensional vectors (longitude, latitude)
#'   specifying the migration trajectory
#' @param start,end Start/end points of the population migration
#' @param overlap Minimum overlap between subsequent spatial boundaries
#' @param snapshots The number of intermediate snapshots (overrides the
#'   \code{overlap} parameter)
#' @param verbose Show the progress of searching through the number of
#'   sufficient snapshots?
#'
#' @return Object of the class \code{slendr_pop}
#'
#' @export
move <- function(pop, trajectory, end, start, overlap = 0.8, snapshots = NULL,
                 verbose = TRUE) {
  check_event_time(c(start, end), pop)
  check_removal_time(start, pop)

  if (!is.null(snapshots))
    if (snapshots <= 0)
      stop("The number of snapshots must be a non-negative integer", call. = FALSE)

  if (!(overlap > 0 & overlap < 1))
    stop("The required overlap between subsequent spatial maps must be a number between 0 and 1",
         call. = FALSE)

  map <- attr(pop, "map")

  # take care of just a single destination point being specified
  if (!is.list(trajectory) & length(trajectory) == 2)
    trajectory <- list(trajectory)

  # get the last available population boundary
  region_start <- pop[nrow(pop), ]
  region_start$time <- start
  sf::st_agr(region_start) <- "constant"

  # prepend the coordinates of the first region to the list of "checkpoints"
  # along the way of the movement
  start_coords <- sf::st_centroid(region_start)

  if (has_crs(map)) {
    source_crs <- "EPSG:4326"
    target_crs <- sf::st_crs(pop)
    start_coords <- sf::st_transform(start_coords, crs = source_crs)
  }

  start_coords <- sf::st_coordinates(start_coords)
  checkpoints <- c(list(as.vector(start_coords)), trajectory)

  traj <- sf::st_linestring(do.call(rbind, checkpoints)) %>% sf::st_sfc()

  if (has_crs(map)) {
    traj <- sf::st_sf(traj, crs = source_crs) %>%
      sf::st_transform(crs = target_crs)
  } else {
    traj <- sf::st_sf(traj)
  }

  if (is.null(snapshots)) {
    n <- 1
    message("Iterative search for the minimum sufficient number of intermediate ",
            "spatial snapshots, starting at ", n, ". This should only take a couple of ",
            "seconds, but if you don't want to wait, you can set `snapshots = N` manually.")
  } else
    n <- snapshots

  # iterate through the number of intermediate spatial boundaries to reach
  # the required overlap between subsequent spatial maps
  repeat {
    traj_segments <- sf::st_segmentize(traj, sf::st_length(traj) / n)

    traj_points <- sf::st_cast(traj_segments, "POINT")
    traj_points_coords <- sf::st_coordinates(traj_points)
    traj_diffs <- diff(traj_points_coords)

    time_slices <- seq(start, end, length.out = nrow(traj_points))[-1]
    traj_diffs <- cbind(traj_diffs, time = time_slices)

    inter_regions <- list()
    inter_regions[[1]] <- region_start
    for (i in seq_len(nrow(traj_diffs))) {
      shifted_region <- sf::st_geometry(inter_regions[[i]]) + traj_diffs[i, c("X", "Y")]
      inter_regions[[i + 1]] <- sf::st_sf(
        data.frame(
          pop = region_start$pop,
          time = traj_diffs[i, "time"],
          stringsAsFactors = FALSE
        ),
        geometry = shifted_region,
        crs = sf::st_crs(inter_regions[[i]])
      )
      sf::st_agr(inter_regions[[i + 1]]) <- "constant"
    }

    if (!is.null(snapshots)) break

    # calculate the overlap between subsequent spatial snapshots
    overlaps <- compute_overlaps(do.call(rbind, inter_regions))

    if (all(overlaps >= overlap)) {
      message("The required ", sprintf("%.1f%%", 100 * overlap),
              " overlap between subsequent spatial maps has been met")
      break
    } else {
      n <- n + 1
      if (verbose)
        message("- testing ", n, " snapshots")
    }
  }

  inter_regions <- rbind(pop, do.call(rbind, inter_regions))
  sf::st_agr(inter_regions) <- "constant"

  result <- copy_attributes(
    inter_regions, pop,
    c("map", "parent", "remove", "intersect", "aquatic", "history")
  )

  attr(result, "history") <- append(attr(result, "history"), list(data.frame(
    pop =  unique(region_start$pop),
    event = "move",
    start = start,
    end = end
  )))

  result
}


#' Expand the population range
#'
#' Expands the spatial population range by a specified distance in a given
#' time-window
#'
#' @param pop Object of the class \code{slendr_pop}
#' @param by How many units of distance to expand by?
#' @param start,end When does the expansion start/end?
#' @param overlap Minimum overlap between subsequent spatial boundaries
#' @param snapshots The number of intermediate snapshots (overrides the
#'   \code{overlap} parameter)
#' @param polygon Geographic region to restrict the expansion to
#'
#' @return Object of the class \code{slendr_pop}
#'
#' @export
expand <- function(pop, by, end, start, overlap = 0.8, snapshots = NULL,
                   polygon = NULL, verbose = TRUE) {
  shrink_or_expand(pop, by, end, start, overlap, snapshots, polygon, verbose)
}


#' Shrink the population range
#'
#' Shrinks the spatial population range by a specified distance in a given
#' time-window
#'
#' @param pop Object of the class \code{slendr_pop}
#' @param by How many units of distance to shrink by?
#' @param start,end When does the boundary shrinking start/end?
#' @param overlap Minimum overlap between subsequent spatial boundaries
#' @param snapshots The number of intermediate snapshots (overrides the
#'   \code{overlap} parameter)
#'
#' @return Object of the class \code{slendr_pop}
#'
#' @export
shrink <- function(pop, by, end, start, overlap = 0.8, snapshots = NULL, verbose = TRUE) {
  shrink_or_expand(pop, -by, end, start, overlap, snapshots, polygon = NULL, verbose)
}


#' Update the population range
#'
#' This function allows a more manual control of spatial map changes in addition to the
#' \code{expand} and \code{move} functions
#'
#' @param pop Object of the class \code{slendr_pop}
#' @param time Time of the change
#' @param center Two-dimensional vector specifying the center of the circular range
#' @param radius Radius of the circular range
#' @param polygon List of vector pairs, defining corners of the polygon range (see also the
#'   \code{region} argument) or a geographic region of the class \code{slendr_region} from which the
#'   polygon coordinates will be extracted
#' @param overlap Minimum required overlap with the previous active population boundary (set to 0 to
#'   disable the check)
#'
#' @return Object of the class \code{slendr_pop}
#'
#' @export
boundary <- function(pop, time, center = NULL, radius = NULL,
                     polygon = NULL, overlap = 0.8) {
  check_event_time(time, pop)
  check_removal_time(time, pop)

  map <- attr(pop, "map")

  # define the population range as a simple geometry object
  # and bind it with the annotation info into an sf object
  if (!is.null(polygon) & inherits(polygon, "slendr_region"))
    geometry <- sf::st_geometry(polygon)
  else
    geometry <- define_boundary(map, center, radius, polygon)

  updated <- sf::st_sf(
    data.frame(pop = unique(pop$pop), time = time, stringsAsFactors = FALSE),
    geometry = geometry
  )

  if (compute_overlaps(rbind(updated, pop[nrow(pop), ])) < overlap)
    stop("Insufficient overlap with the last active spatial boundary (",
         "please adjust the new spatial boundary or adjust the `overlap = `",
         "parameter for less stringent checking)", call. = FALSE)

  combined <- rbind(pop, updated)
  sf::st_agr(combined) <- "constant"

  result <- copy_attributes(
    combined, pop,
    c("map", "parent", "remove", "intersect", "aquatic", "history")
  )

  attr(result, "history") <- append(attr(result, "history"), list(data.frame(
    event = "range",
    time = time
  )))

  result
}


#' Resize the population size
#'
#' @param pop Object of the class \code{slendr_pop}
#' @param N Population size after the change
#' @param how How to change the population size (options are \code{"step"},
#'   \code{"exponential"}, \code{"linear"})
#' @param time Time of the population size change
#' @param start,end Time-window for the population size change (used for
#'   exponential or linear change)
#'
#' @export
resize <- function(pop, N, how = "step", time = NULL, start = NULL, end = NULL) {
  if (N < 1) stop("Only positive population sizes allowed", call. = FALSE)

  if (!how %in% c("step", "exponential", "linear"))
    stop("Only 'step', 'exponential' and 'linear' are allowed as arguments
for the 'how' parameter", call. = FALSE)

  if ((how == "step" & is.null(time)) |
      (how %in% c("exponential", "linear") & is.null(c(start, end))))
    stop("Timing of the population change needs to be specified", call. = FALSE)

  # get the last active population size
  prev_N <- sapply(attr(pop, "history"), function(event) event$N) %>%
    Filter(Negate(is.null), .) %>%
    unlist %>%
    tail(1)

  change <- data.frame(
    pop =  unique(pop$pop),
    event = "resize",
    resize_how = how,
    N = N,
    prev_N = prev_N
  )

  if (how == "step") {
    check_event_time(time, pop)
    check_removal_time(time, pop)
    change$time <- time
  } else {
    check_event_time(c(start, end), pop)
    check_removal_time(start, pop)
    change$tstart <- start
    change$tend <- end
  }

  attr(pop, "history") <- append(attr(pop, "history"), list(change))

  pop
}


#' Define a geneflow event
#'
#' @param from,to Objects of the class \code{slendr_pop}
#' @param rate Scalar value in the range (0, 1] specifying the
#'   proportion of migration over given time period
#' @param start,end Start and end of the geneflow event
#' @param overlap Require spatial overlap between admixing
#'   populations?  (default \code{TRUE})
#'
#' @return Object of the class data.frame
#'
#' @export
geneflow <- function(from, to, rate, start, end, overlap = TRUE) {
  # make sure the population is not removed during the the admixture period
  check_removal_time(start, from)
  check_removal_time(end, from)
  check_removal_time(start, to)
  check_removal_time(end, to)

  from_name <- unique(from$pop)
  to_name <- unique(to$pop)

  if (from$time[1] <= start & from$time[1] <= end &
      to$time[1] <= start & to$time[1] <= end)
    comp <- `<=`
  else if (from$time[1] >= start & from$time[1] >= end &
           to$time[1] >= start & to$time[1] >= end)
    comp <- `>=`
  else
    stop(sprintf("Specified times are not consistent with the assumed direction of time (geneflow %s -> %s in the time window %s-%s)",
                 from_name, to_name, start, end),
         call. = FALSE)

  # get the last specified spatial maps before the geneflow time
  region_from <- intersect_features(from[comp(from$time, start), ] %>% .[nrow(.), ])
  region_to <- intersect_features(to[comp(to$time, start), ] %>% .[nrow(.), ])

  if (nrow(region_from) == 0)
    stop(sprintf("No spatial map defined for %s at/before the time %d",
                 from_name, start),
         call. = FALSE)
  if (nrow(region_to) == 0)
    stop(sprintf("No spatial map defined for %s at/before the time %d",
                 to_name, start),
         call. = FALSE)

  # calculate the overlap of spatial ranges between source and target
  region_overlap <- sf::st_intersection(region_from, region_to)
  area_overlap <- as.numeric(sum(sf::st_area(region_overlap)))

  if (overlap & area_overlap == 0) {
    stop(sprintf("No overlap between population ranges of %s and %s at time %d.

Please check the spatial maps of both populations by running
`plot(%s, %s)` and adjust them accordingly. Alternatively, in case
this makes sense for your model, you can add `overlap = F` which
will instruct slendr to simulate geneflow without spatial overlap
between populations.",
from_name, to_name, start, deparse(substitute(from)),
deparse(substitute(to))), call. = FALSE)
  }

  data.frame(
    from_name = from_name,
    to_name = to_name,
    tstart = start,
    tend = end,
    rate = rate,
    overlap = overlap,
    stringsAsFactors = FALSE
  )
}


#' Define a world map for all spatial operations
#'
#' Defines either an abstract geographic landscape (blank or
#' containing user-defined landscape) or using a real Earth
#' cartographic data from the Natural Earth project
#' (<https://www.naturalearthdata.com>).
#'
#' @param xrange Two-dimensional vector specifying minimum and maximum
#'   horizontal range ("longitude" if using real Earth cartographic
#'   data)
#' @param yrange Two-dimensional vector specifying minimum and maximum
#'   vertical range ("latitude" if using real Earth cartographic data)
#' @param landscape Either "blank" (for blank abstract geography),
#'   "naturalearth" (for real Earth geography) or an object of the
#'   class \code{sf} defining abstract geographic features of the
#'   world
#' @param crs EPSG code of a coordinate reference system to use for
#'   spatial operations. No CRS is assumed by default (\code{NULL}),
#'   implying an abstract landscape not tied to any real-world
#'   geographic region (when \code{landscape = "blank"} or when
#'   \code{landscape} is a custom-defined geographic landscape), or
#'   implying WGS-84 (EPSG 4326) coordinate system when a real Earth
#'   landscape was defined (\code{landscape = "naturalearth"}).
#' @param ne_dir Path to the directory where Natural Earth data was
#'   manually downloaded and unzipped from
#'   <https://www.naturalearthdata.com/downloads/110m-physical-vectors>
#'   (used only when \code{landscape = "naturalearth"})
#'
#' @return Object of the class \code{slendr_map}
#'
#' @export
world <- function(xrange, yrange, landscape = "naturalearth", crs = NULL, ne_dir = NULL) {
  if (inherits(landscape, "sf")) { # a landscape defined by the user
    map <- sf::st_sf(landscape = sf::st_geometry(landscape)) %>%
      set_bbox(xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2])
  } else if (landscape == "blank") { # an empty abstract landscape
    map <- sf::st_sf(geometry = sf::st_sfc()) %>%
      set_bbox(xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2])
  } else if (landscape == "naturalearth") {  # Natural Earth data vector landscape
    scale <- "small"
    type <- "land"
    category <- "physical"

    if (is.null(ne_dir)) {
      map_raw <- rnaturalearth::ne_download(scale, type, category, returnclass = "sf")
    } else {
      map_raw <- rnaturalearth::ne_load(scale, type, category,
                                        returnclass = "sf", destdir = ne_dir)
    }
    sf::st_agr(map_raw) <- "constant"

    ## transform the map (default geographic CRS) into the target CRS
    map_transf <- sf::st_transform(map_raw, crs)

    ## define boundary coordinates in the target CRS
    zoom_bounds <- define_zoom(xrange, yrange, "EPSG:4326")
    zoom_transf <- sf::st_transform(zoom_bounds, crs)

    ## crop the map to the boundary coordinates
    map <- sf::st_crop(map_transf, zoom_transf)
  } else {
    stop("Landscape has to be either 'blank', 'naturalearth' or an object of the class 'sf'",
         call. = FALSE)
  }

  sf::st_agr(map) <- "constant"

  class(map) <- set_class(map, "map")
  attr(map, "xrange") <- xrange
  attr(map, "yrange") <- yrange

  map
}


#' Define a geographic region
#'
#' Creates a geographic region (a polygon) on a given map and gives it
#' a name. This can be used to define objects which can be reused in
#' multiple places in a slendr script (such as \code{region} arguments
#' of \code{population}) without having to repeatedly define polygon
#' coordinates.
#'
#' @param name Name of the geographic region
#' @param map Object of the type \code{sf} which defines the map
#' @param center Two-dimensional vector specifying the center of the
#'   circular range
#' @param radius Radius of the circular range
#' @param polygon List of vector pairs, defining corners of the
#'   polygon range or a geographic region of the class
#'   \code{slendr_region} from which the polygon coordinates will be
#'   extracted (see the \code{region() function})
#'
#' @return Object of the class \code{slendr_region}
#'
#' @export
region <- function(name = NULL, map = NULL, center = NULL, radius = NULL, polygon = NULL) {
  if (is.null(name)) name <- "unnamed region"
  region <- sf::st_sf(
    region = name,
    geometry = define_boundary(map, center, radius, polygon)
  )
  sf::st_agr(region) <- "constant"

  # keep the map as an internal attribute
  attr(region, "map") <- map

  class(region) <- set_class(region, "region")
  region
}


#' Convert between coordinate systems
#'
#' Converts between coordinates on a compiled raster map (i.e. pixel
#' units) and different Geographic Coordinate Systems (CRS).
#'
#' @param x,y Coordinates in two dimensions (if missing, coordinates
#'   are expected to be in the \code{data.frame} specified in the
#'   \code{coords} parameter as columns "x" and "y")
#' @param from,to Either a CRS code accepted by GDAL, a valid integer
#'   EPSG value, an object of class \code{crs}, the value "raster"
#'   (converting from/to pixel coordinates), or "world" (converting
#'   from/to whatever CRS is set for the underlying map)
#' @param coords data.frame-like object with coordinates in columns
#'   "x" and "y"
#' @param model Object of the class \code{slendr_model}
#' @param add Add column coordinates to the input data.frame
#'   \code{coords} (coordinates otherwise returned as a separate
#'   object)?
#'
#' @return Data.frame with converted two-dimensional coordinates
#'
#' @export
convert <- function(from, to, x = NULL, y = NULL, coords = NULL, model = NULL, add = FALSE) {
  if ((is.null(x) | is.null(y)) & is.null(coords))
    stop("Coordinates for conversion are missing", call. = FALSE)

  if ((from == "raster" | to == "raster") & is.null(model))
    stop("Model object needs to be specified for conversion of raster coordinates", call. = FALSE)

  if (add & is.null(coords))
    stop("Converted coordinates can only be added to a provided data.frame", call. = FALSE)

  if (!is.null(coords) & !all(c("x", "y") %in% colnames(coords)))
    stop("Columns 'x' and 'y' must be present in the input data.frame", call. = FALSE)

  if (!is.null(model)) {
    # dimension of the map in the projected CRS units
    bbox <- sf::st_bbox(model$world)
    map_dim <- c(bbox["xmax"] - bbox["xmin"], bbox["ymax"] - bbox["ymin"])

    # dimension of the rasterized map in pixel units
    # (x/y dimensions of PNGs are reversed)
    raster_dim <- dim(png::readPNG(model$maps$path[1]))[2:1]
  }

  if (to == "world") to <- sf::st_crs(model$world)
  if (from == "world") from <- sf::st_crs(model$world)

  if (is.null(coords))
    df <- data.frame(x = x, y = y)
  else
    df <- coords[, c("x", "y")]

  if (from == "raster") {
    # convert pixel coordinates to na sf object in world-based coordinates
    df$x <- bbox["xmin"] + map_dim[1] * df$x / raster_dim[1]
    df$y <- bbox["ymin"] + map_dim[2] * df$y / raster_dim[2]
    point <- sf::st_as_sf(df, coords = c("x", "y"), crs = sf::st_crs(model$world))
  } else {
    # ... otherwise create a formal sf point object from the
    # coordinates already given
    point <- sf::st_as_sf(df, coords = c("x", "y"), crs = from)
  }

  if (to == "raster") {
    point_coords <- sf::st_transform(point, crs = sf::st_crs(model$world)) %>%
      sf::st_coordinates()
    newx <- abs((point_coords[, "X"] - bbox["xmin"])) / map_dim[1] * raster_dim[1]
    newy <- abs((point_coords[, "Y"] - bbox["ymin"])) / map_dim[2] * raster_dim[2]
    new_point <- data.frame(newx = round(as.vector(newx)), newy = round(as.vector(newy)))
  } else {
    new_point <- sf::st_transform(point, crs = to) %>% sf::st_coordinates()
    colnames(new_point) <- c("newx", "newy")
  }

  if (nrow(new_point) == 1)
    return(as.vector(unlist(new_point)))

  if (add) new_point <- cbind(coords, new_point)

  new_point
}


#' Combine two \code{slendr_region} objects into a single geographic
#' region
#'
#' @param x Object of the class \code{slendr_region}
#' @param y Object of the class \code{slendr_region}
#' @param name Name of the resulting geographic region
#'
#' @return Object of the class \code{slendr_region}
#'
#' @export
join <- function(x, y, name = NULL) {
  if (!inherits(x, "slendr_region")) x <- region(polygon = x)
  if (!inherits(y, "slendr_region")) y <- region(polygon = y)

  result <- sf::st_union(x, y)
  result$region.1 <- NULL
  if (is.null(name))
    result$region <- sprintf("(%s and %s)", x$region, y$region)
  else
    result$region <- name
  attrs <- if (!is.null(attr(x, "map"))) "map" else NULL
  result <- copy_attributes(result, x, attrs)
  sf::st_agr(result) <- "constant"
  result
}


#' Generate the overlap of two \code{slendr_region} objects
#'
#' @inheritParams join
#'
#' @return Object of the class \code{slendr_region}
#'
#' @export
overlap <- function(x, y, name = NULL) {
  if (!inherits(x, "slendr_region")) x <- region(polygon = x)
  if (!inherits(y, "slendr_region")) y <- region(polygon = y)

  result <- sf::st_intersection(x, y)
  if (nrow(result) == 0) stop("No region left after intersection", call. = FALSE)
  result$region.1 <- NULL
  if (is.null(name))
    result$region <- sprintf("(overlap of %s and %s)", x$region, y$region)
  else
    result$region <- name
  attrs <- if (!is.null(attr(x, "map"))) "map" else NULL
  result <- copy_attributes(result, x, attrs)
  sf::st_agr(result) <- "constant"
  result
}


#' Generate the difference between two \code{slendr_region} objects
#'
#' @inheritParams join
#'
#' @return Object of the class \code{slendr_region}
#'
#' @export
subtract <- function(x, y, name = NULL) {
  if (!inherits(x, "slendr_region")) x <- region(polygon = x)
  if (!inherits(y, "slendr_region")) y <- region(polygon = y)

  result <- sf::st_difference(x, y)
  if (nrow(result) == 0) stop("No region left after subtraction", call. = FALSE)
  result$region.1 <- NULL
  if (is.null(name))
    result$region <- sprintf("(%s minus %s)", x$region, y$region)
  else
    result$region <- name
  attrs <- if (!is.null(attr(x, "map"))) "map" else NULL
  result <- copy_attributes(result, x, attrs)
  sf::st_agr(result) <- "constant"
  result
}


#' Calculate the distance between a pair of spatial boundaries
#'
#' @param x,y Objects of the class \code{slendr}
#' @param measure How to measure distance? This can be either \code{'border'}
#'   (distance between the borders of \code{x} and \code{y}) or \code{'center'}
#'   (distance between their centroids).
#' @param time Time closest to the spatial maps of \code{x} and \code{y} if they
#'   represent \code{slendr_pop} population boundaries (ignored for general
#'   \code{slendr_region} objects)
#'
#' @return If the coordinate reference system was specified, a distance in
#'   projected units (i.e. meters) is returned. Otherwise the function returns a
#'   normal Euclidian distance.
#'
#' @export
distance <- function(x, y, measure, time = NULL) {
  if (!measure %in% c("center", "border"))
    stop("Unknown distance measure method provided (must be either 'center' or 'border')", call. = FALSE)

  if ((inherits(x, "slendr_pop") | inherits(y, "slendr_pop")) & is.null(time))
    stop("Time of the nearest spatial snapshot must be be given when calculating
distance between population boundaries", call. = FALSE)

  if (inherits(x, "slendr_pop")) x <- x[which.min(abs(x$time - time)), ]
  if (inherits(y, "slendr_pop")) y <- y[which.min(abs(y$time - time)), ]

  if (measure == "center") {
    x <- sf::st_centroid(x)
    y <- sf::st_centroid(y)
  }

  as.vector(sf::st_distance(x, y))
}


#' Return the dimensions of the world map
#'
#' @param map Object of the type \code{slendr_map}
#'
#' @return If the coordinate reference system was specified, the function
#'   returns a two-dimensional vector of the world dimensions in the projected
#'   units. Otherwise the dimensions of the two-dimensional abstract plane are
#'   returned.
#'
#' @export
dimension <- function(map) {
  if (!inherits(map, "slendr_map"))
    stop("Incorrect input type. Object of the type 'slendr_map' expected", call. = FALSE)
  c(as.vector(diff(sf::st_bbox(map)[c("xmin", "xmax")])),
    as.vector(diff(sf::st_bbox(map)[c("ymin", "ymax")])))
}


# Internal implementation of expand() and shrink() functions
shrink_or_expand <- function(pop, by, end, start, overlap, snapshots, polygon, verbose) {
  check_event_time(c(start, end), pop)
  check_removal_time(start, pop)

  map <- attr(pop, "map")

  # get the last available population boundary
  region_start <- pop[nrow(pop), ]
  sf::st_agr(region_start) <- "constant"

  if (is.null(snapshots)) {
    n <- 1
    message("Iterative search for the minimum sufficient number of intermediate
spatial snapshots, starting at ", n, ". This should only take a couple of
seconds, but if you don't want to wait, you can set `snapshots = N` manually.")
  } else
    n <- snapshots

  # iterate through the number of intermediate spatial boundaries to reach
  # the required overlap between subsequent spatial maps
  repeat {
    times <- seq(start, end, length.out = n + 1)[-1]

    # generate intermediate spatial maps, starting from the last one
    inter_regions <- list()
    inter_regions[[1]] <- region_start
    for (i in seq_along(times)) {
      exp_region <- sf::st_buffer(inter_regions[[1]], dist = i * (by / n))
      exp_region$time <- times[i]
      sf::st_agr(exp_region) <- "constant"

      # restrict the next spatial boundary to the region of interest
      if (!is.null(polygon)) {
        if (!inherits(polygon, "slendr_region"))
          polygon <- region(polygon = polygon, map = map)
        exp_region <- sf::st_intersection(exp_region, polygon)
        exp_region$region <- NULL
        sf::st_agr(exp_region) <- "constant"
      }

      inter_regions[[i + 1]] <- exp_region
    }

    if (!is.null(snapshots)) break

    overlaps <- compute_overlaps(do.call(rbind, inter_regions))

    if (all(overlaps >= overlap)) {
      message("The required ", sprintf("%.1f%%", 100 * overlap),
              " overlap between subsequent spatial maps has been met")
      break
    } else {
      n <- n + 1
      if (verbose)
        message("- Increasing to ", n, " snapshots")
    }
  }

  all_maps <- do.call(rbind, inter_regions) %>% rbind(pop, .)
  sf::st_agr(all_maps) <- "constant"

  result <- copy_attributes(
    all_maps, pop,
    c("map", "parent", "remove", "intersect", "aquatic", "history")
  )

  attr(result, "history") <- append(attr(result, "history"), list(data.frame(
    pop =  unique(region_start$pop),
    event = "expand",
    start = start,
    end = end
  )))

  result
}
