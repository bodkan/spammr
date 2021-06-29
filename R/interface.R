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
#' @param offspring_dist Standard deviation of the normal distribution of the
#'   parent-offspring distance
#' @param aquatic Is the species aquatic (\code{FALSE} by default, i.e.
#'   terrestrial species)?
#'
#' @return Object of the class \code{slendr_pop}
#'
#' @export
population <- function(name, time, N, parent = "ancestor", map = NULL,
                       center = NULL, radius = NULL, polygon = NULL,
                       remove = NULL, intersect = TRUE,
                       competition_dist = NULL, mate_dist = NULL, offspring_dist = NULL,
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
    data.frame(pop = name, time = time, N = N, stringsAsFactors = FALSE),
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

  attr(pop, "competition_dist") <- competition_dist
  attr(pop, "mate_dist") <- mate_dist
  attr(pop, "offspring_dist") <- offspring_dist

  attr(pop, "aquatic") <- aquatic

  class(pop) <- set_class(pop, "pop")

  pop
}


#' Update the population map or one of its parameters
#'
#' This function allows a more manual control of spatial map changes
#' in addition to the \code{expand} and \code{move} functions
#'
#' @param pop Object of the class \code{slendr_pop}
#' @param time Time of the change
#' @param N Number of individuals
#' @param center Two-dimensional vector specifying the center of the
#'   circular range
#' @param radius Radius of the circular range
#' @param polygon List of vector pairs, defining corners of the
#'   polygon range (see also the \code{region} argument) or a
#'   geographic region of the class \code{slendr_region} from which
#'   the polygon coordinates will be extracted
#'
#' @return Object of the class \code{slendr_pop}
#'
#' @export
change <- function(pop, time, N = NULL,
                   center = NULL, radius = NULL, polygon = NULL) {
  check_event_time(start, pop)
  check_removal_time(start, pop)

  map <- attr(pop, "map")

  # define the population range as a simple geometry object
  #or reuse the old one
  if (!is.null(polygon) & inherits(polygon, "slendr_region"))
      polygon <- sf::st_geometry(polygon)

  if ((!is.null(center) & !is.null(radius)) | !is.null(polygon))
    boundary <- define_boundary(map, center, radius, NULL)
  else
    boundary <- sf::st_geometry(pop[nrow(pop), ])

  if (is.null(N)) N <- pop[nrow(pop), ]$N

  updated <- sf::st_sf(
    data.frame(pop = unique(pop$pop), time = time, N = N, stringsAsFactors = FALSE),
    geometry = boundary
  )

  combined <- rbind(pop, updated)
  sf::st_agr(combined) <- "constant"

  result <- copy_attributes(
    combined, pop,
    c("map", "parent", "remove", "intersect", "competition_dist",
      "mate_dist", "offspring_dist", "aquatic")
  )

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
seconds but if you don't want to wait, you can set `snapshots = N` manually.")
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

    # calculate the overlap between subsequent spatial snapshots
    overlaps <- sapply(
      seq_along(inter_regions)[-1], function(i) {
        a <- inter_regions[[i - 1]]
        b <- inter_regions[[i]]
        sf::st_area(sf::st_intersection(a, b)) / sf::st_area(b)
      }
    )

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
    c("map", "parent", "remove", "intersect", "competition_dist",
      "mate_dist", "offspring_dist", "aquatic")
  )

  result
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
    stop("The required overlap between subsequent spatial maps must be a number
between 0 and 1", call. = FALSE)

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
    message("Iterative search for the minimum sufficient number of intermediate
spatial snapshots, starting at ", n, ". This should only take a couple of
seconds but if you don't want to wait, you can set `snapshots = N` manually.")
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
          N = region_start$N,
          stringsAsFactors = FALSE
        ),
        geometry = shifted_region,
        crs = sf::st_crs(inter_regions[[i]])
      )
      sf::st_agr(inter_regions[[i + 1]]) <- "constant"
    }

    if (!is.null(snapshots)) break

    # calculate the overlap between subsequent spatial snapshots
    overlaps <- sapply(
      seq_along(inter_regions)[-1], function(i) {
        a <- inter_regions[[i - 1]]
        b <- inter_regions[[i]]
        intersection <- sf::st_intersection(a, b)
        if (nrow(intersection) == 0) return(0)
        sf::st_area(intersection) / sf::st_area(b)
      }
    )

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
    c("map", "parent", "remove", "intersect", "competition_dist",
      "mate_dist", "offspring_dist", "aquatic")
  )

  result
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
    stop(sprintf("Specified times are not consistent with the assumed direction of
time (geneflow %s -> %s in the time window %s-%s)", from_name, to_name, start, end),
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


#' Print a summary of a \code{slendr} object
#'
#' Prints a short summary of any \code{slendr} object.
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
#' snapshots can be printed using the \code{full} argument.
#'
#' @param x Object of a class \code{slendr}
#'
#' @export
print.slendr <- function(x, ...) {
  if (any(grepl("slendr_pop", class(x))))
    type <- "population"
  else if (any(grepl("slendr_region", class(x))))
    type <- "region"
  else if (any(grepl("slendr_map", class(x))))
    type <- "map"
  else if (any(grepl("slendr_model", class(x))))
    type <- "model"
  else
    stop("Unknown object")

  header <- sprintf("slendr '%s' object", type)
  sep <- paste(rep("-", nchar(header)), collapse = "")

  cat(header, "\n")
  cat(sep, "\n")

  if (type == "region")
    cat("name:", x$region, "\n\n")

  if (type == "population") {
    cat("name:", unique(x$pop), "\n")
    parent <- attr(x, "parent")
    cat("split from: ")
    if (is.character(parent) && parent == "ancestor")
      cat("[this is an ancestral population]\n")
    else {
      cat(parent$pop[1], "\n")
      cat("split time:", x$time[1], "\n")
    }
    cat("removed at: ")
    if (attr(x, "remove") == -1)
      cat("[will not be removed]\n")
    else
      cat((attr(x, "remove")), "\n")

    if (attr(x, "aquatic") == FALSE)
      cat("habitat: terrestrial\n")
    else
      cat("habitat: aquatic\n")

    # pretty print the raw sf data as a simplified table
    cat("snapshots:\n")
    # TODO: as.data.frame(x) is giving a very strange error here (on
    # the intersected ANA expansion sf object in plotting) - a bug
    # in sf conversion? manual creation of the data.frame does not
    # have this issue
    snapshots_df <- unique(data.frame(
      pop = x$pop,
      time = x$time,
      N = x$N,
      stringsAsFactors = FALSE
    ))
    snapshots_df$`#` <- seq_len(nrow(snapshots_df))
    # determine which maps over time are new and which are re-used from the
    # previous time point (we do this because the raw spatial geometry
    # representation is hard to read and not useful for seeing what changes
    # when)
    snapshots_df <- snapshots_df[, c("#", "time", "N")]
    print(snapshots_df, row.names = FALSE)
    cat("\n")
  }

  if (type %in% c("map", "region", "population")) {
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
  } else if (type == "model") {
    cat("populations:", paste0(x$splits$pop, collapse = ", "), "\n")
    cat("geneflow events: ")
    if (!is.null(x$geneflows))
      cat(nrow(x$geneflows), "\n")
    else
      cat("[no geneflow]\n")
    cat("generation time:", x$generation_time, "\n")
    cat("time direction:", get_time_direction(tail(x$populations, 1)[[1]]), "\n")
    cat("number of spatial maps:", nrow(x$maps), "\n")
    cat("resolution:", x$resolution, "distance unit per pixel\n\n")
    cat("configuration files in:", normalizePath(x$config$directory), "\n\n")
    cat(
"A detailed model specification can be found in `$splits`, `$geneflows`,
`$maps`, `$populations`, and other components of the model object (for
a complete list see `names(<model object>)`). You can also examine
the serialized configuration files in the model directory.\n")
  } else {
    stop("Unknown object type", call. = FALSE)
  }
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
