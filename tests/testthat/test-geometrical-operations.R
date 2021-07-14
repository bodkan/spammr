test_that("circular ranges need two parameters", {
  map <- readRDS("map.rds")
  expect_error(
    population("pop", parent = "ancestor", time = 1, N = 1, map = map, center = 1),
    "Missing radius argument when defining a circular population range"
  )
  expect_error(
    population("pop", parent = "ancestor", time = 1, N = 1, map = map, radius = 1),
    "Missing center argument when defining a circular population range"
  )
})

test_that("polygon boundary can be specified in two ways", {
  map <- readRDS("map.rds")
  coords <- list(c(-18, 20), c(40, 20), c(30, 33), c(20, 32), c(10, 35), c(-8, 35))
  africa <- region("Africa", map, polygon = coords)
  afr1 <- population("AFR", parent = "ancestor", time = 1, N = 1, map = map, polygon = africa)
  afr2 <- population("AFR", parent = "ancestor", time = 1, N = 1, map = map, polygon = coords)
  expect_equal(afr1, afr2)
})

test_that("polygon ranges need at least three corners", {
  map <- readRDS("map.rds")
  coords <- list(c(-18, 20), c(40, 20))
  expect_error(
    region("Africa", map, polygon = coords),
    "Polygon range needs to have at least three corners"
  )
})

test_that("polygon ranges need at least three corners", {
  map <- readRDS("map.rds")
  coords <- list(c(-18, 20), c(40, 20), c(50, 50))
  expect_error(
    population("AFR", parent = "ancestor", time = 1, N = 1, map = map, center = 1, radius = 1, polygon = coords),
    "Either a circular range (center and radius) or the corners of
a polygon need to be specified, but not both", fixed = TRUE
  )
})

# set operations

test_that("union of a region with itself gives the same region", {
  map <- world(xrange = c(-15, 60), yrange = c(-20, 65), landscape = "blank")
  coords <- list(c(-18, 20), c(40, 20), c(30, 33), c(20, 32), c(10, 35), c(-8, 35))
  africa <- region("Africa", map, polygon = coords)
  expect_true(all(africa == join(africa, africa, name = "Africa")))
})

test_that("intersection of a region with itself gives the same region", {
  map <- world(xrange = c(-15, 60), yrange = c(-20, 65), landscape = "blank")
  coords <- list(c(-18, 20), c(40, 20), c(30, 33), c(20, 32), c(10, 35), c(-8, 35))
  africa <- region("Africa", map, polygon = coords)
  expect_true(all(africa == overlap(africa, africa, name = "Africa")))
})

test_that("subtraction of a region from itself gives an empty region", {
  map <- world(xrange = c(-15, 60), yrange = c(-20, 65), landscape = "blank")
  coords <- list(c(-18, 20), c(40, 20), c(30, 33), c(20, 32), c(10, 35), c(-8, 35))
  africa <- region("Africa", map, polygon = coords)
  expect_true(nrow(subtract(africa, africa)) == 0)
})

test_that("subtraction of a non-overlapping region does not change the result", {
  map <- world(xrange = c(-15, 60), yrange = c(-20, 65), landscape = "blank")
  r1 <- region("r1", map, center = c(0, 0), radius = 5)
  r2 <- region("r2", map, center = c(10, 10), radius = 5)
  expect_true(all(subtract(r1, r2, name = "r1") == r1))
})

test_that("nonoverlapping join gives empty intersection", {
  map <- world(xrange = c(-15, 60), yrange = c(-20, 65), landscape = "blank")
  r1 <- region("r1", map, center = c(0, 0), radius = 5)
  r2 <- region("r2", map, center = c(10, 10), radius = 5)
  expect_true(nrow(overlap(r1, r2, "")) == 0)
})

test_that("custom landscapes can be specified", {
  r1 <- region(center = c(0, 0), radius = 10)
  r2 <- region(center = c(10, 10), radius = 10)
  r3 <- region(center = c(-10, -10), radius = 10)
  comb <- r1 %>% join(r2) %>% join(r3)
  map <- world(xrange = c(-20, 20), yrange = c(-20, 20), landscape = comb)
  # is the geometry component of the builtin map equal to the are of
  # the original regions?
  expect_true(sf::st_area(map$landscape) == sf::st_area(comb))
})

test_that("boundaries are either circles or polygons", {
  expect_error(region(), "Either a circular range or a polygon range must be specified")
})



# distances between nonoverlapping ranges ---------------------------------

map <- world(xrange = c(0, 100), yrange = c(0, 100), landscape = "blank")
x1 <- 0; x2 <- 100; radius <- 10
p1 <- population("p1", map = map, time = 1, center = c(x1, 50), radius = radius, N = 1)
p2 <- population("p2", map = map, time = 1, center = c(x2, 50), radius = radius, N = 1)

test_that("distance calculation method must be provided", {
  expect_error(distance(p1, p2, measure = "blah"),
               "Unknown distance measure method")
  expect_silent(distance(p1, p2, measure = "center", 0))
})

test_that("time must be given when calculating population distances", {
  expect_error(distance(p1, p2, measure = "border"),
               "Time of the nearest spatial snapshot must be be given")
})

test_that("distance between borders of geographic regions is correctly calculated", {
  expect_true(distance(p1, p2, measure = "border", time = 0) == (x2 - x1) - 2*radius)
})

test_that("distance between centers of geographic regions is correctly calculated", {
  expect_true(all.equal(distance(p1, p2, measure = "center", time = 0), 100))
})


# distances between overlapping ranges ------------------------------------

radius <- 60
p1 <- population("p1", map = map, time = 1, center = c(0, 50), radius = radius, N = 1)
p2 <- population("p2", map = map, time = 1, center = c(100, 50), radius = radius, N = 1)

test_that("distance between population borders is zero", {
  expect_true(distance(p1, p2, measure = "border", time = 0) == 0)
})

test_that("distance between population centers is still the same", {
  expect_true(all.equal(distance(p1, p2, measure = "center", time = 0), 100))
})

# distances can be calculated even for generic regions --------------------
r1 <- region(map = map, center = c(0, 50), radius = 60)
r2 <- region(map = map, center = c(100, 50), radius = 60)

test_that("distance between borders of regions is zero", {
  expect_true(distance(r1, r2, measure = "border") == 0)
})

test_that("distance between centers of regions is the same", {
  expect_true(distance(r1, r2, measure = "center") == 100)
})

test_that("it is possible to calculate distances betwewen populations and regions", {
  expect_true(distance(p1, r2, measure = "border", time = 0) == 0)
  expect_true(distance(r1, p2, measure = "border", time = 0) == 0)
  expect_true(distance(p1, r2, measure = "center", time = 0) == 100)
  expect_true(distance(r1, p2, measure = "center", time = 0) == 100)
})


test_that("world dimensions are correctly calculated", {
  map1 <- world(xrange = c(0, 100), yrange = c(0, 100), landscape = "blank")
  map2 <- readRDS("map.rds"); bbox <- sf::st_bbox(map2)
  expect_equal(dimension(map1), c(100, 100))
  expect_equal(dimension(map2), as.vector(c(bbox[3] - bbox[1], bbox[4] - bbox[2])))
})
