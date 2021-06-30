map <- readRDS("map.rds")

test_that("move requires overlap between 0 and 1", {
  p <- population(name = "pop2", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  expect_error(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, overlap = -1),
               "The required overlap between subsequent spatial maps must be a number between 0 and 1", fixed = TRUE)
  expect_error(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, overlap = 0),
               "The required overlap between subsequent spatial maps must be a number between 0 and 1", fixed = TRUE)
  expect_silent(suppressMessages(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, overlap = 0.3)))
})

test_that("move requires a positive snapshot number", {
  p <- population(name = "pop2", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  expect_error(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, snapshots = -1),
               "The number of snapshots must be a non-negative integer", fixed = TRUE)
  expect_error(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, snapshots = 0),
               "The number of snapshots must be a non-negative integer", fixed = TRUE)
  expect_silent(move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, snapshots = 3, verbose = FALSE))
})

test_that("move generates the correct number of snapshots", {
  p <- population(name = "pop2", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  nsnapshots <- 3
  moved_p <- move(p, trajectory = list(c(24, 25)), start = 29000, end = 25000, snapshots = nsnapshots, verbose = FALSE)
  expect_true(nrow(moved_p) == nsnapshots + 2)
})

test_that("only positive population sizes allowed", {
  p <- population(name = "pop", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  expect_error(resize(p, N = -1, time = 1000),
               "Only positive population sizes allowed")
})

test_that("only three values of population size changes allowed", {
  p <- population(name = "pop", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  expect_error(resize(p, N = 10, how = "asdf"), "Only 'step', 'exponential' and 'linear' are allowed")
})

test_that("no overlap with a manually specified boundary is caught", {
  map <- readRDS("map.rds")
  africa <- region("Africa", map, polygon = list(c(-18, 20), c(40, 20), c(30, 33),
                                                 c(20, 32), c(10, 35), c(-8, 35)))
  europe <- region("Europe", map, polygon = list(c(-8, 35), c(-5, 36), c(10, 38), c(20, 35), c(25, 35),
                                                 c(33, 45), c(20, 58), c(-5, 60), c(-15, 50)))
  yam <- population(name = "YAM", time = 7000, N = 500, map = map,
                    polygon = list(c(26, 50), c(38, 49), c(48, 50), c(48, 56), c(38, 59), c(26, 56))) %>%
    move(trajectory = list(c(15, 50)), start = 5000, end = 3000, snapshots = 10)
  expect_error(boundary(yam, time = 2800, polygon = africa), "Insufficient overlap")
  expect_silent(boundary(yam, time = 2800, polygon = africa, overlap = 0))
})
