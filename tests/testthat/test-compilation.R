test_that("'competition_dist' must be specified in compile() if missing", {
  map <- readRDS("map.rds")
  p <- population(mate_dist = 10, offspring_dist = 10, name = "pop1", parent = "ancestor", N = 700, time = 40000, radius = 600000, center = c(10, 25), map = map)
  expect_error(compile(populations = list(p), generation_time = 30, resolution = 11e3, dir = tempfile(), overwrite = TRUE),
               "Parameter 'competition_dist' missing", fixed = TRUE)
  expect_silent(compile(competition_dist = 50e3, populations = list(p), generation_time = 30, resolution = 10e3, dir = tempfile(), overwrite = TRUE))
})

test_that("'mate_dist' must be specified in compile() if missing", {
  map <- readRDS("map.rds")
  p <- population(competition_dist = 10, offspring_dist = 10, name = "pop1", parent = "ancestor", N = 700, time = 40000, radius = 600000, center = c(10, 25), map = map)
  expect_error(compile(populations = list(p), generation_time = 30, resolution = 10e3, dir = tempfile(), overwrite = TRUE),
               "Parameter 'mate_dist' missing", fixed = TRUE)
  expect_silent(compile(mate_dist = 50e3, populations = list(p), generation_time = 30, resolution = 10e3, dir = tempfile(), overwrite = TRUE))
})

test_that("'offspring_dist' must be specified in compile() if missing", {
  map <- readRDS("map.rds")
  p <- population(competition_dist = 10, mate_dist = 10, name = "pop1", parent = "ancestor", N = 700, time = 40000, radius = 600000, center = c(10, 25), map = map)
  expect_error(compile(populations = list(p), generation_time = 30, resolution = 10e3, dir = tempfile(), overwrite = TRUE),
               "Parameter 'offspring_dist' missing", fixed = TRUE)
  expect_silent(compile(offspring_dist = 50e3, populations = list(p), generation_time = 30, resolution = 10e3, dir = tempfile(), overwrite = TRUE))
})

test_that("'competition_dist', 'mate_dist', and 'offspring_dist' do not have to be specified in compile() if already present", {
  map <- readRDS("map.rds")
  p <- population(competition_dist = 10, mate_dist = 10, offspring_dist = 10, name = "pop1", parent = "ancestor", N = 700, time = 40000, radius = 600000, center = c(10, 25), map = map)
  expect_silent(compile(populations = list(p), generation_time = 30, resolution = 10e3, dir = tempfile(), overwrite = TRUE))
})
