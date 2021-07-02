test_that("forward and backward time model objects are equivalent", {
  map <- world(xrange = c(0, 100), yrange = c(0, 100), landscape = "blank")

  # forward direction
  p1 <- population(name = "p1",              time = 1, N = 1, center = c(1, 1), radius = 1, map = map)
  p2 <- population(name = "p2", parent = p1, time = 2, N = 1, center = c(1, 1), radius = 1)
  p3 <- population(name = "p3", parent = p1, time = 3, N = 1, center = c(1, 1), radius = 1)
  p4 <- population(name = "p4", parent = p3, time = 4, N = 1, center = c(1, 1), radius = 1)
  p5 <- population(name = "p5", parent = p2, time = 5, N = 1, center = c(1, 1), radius = 1)

  geneflows <- geneflow(p1, p2, start = 2, end = 3, rate = 0.5, overlap = FALSE)

  forward <- compile(
    dir = file.path(tempdir(), "tmp-forward"),
    populations = list(p1, p2, p3, p4, p5), geneflow = geneflows,
    generation_time = 1, resolution = 1,
    overwrite = TRUE,
    competition_dist = 1, mate_dist = 1, offspring_dist = 1,
    sim_length = 5
  )

  # backward direction
  p1 <- population(name = "p1",              time = 5, N = 1, center = c(1, 1), radius = 1, map = map)
  p2 <- population(name = "p2", parent = p1, time = 4, N = 1, center = c(1, 1), radius = 1)
  p3 <- population(name = "p3", parent = p1, time = 3, N = 1, center = c(1, 1), radius = 1)
  p4 <- population(name = "p4", parent = p3, time = 2, N = 1, center = c(1, 1), radius = 1)
  p5 <- population(name = "p5", parent = p2, time = 1, N = 1, center = c(1, 1), radius = 1)

  geneflows <- geneflow(p1, p2, start = 4, end = 3, rate = 0.5, overlap = FALSE)

  backward <- compile(
    dir = file.path(tempdir(), "tmp-backward"),
    populations = list(p1, p2, p3, p4, p5), geneflow = geneflows,
    generation_time = 1, resolution = 1,
    overwrite = TRUE,
    competition_dist = 1, mate_dist = 1, offspring_dist = 1
  )

  expect_true(all.equal(forward$splits[, grep("_orig", colnames(forward$splits), value = TRUE, invert = TRUE)],
                        backward$splits[, grep("_orig", colnames(backward$splits), value = TRUE, invert = TRUE)]))
  expect_true(all.equal(forward$geneflow[, grep("_orig", colnames(forward$geneflow), value = TRUE, invert = TRUE)],
                        backward$geneflow[, grep("_orig", colnames(forward$geneflow), value = TRUE, invert = TRUE)]))
  expect_true(all.equal(forward$maps[, c("pop", "pop_id", "time_gen")],
                        backward$maps[, c("pop", "pop_id","time_gen")]))

  components <- c("generation_time", "resolution", "world")
  expect_true(all(sapply(components, function(i) all.equal(forward[[i]], backward[[i]]))))
})


test_that("forward and backward models yield the same simulation result", {
  map <- readRDS("map.rds")

  # forward simulation ------------------------------------------------------
  p1 <- population(name = "pop1", parent = "ancestor", N = 100, time = 1, radius = 600000, center = c(10, 25), map = map)
  p2 <- population(name = "pop2", parent = p1, time = 10, N = 100, center = c(10, 25), radius = 300000) %>%
    move(trajectory = list(c(25, 25), c(40, 30), c(40, 40), c(50, 50)), start = 10, end = 200, snapshots = 30)
  p3 <- population(name = "pop3", parent = p2, time = 200, N = 100, polygon = list(c(-10, 50), c(10, 50), c(20, 53), c(40, 60), c(40, 70), c(-10, 65)))
  p4 <- population(name = "pop4", parent = p2, time = 200, N = 100, polygon = list(c(-10, 35), c(20, 37), c(25, 40), c(30, 45), c(10, 50), c(-10, 45)))
  p5 <- population(name = "pop5", parent = p1, time = 300, N = 100, center = c(10, 25), radius = 300000) %>%
    move(trajectory = list(c(-5, 33), c(-5, 40)), start = 300, end = 400, snapshots = 20)

  geneflow <- list(
    geneflow(from = p5, to = p4, rate = 0.2, start = 310, end = 350, overlap = F),
    geneflow(from = p5, to = p3, rate = 0.3, start = 340, end = 480, overlap = F)
  )

  forward <- compile(
    dir = file.path(tempdir(), "tmp-forward"),
    populations = list(p1, p2, p3, p4, p5),
    geneflow = geneflow,
    generation_time = 1,
    resolution = 10000,
    overwrite = TRUE,
    competition_dist = 100e3, mate_dist = 100e3, offspring_dist = 100e3,
    sim_length = 480
  )

  # backward simulation -----------------------------------------------------
  p1 <- population(name = "pop1", parent = "ancestor", N = 100, time = 480, radius = 600000, center = c(10, 25), map = map)
  p2 <- population(name = "pop2", parent = p1, time = 471, N = 100, center = c(10, 25), radius = 300000) %>%
    move(trajectory = list(c(25, 25), c(40, 30), c(40, 40), c(50, 50)), start = 471, end = 281, snapshots = 30)
  p3 <- population(name = "pop3", parent = p2, time = 281, N = 100, polygon = list(c(-10, 50), c(10, 50), c(20, 53), c(40, 60), c(40, 70), c(-10, 65)))
  p4 <- population(name = "pop4", parent = p2, time = 281, N = 100, polygon = list(c(-10, 35), c(20, 37), c(25, 40), c(30, 45), c(10, 50), c(-10, 45)))
  p5 <- population(name = "pop5", parent = p1, time = 181, N = 100, center = c(10, 25), radius = 300000) %>%
    move(trajectory = list(c(-5, 33), c(-5, 40)), start = 181, end = 81, snapshots = 20)

  geneflow <- list(
    geneflow(from = p5, to = p4, rate = 0.2, start = 171, end = 131, overlap = F),
    geneflow(from = p5, to = p3, rate = 0.3, start = 141, end = 1, overlap = F)
  )

  backward <- compile(
    dir = file.path(tempdir(), "tmp-backward"),
    populations = list(p1, p2, p3, p4, p5),
    geneflow = geneflow,
    generation_time = 1,
    resolution = 10000,
    overwrite = TRUE,
    competition_dist = 100e3, mate_dist = 100e3, offspring_dist = 100e3
  )

  expect_true(all.equal(forward$splits[, grep("_orig", colnames(forward$splits), value = TRUE, invert = TRUE)],
                        backward$splits[, grep("_orig", colnames(backward$splits), value = TRUE, invert = TRUE)]))
  expect_true(all.equal(forward$geneflow[, grep("_orig", colnames(forward$geneflow), value = TRUE, invert = TRUE)],
                        backward$geneflow[, grep("_orig", colnames(forward$geneflow), value = TRUE, invert = TRUE)]))
  expect_true(all.equal(forward$maps[, c("pop", "pop_id", "time_gen")],
                        backward$maps[, c("pop", "pop_id","time_gen")]))

  components <- c("generation_time", "resolution", "world")
  expect_true(all(sapply(components, function(i) all.equal(forward[[i]], backward[[i]]))))

  # simulation runs are the same
  slim(forward, seq_length = 1, recomb_rate = 0, save_locations = TRUE, method = "batch", seed = 123)
  slim(backward, seq_length = 1, recomb_rate = 0, save_locations = TRUE, method = "batch", seed = 123)

  # make sure the scripts are the same
  f_script <- file.path(forward$directory, "script.slim") %>%
    readLines %>% grep("MODEL_DIR|OUTPUT_PREFIX", ., value = TRUE, invert = TRUE)
  b_script <- file.path(backward$directory, "script.slim") %>%
    readLines %>% grep("MODEL_DIR|OUTPUT_PREFIX", ., value = TRUE, invert = TRUE)
  expect_equal(f_script, b_script)

  # make sure that the simulated location data is the same
  f_loc <- data.table::fread(file.path(forward$directory, "output_ind_locations.tsv.gz"))
  b_loc <- data.table::fread(file.path(backward$directory, "output_ind_locations.tsv.gz"))

  expect_equal(f_loc, b_loc)
})

test_that("move preceding population split results in an error", {
  map <- world(xrange = c(0, 100), yrange = c(0, 100), landscape = "blank")
  p1 <- population(name = "p1", time = 1, N = 1, center = c(1, 1), radius = 10, map = map)
  p2 <- population(name = "p2", parent = p1, time = 10, N = 1, center = c(20, 50), radius = 20)
  expect_error(move(p2, trajectory = c(100, 100), start = 3, end = 20, snapshots = 5),
               "The specified event time")
})

test_that("overlaps in time result in an error (forward time)", {
  map <- world(xrange = c(0, 100), yrange = c(0, 100), landscape = "blank")
  p1 <- population(name = "p1", time = 1, N = 1, center = c(1, 1), radius = 10, map = map)
  p2 <- population(name = "p2", parent = p1, time = 10, N = 1, center = c(20, 50), radius = 20) %>%
    move(trajectory = c(100, 100), start = 10, end = 20, snapshots = 5)
  # repeating the same move
  expect_error(move(p2, trajectory = c(100, 100), start = 10, end = 20, snapshots = 5),
               "The specified event time")
  # overlapping time window of the second move
  expect_error(move(p2, trajectory = c(100, 100), start = 5, end = 30, snapshots = 5),
               "The specified event time")
  # overlapping time window of the following range expansion
  expect_error(expand(p2, by = 20, start = 5, end = 30, snapshots = 5),
               "The specified event time" )
  # overlapping time window of the following range change
  expect_error(boundary(p2, time = 8, center = c(50, 10), radius = 8),
               "The specified event time")
})

test_that("overlaps in time result in an error (backward time)", {
  map <- world(xrange = c(0, 100), yrange = c(0, 100), landscape = "blank")
  p1 <- population(name = "p1", time = 100, N = 1, center = c(1, 1), radius = 10, map = map)
  p2 <- population(name = "p2", parent = p1, time = 40, N = 1, center = c(20, 50), radius = 20) %>%
    move(trajectory = c(100, 100), start = 30, end = 10, snapshots = 5)
  # repeating the same move
  expect_error(move(p2, trajectory = c(100, 100), start = 30, end = 10, snapshots = 5),
               "The specified event time")
  # overlapping time window of the second move
  expect_error(move(p2, trajectory = c(100, 100), start = 25, end = 5, snapshots = 5),
              "The specified event time")
  # overlapping time window of the following range expansion
  expect_error(expand(p2, by = 20, start = 25, end = 5, snapshots = 5),
               "The specified event time")
  # overlapping time window of the following range change
  expect_error(boundary(p2, time = 15, center = c(50, 10), radius = 8),
               "The specified event time")
})

test_that("times are correctly specified for different population size change methods", {
  map <- world(xrange = c(0, 100), yrange = c(0, 100), landscape = "blank")
  p <- population(name = "pop", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  expect_error(resize(p, N = 10, how = "step"), "Timing of the population change needs to be specified")
  expect_silent(resize(p, N = 10, how = "step", time = 123))
  expect_error(resize(p, N = 10, how = "linear"), "Timing of the population change needs to be specified")
  expect_silent(resize(p, N = 10, how = "linear", start = 123, end = 10))
  expect_error(resize(p, N = 10, how = "exponential"), "Timing of the population change needs to be specified")
  expect_silent(resize(p, N = 10, how = "exponential", start = 123, end = 10))
})

test_that("resizing of populations is consistent with established population dynamics (forward time)", {
  map <- world(xrange = c(0, 100), yrange = c(0, 100), landscape = "blank")
  p1 <- population(name = "pop1", map = map, time = 1, N = 500, center = c(10, 25), radius = 300000)
  p2 <- population(name = "pop2", parent = p1, time = 10, N = 500, center = c(10, 25), radius = 300000)
  expect_error(resize(p2, N = 10, how = "step", time = 5), "The specified event time")
})

test_that("resizing of populations is consistent with established population dynamics (backward time)", {
  map <- world(xrange = c(0, 100), yrange = c(0, 100), landscape = "blank")
  p1 <- population(name = "pop1", map = map, time = 30000, N = 500, center = c(10, 25), radius = 300000)
  p2 <- population(name = "pop2", parent = p1, time = 25000, N = 500, center = c(10, 25), radius = 300000)
  expect_error(resize(p2, N = 10, how = "step", time = 28000), "The specified event time")
})
