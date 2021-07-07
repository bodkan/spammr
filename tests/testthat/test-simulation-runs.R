# Run a given forward/backward simulation with different combinations
# of start-end and burnin parameters
run_sim <- function(direction, start, burnin, sim_length = NULL) {
  map <- world(xrange = c(0, 10), yrange = c(0, 10), landscape = "blank")
  p <- population("pop", N = 5, time = start, map = map, center = c(5, 5), radius = 1)
  model <- compile(p, direction = direction, sim_length = sim_length, dir = tempdir(), generation_time = 1, resolution = 1, competition_dist = 10, mate_dist = 10, dispersal_dist = 10, overwrite = TRUE)
  slim(model, method = "batch", burnin = burnin, seq_length = 1, recomb_rate = 0, save_locations = TRUE)
  model
}

# forward simulations

test_that("Forward simulation from generation 1 has the correct length without burnin", {
  direction <- "forward"; start <- 1; sim_length <- 10; burnin <- 0
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(min(inds$t) == start)
  expect_true(max(inds$t) == (start + sim_length - 1))
  expect_true(length(unique(inds$t)) == sim_length)
})

test_that("Forward simulation from generation 1 has the correct length with burnin", {
  direction <- "forward"; start <- 1; sim_length <- 10; burnin <- 100
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(min(inds$t) == start)
  expect_true(max(inds$t) == (start + sim_length - 1))
  expect_true(length(unique(inds$t)) == sim_length)
})

test_that("Forward simulation from generation > 1 has the correct length without burnin", {
  direction <- "forward"; start <- 8; sim_length <- 10; burnin <- 0
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(min(inds$t) == start)
  expect_true(max(inds$t) == (start + sim_length - 1))
  expect_true(length(unique(inds$t)) == sim_length)
})

test_that("Forward simulation from generation > 1 has the correct length with burnin", {
  direction <- "forward"; start <- 8; sim_length <- 10; burnin <- 100
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(min(inds$t) == start)
  expect_true(max(inds$t) == (start + sim_length - 1))
  expect_true(length(unique(inds$t)) == sim_length)
})

# backward simulations

test_that("Backward simulation has the correct length without burnin", {
  direction <- "backward"; start <- 10; burnin <- 0
  result <- run_sim(direction, start, burnin)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(max(inds$t) == start)
  expect_true(min(inds$t) == 1)
  expect_true(length(unique(inds$t)) == start)
})

test_that("Backward simulation has the correct length with burnin", {
  direction <- "backward"; start <- 10; burnin <- 100
  result <- run_sim(direction, start, burnin)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(max(inds$t) == start)
  expect_true(min(inds$t) == 1)
  expect_true(length(unique(inds$t)) == start)
})

test_that("Backward simulation of limited length has the correct length without burnin", {
  direction <- "backward"; start <- 10; sim_length = 7; burnin <- 0
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(max(inds$t) == start)
  expect_true(min(inds$t) == start - sim_length + 1)
  expect_true(length(unique(inds$t)) == sim_length)
})

test_that("Backward simulation of limited length has the correct length with burnin", {
  direction <- "backward"; start <- 10; sim_length = 7; burnin <- 100
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(max(inds$t) == start)
  expect_true(min(inds$t) == start - sim_length + 1)
  expect_true(length(unique(inds$t)) == sim_length)
})
