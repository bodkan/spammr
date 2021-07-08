# Run a given forward/backward simulation with different combinations
# of start-end and burnin parameters
run_sim <- function(direction, start, burnin, sim_length = NULL, gen_time = 1) {
  map <- world(xrange = c(0, 10), yrange = c(0, 10), landscape = "blank")
  p <- population("pop", N = 5, time = start, map = map, center = c(5, 5), radius = 1)
  model <- compile(p, direction = direction, sim_length = sim_length, dir = tempdir(), generation_time = gen_time, resolution = 1, competition_dist = 10, mate_dist = 10, dispersal_dist = 10, overwrite = TRUE)
  slim(model, method = "batch", burnin = burnin, seq_length = 1, recomb_rate = 0, save_locations = TRUE)
  model
}


# forward simulations - generation time = 1 -------------------------------

test_that("Forward simulation from generation 1 has the correct length without burnin", {
  direction <- "forward"; start <- 1; sim_length <- 5; burnin <- 0
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(min(inds$time) == start)
  expect_true(max(inds$time) == (start + sim_length - 1))
  expect_true(length(unique(inds$time)) == sim_length)
})

test_that("Forward simulation from generation 1 has the correct length with burnin", {
  direction <- "forward"; start <- 1; sim_length <- 5; burnin <- 20
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(min(inds$time) == start)
  expect_true(max(inds$time) == (start + sim_length - 1))
  expect_true(length(unique(inds$time)) == sim_length)
})

test_that("Forward simulation from generation > 1 has the correct length without burnin", {
  direction <- "forward"; start <- 8; sim_length <- 5; burnin <- 0
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(min(inds$time) == start)
  expect_true(max(inds$time) == (start + sim_length - 1))
  expect_true(length(unique(inds$time)) == sim_length)
})

test_that("Forward simulation from generation > 1 has the correct length with burnin", {
  direction <- "forward"; start <- 8; sim_length <- 5; burnin <- 100
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(min(inds$time) == start)
  expect_true(max(inds$time) == (start + sim_length - 1))
  expect_true(length(unique(inds$time)) == sim_length)
})

# backward simulations - generation time = 1 ------------------------------

test_that("Backward simulation has the correct length without burnin", {
  direction <- "backward"; start <- 5; burnin <- 0
  result <- run_sim(direction, start, burnin)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(max(inds$time) == start)
  expect_true(min(inds$time) == 0)
  expect_true(length(unique(inds$time)) == start + 1)
})

test_that("Backward simulation has the correct length with burnin", {
  direction <- "backward"; start <- 5; burnin <- 20
  result <- run_sim(direction, start, burnin)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(max(inds$time) == start)
  expect_true(min(inds$time) == 0)
  expect_true(length(unique(inds$time)) == start + 1)
})

test_that("Backward simulation of limited length has the correct length without burnin", {
  direction <- "backward"; start <- 5; sim_length = 3; burnin <- 0
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(max(inds$time) == start)
  expect_true(min(inds$time) == start - sim_length + 1)
  expect_true(length(unique(inds$time)) == sim_length)
})

test_that("Backward simulation of limited length has the correct length with burnin", {
  direction <- "backward"; start <- 5; sim_length = 3; burnin <- 20
  result <- run_sim(direction, start, burnin, sim_length)
  inds <- fread(file.path(result$directory, "output_ind_locations.tsv.gz")); inds
  expect_true(max(inds$time) == start)
  expect_true(min(inds$time) == start - sim_length + 1)
  expect_true(length(unique(inds$time)) == sim_length)
})
