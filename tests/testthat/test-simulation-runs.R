# Run a given forward/backward simulation with different combinations
# of start-end and burnin parameters
run_sim <- function(direction, start, burnin, gen_time, sim_length = NULL) {
  map <- world(xrange = c(0, 10), yrange = c(0, 10), landscape = "blank")
  p <- population("pop", N = 5, time = start, map = map, center = c(5, 5), radius = 1)
  model <- compile(p, direction = direction, sim_length = sim_length, dir = tempdir(), generation_time = gen_time, resolution = 1, competition_dist = 10, mate_dist = 10, dispersal_dist = 10, overwrite = TRUE)
  slim(model, method = "batch", burnin = burnin, seq_length = 1, recomb_rate = 0, save_locations = TRUE)
  fread(file.path(model$directory, "output_ind_locations.tsv.gz"))
}

# forward simulations - generation time = 1 -------------------------------

test_that("Forward simulation from generation 1 has the correct length without burnin", {
  direction <- "forward"; start <- 1; sim_length <- 5; burnin <- 0; gen_time <- 1
  result <- run_sim(direction, start, burnin, gen_time, sim_length)
  expect_true(min(result$time) == start)
  expect_true(max(result$time) == start + sim_length)
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

test_that("Forward simulation from generation 1 has the correct length with burnin", {
  direction <- "forward"; start <- 1; sim_length <- 5; burnin <- 20; gen_time <- 1
  result <- run_sim(direction, start, burnin, gen_time, sim_length)
  expect_true(min(result$time) == start)
  expect_true(max(result$time) == start + sim_length)
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

test_that("Forward simulation from generation > 1 has the correct length without burnin", {
  direction <- "forward"; start <- 8; sim_length <- 5; burnin <- 0; gen_time <- 1
  result <- run_sim(direction, start, burnin, gen_time, sim_length)
  expect_true(min(result$time) == start)
  expect_true(max(result$time) == start + sim_length)
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

test_that("Forward simulation from generation > 1 has the correct length with burnin", {
  direction <- "forward"; start <- 8; sim_length <- 5; burnin <- 100; gen_time <- 1
  result <- run_sim(direction, start, burnin, gen_time, sim_length)
  expect_true(min(result$time) == start)
  expect_true(max(result$time) == start + sim_length)
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

# backward simulations - generation time = 1 ------------------------------

test_that("Backward simulation has the correct length without burnin", {
  direction <- "backward"; start <- 5; burnin <- 0; gen_time <- 1
  result <- run_sim(direction, start, burnin, gen_time)
  expect_true(max(result$time) == start)
  expect_true(min(result$time) == 0)
  expect_true(length(unique(result$time)) == round(start / gen_time) + 1)
})

test_that("Backward simulation has the correct length with burnin", {
  direction <- "backward"; start <- 5; burnin <- 20; gen_time <- 1
  result <- run_sim(direction, start, burnin, gen_time)
  expect_true(max(result$time) == start)
  expect_true(min(result$time) == 0)
  expect_true(length(unique(result$time)) == round(start / gen_time) + 1)
})

test_that("Backward simulation of limited length has the correct length without burnin", {
  direction <- "backward"; start <- 5; sim_length = 3; burnin <- 0; gen_time <- 1
  result <- run_sim(direction, start, burnin, gen_time, sim_length)
  expect_true(max(result$time) == start)
  expect_true(min(result$time) == start - sim_length)
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

test_that("Backward simulation of limited length has the correct length with burnin", {
  direction <- "backward"; start <- 5; sim_length = 3; burnin <- 20; gen_time <- 1
  result <- run_sim(direction, start, burnin, gen_time, sim_length)
  expect_true(max(result$time) == start)
  expect_true(min(result$time) == start - sim_length)
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

# forward simulations - generation time > 1 -------------------------------

test_that("Forward simulation from generation 1 has the correct length without burnin", {
  direction <- "forward"; start <- 20; sim_length <- 50; burnin <- 0; gen_time <- 20
  result <- run_sim(direction, start, burnin, gen_time, sim_length); result
  expect_true(min(result$time) == start)
  expect_true(max(result$time) == start + gen_time * round(sim_length / gen_time))
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

test_that("Forward simulation from generation 1 has the correct length with burnin", {
  direction <- "forward"; start <- 20; sim_length <- 50; burnin <- 100; gen_time <- 20
  result <- run_sim(direction, start, burnin, gen_time, sim_length); result
  expect_true(min(result$time) == start)
  expect_true(max(result$time) == start + gen_time * round(sim_length / gen_time))
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

test_that("Forward simulation from generation > 1 has the correct length without burnin", {
  direction <- "forward"; start <- 200; sim_length <- 50; burnin <- 0; gen_time <- 20
  result <- run_sim(direction, start, burnin, gen_time, sim_length); result
  expect_true(min(result$time) == start)
  expect_true(max(result$time) == start + gen_time * round(sim_length / gen_time))
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

test_that("Forward simulation from generation > 1 has the correct length with burnin", {
  direction <- "forward"; start <- 200; sim_length <- 50; burnin <- 100; gen_time <- 20
  result <- run_sim(direction, start, burnin, gen_time, sim_length); result
  expect_true(min(result$time) == start)
  expect_true(max(result$time) == start + gen_time * round(sim_length / gen_time))
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

# backward simulations - generation time > 1 ------------------------------

test_that("Backward simulation has the correct length without burnin", {
  direction <- "backward"; start <- 200; burnin <- 0; gen_time <- 20
  result <- run_sim(direction, start, burnin, gen_time)
  expect_true(max(result$time) == start)
  expect_true(min(result$time) == 0)
  expect_true(length(unique(result$time)) == round(start / gen_time) + 1)
})

test_that("Backward simulation has the correct length with burnin", {
  direction <- "backward"; start <- 200; burnin <- 300; gen_time <- 20
  result <- run_sim(direction, start, burnin, gen_time)
  expect_true(max(result$time) == start)
  expect_true(min(result$time) == 0)
  expect_true(length(unique(result$time)) == round(start / gen_time) + 1)
})

test_that("Backward simulation of limited length has the correct length without burnin", {
  direction <- "backward"; start <- 200; sim_length = 50; burnin <- 0; gen_time <- 20
  result <- run_sim(direction, start, burnin, gen_time, sim_length)
  expect_true(max(result$time) == start)
  expect_true(min(result$time) == start - gen_time * round(sim_length / gen_time))
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

test_that("Backward simulation of limited length has the correct length with burnin", {
  direction <- "backward"; start <- 5; sim_length = 3; burnin <- 20; gen_time <- 1
  result <- run_sim(direction, start, burnin, gen_time, sim_length)
  expect_true(max(result$time) == start)
  expect_true(min(result$time) == start - gen_time * round(sim_length / gen_time))
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})


# backward simulations - generation time > 1 (times not multiples) --------

test_that("Backward simulation has the correct length without burnin", {
  direction <- "backward"; start <- 200; burnin <- 0; gen_time <- 30
  result <- run_sim(direction, start, burnin, gen_time)
  shifted_start <- round(start / gen_time) * gen_time
  expect_true(max(result$time) == shifted_start)
  expect_true(min(result$time) == 0)
  expect_true(length(unique(result$time)) == round(start / gen_time) + 1)
})

test_that("Backward simulation has the correct length with burnin", {
  direction <- "backward"; start <- 200; burnin <- 300; gen_time <- 30
  result <- run_sim(direction, start, burnin, gen_time)
  shifted_start <- round(start / gen_time) * gen_time
  expect_true(max(result$time) == shifted_start)
  expect_true(min(result$time) == 0)
  expect_true(length(unique(result$time)) == round(start / gen_time) + 1)
})

test_that("Backward simulation of limited length has the correct length without burnin", {
  direction <- "backward"; start <- 200; sim_length = 50; burnin <- 0; gen_time <- 30
  result <- run_sim(direction, start, burnin, gen_time, sim_length)
  shifted_start <- round(start / gen_time) * gen_time
  expect_true(max(result$time) == shifted_start)
  expect_true(min(result$time) == shifted_start - gen_time * round(sim_length / gen_time))
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})

test_that("Backward simulation of limited length has the correct length with burnin", {
  direction <- "backward"; start <- 5; sim_length = 3; burnin <- 20; gen_time <- 30
  result <- run_sim(direction, start, burnin, gen_time, sim_length)
  shifted_start <- round(start / gen_time) * gen_time
  expect_true(max(result$time) == shifted_start)
  expect_true(min(result$time) == shifted_start - gen_time * round(sim_length / gen_time))
  expect_true(length(unique(result$time)) == round(sim_length / gen_time) + 1)
})
