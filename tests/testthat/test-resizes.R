map <- readRDS("map.rds")

run_sim <- function(pop, direction, sim_length = NULL) {
  model_dir <- "/tmp/test_resize"
  unlink(model_dir, recursive = TRUE, force = TRUE)

  model <- compile(
    populations = list(pop), generation_time = 1,
    resolution = 10e3, competition_dist = 130e3, mate_dist = 100e3, dispersal_dist = 70e3, # how far will offspring end up from their parents
    dir = model_dir, direction = direction, sim_length = sim_length
  )

  slim(model, seq_length = 1, recomb_rate = 0, save_locations = TRUE, method = "batch")

  dt <- fread(file.path(model$directory, "output_ind_locations.tsv.gz"))

  dt[, .N, by = .(gen, time, pop)]
}

calculate_exp_sizes <- function(N1, N2, t1, t2) {
  r <- log(N2 / N1) / abs(t2 - t1)
  N <- rev(round(N1 * exp(r * (abs(t2 - t1) - seq(0, abs(t2 - t1))))))
  N
}

test_that("Population size stays constant when specified", {
  pop <- population("pop", time = 1000, N = 100, map = map, center = c(20, 50), radius = 500e3)
  res <- run_sim(pop, "backward")
  start_N <- attr(pop, "history")[[1]]$N
  expect_true(all(res$N == start_N))
})

test_that("Population size is increased correctly in single step", {
  pop <- population("pop", time = 1000, N = 100, map = map, center = c(20, 50), radius = 500e3)
  res <- run_sim(resize(pop, N = 50, time = 900, how = "step"), "backward")
  start_N <- attr(pop, "history")[[1]]$N
  expect_true(res[time > 900, all(N == start_N)])
  expect_true(res[time <= 900, all(N == 50)])
})

test_that("Population size is increased correctly in two steps", {
  pop <- population("pop", time = 1000, N = 100, map = map, center = c(20, 50), radius = 500e3)
  res <- run_sim(
    resize(pop, N = 50, time = 900, how = "step") %>%
      resize(pop, N = 500, time = 300, how = "step"), "backward")
  start_N <- attr(pop, "history")[[1]]$N
  expect_true(res[time > 900, all(N == start_N)])
  expect_true(res[time <= 900 & time > 300, all(N == 50)])
  expect_true(res[time < 300, all(N == 500)])
})

test_that("Simulated exponential growth matches theoretical expectations", {
  N1 <- 100; N2 <- 1000
  t1 <- 900; t2 <- 600

  pop <- population("pop", time = 1000, N = N1, map = map, center = c(20, 50), radius = 500e3)
  res <- run_sim(resize(pop, N = N2, time = t1, end = t2, how = "exponential"), "backward")
  expected_N <- calculate_exp_sizes(N1, N2, t1, t2)

  expect_true(res[1, N == N1]) # initial size
  expect_true(res[N > N1, time[1] == t1 - 1]) # start of the growth
  expect_true(res[.N, N == N2]) # final size

    # sizes match expectations
  expect_true(res[time <= t1 & time >= t2, all(N == expected_N)])
})

test_that("Simulated exponential shrinking matches theoretical expectations", {
  N1 <- 1000; N2 <- 300
  t1 <- 900; t2 <- 600

  pop <- population("pop", time = 1000, N = N1, map = map, center = c(20, 50), radius = 500e3)
  res <- run_sim(resize(pop, N = N2, time = t1, end = t2, how = "exponential"), "backward")
  expected_N <- calculate_exp_sizes(N1, N2, t1, t2)

  expect_true(res[1, N == N1]) # initial size
  expect_true(res[N < N1, time[1] == t1 - 1]) # start of the growth
  expect_true(res[.N, N == N2]) # final size

  # sizes match expectations
  expect_true(res[time <= t1 & time >= t2, all(N == expected_N)])
})

test_that("Multiple resize event types are allowed (backward model)", {
  N1 <- 100; N2 <- 1000; N3 <- 80
  t1 <- 900; t2 <- 600; t3 <- 200

  pop <- population("pop", time = 1000, N = N1, map = map, center = c(20, 50), radius = 500e3) %>%
    resize(N = N2, time = t1, end = t2, how = "exponential") %>%
    resize(N = N3, time = t3, how = "step")

  res <- run_sim(pop, "backward")
  expected_N <- calculate_exp_sizes(N1, N2, t1, t2)

  expect_true(res[1, N == N1]) # initial size
  expect_true(res[N > N1, time[1] == t1 - 1]) # start of the exponential growth
  expect_true(res[N == N3, time[1] == t3]) # time of the step resize event
  expect_true(res[.N, N == N3]) # final size

  # sizes match expectations
  expect_true(res[time <= t1 & time >= t2, all(N == expected_N)])
})

test_that("Multiple resize event types are allowed (forward model)", {
  N1 <- 1000; N2 <- 300; N3 <- 80
  t1 <- 2; t2 <- 600; t3 <- 900

  pop <- population("pop", time = 1, N = N1, map = map, center = c(20, 50), radius = 500e3) %>%
    resize(N = N2, time = t1, end = t2, how = "exponential") %>%
    resize(N = N3, time = t3, how = "step")

  res <- run_sim(pop, "forward", sim_length = 1000)
  expected_N <- calculate_exp_sizes(N1, N2, t1, t2)

  expect_true(res[1, N == N1]) # initial size
  expect_true(res[N < N1, time[1] == t1 + 1]) # start of the exponential growth
  expect_true(res[N == N3, time[1] == t3]) # time of the step resize event
  expect_true(res[.N, N == N3]) # final size

  # sizes match expectations
  expect_true(res[time <= t1 & time >= t2, all(N == expected_N)])
})
