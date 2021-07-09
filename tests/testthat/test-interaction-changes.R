map <- readRDS("map.rds")

pop <- population("pop", time = 1000, N = 10, map = map, center = c(0, 0), radius = 500e3) %>%
  move(trajectory = c(10, 10), start = 900, end = 700, snapshots = 3)

test_that("temporal consistency of interaction parameter changes is enforced", {
  expect_error(dispersal(pop, time = 950, competition_dist = 100),
               "The new event (.*) pre-dates the last specified active event (.*)")
  expect_silent(dispersal(pop, time = 50, competition_dist = 100))
})

test_that("at least one interaction parameter is specified", {
  expect_error(dispersal(pop, time = 1000),
               "At least one spatial interaction parameter must be specified")
})

test_that("interaction parameter must be positive, non-zero values", {
  msg <- "Spatial interaction parameters can only have positive"
  expect_error(dispersal(pop, time = 1000, competition_dist = -100), msg)
  expect_error(dispersal(pop, time = 1000, mate_dist = -100), msg)
  expect_error(dispersal(pop, time = 1000, dispersal_dist = -100), msg)
})

test_that("interaction parameter change is correctly recorded", {
  x1 <- dispersal(pop, time = 100, competition_dist = 100)
  x2 <- dispersal(pop, time = 100, mate_dist = 50)
  x3 <- dispersal(pop, time = 100, dispersal_dist = 20)
  x4 <- dispersal(pop, time = 100, competition_dist = 50, dispersal_dist = 10)

  hist1 <- attr(x1, "history") %>% .[[length(.)]]
  expect_true(hist1$pop == pop$pop[1])
  expect_true(hist1$time == 100)
  expect_true(hist1$event == "dispersal")
  expect_true(hist1$competition_dist == 100)
  expect_true(is.na(hist1$mate_dist))
  expect_true(is.na(hist1$dispersal_dist))

  hist2 <- attr(x2, "history") %>% .[[length(.)]]
  expect_true(hist2$mate_dist == 50)

  hist3 <- attr(x3, "history") %>% .[[length(.)]]
  expect_true(hist3$dispersal_dist == 20)

  hist4 <- attr(x4, "history") %>% .[[length(.)]]
  expect_true(hist4$competition_dist == 50 && hist4$dispersal_dist == 10)
})

