test_that("non-overlapping geneflow leads to error", {
  map <- readRDS("map.rds")
  pop1 <- population("pop1", parent = "ancestor", N = 100, time = 10000,
                      center = c(0, 10), radius = 1, map = map, intersect = FALSE)
  pop2 <- population("pop2", parent = "ancestor", N = 100, time = 10000,
                     center = c(0, -10), radius = 1, map = map, intersect = FALSE)
  expect_error(geneflow(from = pop1, to = pop2, start = 1000, end = 0, rate = 0.1),
               "No overlap between population ranges of pop1 and pop2 at time 1000.")
})

test_that("non-overlapping geneflow passes if the check is explicitly turned off", {
  map <- readRDS("map.rds")
  pop1 <- population("pop1", parent = "ancestor", N = 100, time = 10000,
                     center = c(0, 10), radius = 1, map = map, intersect = FALSE)
  pop2 <- population("pop2", parent = "ancestor", N = 100, time = 10000,
                     center = c(0, -10), radius = 1, map = map, intersect = FALSE)
  expect_s3_class(
    geneflow(from = pop1, to = pop2, start = 1000, end = 0, rate = 0.1, overlap = FALSE),
    "data.frame"
  )
})

test_that("populations must be present for them to mix", {
  map <- readRDS("map.rds")
  pop1 <- population("pop1", parent = "ancestor", N = 100, time = 100,
                      center = c(0, 10), radius = 1, map = map, intersect = FALSE)
  pop2 <- population("pop2", parent = "ancestor", N = 100, time = 100,
                     center = c(0, -10), radius = 1, map = map, intersect = FALSE)
  expect_error(geneflow(from = pop1, to = pop2, start = 1000, end = 0, rate = 0.1),
               "Specified times are not consistent with the assumed direction of\ntime (geneflow pop1 -> pop2 in the time window 1000-0)", fixed = TRUE)
})
