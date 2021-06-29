test_that("anything but the three allowed landsccape sources leads to error", {
  error_msg <- "Landscape has to be either 'blank', 'naturalearth' or an object of the class 'sf'"
  expect_error(world(xrange = c(-15, 60), yrange = c(20, 65), landscape = "blank123"), error_msg, fixed = TRUE)
  expect_error(world(xrange = c(-15, 60), yrange = c(20, 65), landscape = 123), error_msg, fixed = TRUE)
})

test_that("blank abstract landscape is allowed", {
  expect_silent(world(xrange = c(-15, 60), yrange = c(20, 65), landscape = "blank"))
})

test_that("user-defined abstract landscape is allowed", {
  xrange <- c(-15, 60)
  yrange <- c(20, 65)
  polygon <- slendr:::create_polygon(list(c(-10, 30), c(50, 30), c(40, 50), c(0, 40)))
  expect_silent(world(xrange = xrange, yrange = yrange, landscape = polygon))
})
