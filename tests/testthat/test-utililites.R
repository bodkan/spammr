test_that("distances beyond the world dimension throw and error", {
  map <- readRDS("map.rds")
  xrange <- diff(sf::st_bbox(map)[c("xmin", "xmax")])
  yrange <- diff(sf::st_bbox(map)[c("ymin", "ymax")])
  error_msg <- "larger than the overall world size"
  expect_error(check_resolution(map, xrange * 10), error_msg)
  expect_error(check_resolution(map, yrange * 10), error_msg)
})

test_that("distances beyond the world dimension throw and error", {
  map <- readRDS("map.rds")
  xrange <- diff(sf::st_bbox(map)[c("xmin", "xmax")])
  yrange <- diff(sf::st_bbox(map)[c("ymin", "ymax")])
  expect_silent(check_resolution(map, xrange / 10))
  expect_silent(check_resolution(map, yrange / 10))
})

test_that("binaries are found for all methods of running SLiM", {
  skip_on_os("windows")
  skip_on_os("linux")
  expect_equal(get_binary("gui"), "open -a SLiMgui")
  expect_equal(get_binary("batch"), as.character(Sys.which("slim")))
})
