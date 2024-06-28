test_that("check_geq2", {
  x <- 1
  expect_silent(check_geq2(x, 1))
  expect_error(check_geq2(x, 2), "x must be >= 2.")
  y <- "a"
  expect_error(check_geq2(y, 1), "y must be numeric, not character.")
})

test_that("check_in_range2", {
  x <- 0.01
  expect_silent(check_in_range2(x))
  x <- -1
  expect_error(check_in_range2(x), "x = -1")
  x <- 2
  expect_error(check_in_range2(x), "x = 2")
  y <- "a"
  expect_error(check_in_range2(y), "y must be numeric, not character.")
})