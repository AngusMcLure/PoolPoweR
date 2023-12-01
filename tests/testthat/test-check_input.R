test_that("check_geq()", {
  expect_silent(check_geq("max_s", 1)) # > 1
  # > 0
  expect_silent(check_geq("pool_size", 1))
  expect_error(check_geq("pool_size", "chr"), "chr is a character.")
  expect_error(check_geq("pool_size", -1), "-1 is < 0")
  expect_error(check_geq("max_s", 0), "0 is < 1")
  expect_error(check_geq("prevalence", 0.05),
               "Needs to be one of the accepted_args")
})

test_that("check_in_range()", {
  expect_silent(check_in_range("prevalence", 0.05))
  expect_error(check_in_range("prevalence", -1), "-1 is < 0")
  expect_error(check_in_range("prevalence", 1.1), "1.1 is > 1")
  expect_error(check_in_range("pool_size", 1.1),
               "Needs to be one of the accepted_args")
  expect_silent(check_in_range("sensitivity", 1))
  expect_silent(check_in_range("specificity", 0))
})

test_that("check_rho()", {
  expect_silent(check_rho(0))
  expect_silent(check_rho(1))
  expect_silent(check_rho(NA))
  expect_error(check_rho(-1), "-1 is < 0")
  expect_error(check_rho(2), "2 is > 1")
  expect_error(check_rho("chr"), "chr is a character.")
})

test_that("check_form()", {
  e <- "form must be one of 'beta', 'logitnorm', 'cloglognorm', or 'discrete'."
  expect_silent(check_form("beta"))
  expect_error(check_form("binomial"), e)
  expect_error(check_form(1), e)
})

test_that("check_scale()", {
  expect_silent(check_scale(T))
  expect_error(check_scale("chr"), "chr is not TRUE/FALSE")
  expect_error(check_scale(10), "10 is not TRUE/FALSE")
})
