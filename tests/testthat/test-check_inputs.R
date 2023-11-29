# Not using cli::test_that_cli for cleaner output
test_that("check_geq()", {
  expect_silent(check_geq("max_s", 1)) # > 1
  # > 0
  expect_silent(check_geq("pool_size", 1))
  expect_snapshot(check_geq("pool_size", "chr"))
  expect_snapshot(check_geq("pool_size", -1))
  expect_snapshot(check_geq("max_s", 0))
  expect_error(check_geq("prevalence", 0.05),
               "Needs to be one of the accepted_args")
})

test_that("check_in_range()", {
  expect_silent(check_in_range("prevalence", 0.05))
  expect_snapshot(check_in_range("prevalence", -1))
  expect_snapshot(check_in_range("prevalence", 1.1))
  expect_error(check_in_range("pool_size", 1.1),
               "Needs to be one of the accepted_args")
  expect_silent(check_in_range("sensitivity", 1))
  expect_silent(check_in_range("specificity", 0))
})
