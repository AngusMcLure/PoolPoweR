# Not using cli::test_that_cli for cleaner output
test_that("check_geq()", {
  expect_silent(check_geq("max_s", 1)) # > 1
  # > 0
  expect_silent(check_geq("pool_size", 1))
  expect_snapshot(check_geq("pool_size", "chr"))
  expect_snapshot(check_geq("pool_size", -1))
  expect_snapshot(check_geq("max_s", 0))
})
