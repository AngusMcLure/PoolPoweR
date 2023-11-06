test_that("fi_pool() works with expected ranges of pool_size and prevalence", {
  # Some nasty floating point errors
  expect_equal(fi_pool(pool_size = 10, prevalence = 0.85), 2.562891e-05, tolerance = 1e-6)
  expect_equal(fi_pool(pool_size = 5, prevalence = 0.85), 0.08438141, tolerance = 1e-6)
  expect_equal(fi_pool(pool_size = 5, prevalence = 0.65), 1.077534, tolerance = 1e-6)
})
