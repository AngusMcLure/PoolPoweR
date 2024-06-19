perfect_design <- sample_design(
  pool_size = 10, pool_number = NULL, sensitivity = 1, specificity = 1
)

test_that("sample_design constructor", {
  expect_equal(class(perfect_design), "sample_design")
  expect_equal(perfect_design$pool_size, 10)
  expect_equal(perfect_design$pool_number, NULL)
  expect_equal(perfect_design$sensitivity, 1)
  expect_equal(perfect_design$specificity, 1)
})

test_that("", {
  act <- fi_pool(perfect_design, prevalence = 0.01)
  expect_equal(act, 965.0332, tolerance = 1e-4)
})
