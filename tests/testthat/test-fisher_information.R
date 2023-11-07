# Tolerance to address floating point precision errors
test_that("fi_pool_imperfect() works with expected input ranges", {
  # Tests mainly to ensure same outputs when refactoring fi_pool internals
  expect_equal(
    fi_pool_imperfect(
      pool_size = 10,
      prevalence = 0.7,
      sensitivity = 0.8,
      specificity = 0.9
    ),
    1.186457e-07,
    tolerance = 1e-6
  )
  expect_equal(
    fi_pool_imperfect(
      pool_size = 20,
      prevalence = 0.55,
      sensitivity = 0.9,
      specificity = 0.6
    ),
    7.376202e-11,
    tolerance = 1e-6
  )
  expect_true(is.nan(fi_pool_imperfect(
    pool_size = 10,
    prevalence = 1,
    sensitivity = 1,
    specificity = 1
  )))
})
