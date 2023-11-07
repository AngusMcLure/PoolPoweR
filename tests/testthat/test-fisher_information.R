# Tolerance to address floating point precision errors
test_that("fi_pool() works with expected ranges of pool_size and prevalence", {
  expect_equal(fi_pool(pool_size = 10, prevalence = 0.85), 2.562891e-05, tolerance = 1e-6)
  expect_equal(fi_pool(pool_size = 5, prevalence = 0.85), 0.08438141, tolerance = 1e-6)
  expect_equal(fi_pool(pool_size = 5, prevalence = 0.65), 1.077534, tolerance = 1e-6)
  # Should NaN be returned when prevalence == 1?
  expect_true(is.nan(fi_pool(pool_size = 10, prevalence = 1)))
})

test_that("fi_pool_imperfect() works with expected input ranges", {
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
})

test_that("fi_ratio() works with expected input ranges", {
  expect_equal(
    fi_ratio(
      s = 10, p = 0.9, sensitivity = 1, specificity = 1, perfect_reference = T
    ),
    111111111
  )
  expect_equal(
    fi_ratio(
      s = 10, p = 0.7, sensitivity = 0.85, specificity = 0.9, perfect_reference = T
    ),
    278609767
  )
  expect_equal(
    fi_ratio(
      s = 10, p = 0.7, sensitivity = 0.85, specificity = 0.9, perfect_reference = F
    ),
    140419323
  )
})
