# Tolerance to address floating point precision errors
test_that("fi_pool() works with expected input ranges", {
  # Tests mainly to ensure same outputs when refactoring fi_pool internals
  # This one has reasonable params
  expect_equal(fi_pool(
    pool_size = 10, prevalence = 0.01, sensitivity = 0.95, specificity = 0.99),
    820.1759, tolerance = 1e-4
  )
  # These ones do not really
  expect_equal(
    fi_pool(
      pool_size = 10,
      prevalence = 0.7,
      sensitivity = 0.8,
      specificity = 0.9
    ),
    1.186457e-07,
    tolerance = 1e-6
  )
  expect_equal(
    fi_pool(
      pool_size = 20,
      prevalence = 0.55,
      sensitivity = 0.9,
      specificity = 0.6
    ),
    7.376202e-11,
    tolerance = 1e-6
  )
  expect_true(is.nan(fi_pool(
    pool_size = 10,
    prevalence = 1,
    sensitivity = 1,
    specificity = 1
  )))
})

test_that("fi_pool_cluster() outputs a 2x2 matrix for input vectors of length 2", {
  ## Reasonable params
  expect_true(all.equal(fi_pool_cluster(
    pool_size = 10, pool_number = 5, prevalence = 0.01, correlation = 0.05,
    sensitivity = 0.95, specificity = 0.99),
    matrix(c(1880.3484, -125.47514, -125.4751, 23.71574), nrow = 2), tolerance = 1e-5
    ))
  expect_true(all.equal(fi_pool_cluster(
    pool_size = c(1, 2), pool_number = c(5, 10), prevalence = 0.01, correlation = 0.05,
    sensitivity = 0.95, specificity = 0.99),
    matrix(c(926.41807, -23.055960, -23.055960, 9.535592), nrow = 2), tolerance = 1e-6
    ))
})

test_that("fi_pool_cluster() fails when integral is divergent", {
  expect_error(
    fi_pool_cluster(
      pool_size = c(10, 20),
      pool_number = c(10, 100),
      prevalence = 0.9,
      sensitivity = 1,
      specificity = 1,
      correlation = 0.75,
      form = "beta",
      real_scale = FALSE
    ),
    "the integral is probably divergent"
  )
  expect_error(
    fi_pool_cluster(
      pool_size = c(10, 20),
      pool_number = c(10, 100),
      prevalence = 0.9,
      sensitivity = 1,
      specificity = 1,
      correlation = 0.5,
      form = "beta",
      real_scale = FALSE
    ),
    "the integral is probably divergent"
  )
})

test_that("fi_pool_cluster() fails when likelihoods or derivatives do not add up", {
  # SLOW (doesn't scale with pool size, number, and vector length)
  expect_error(
    fi_pool_cluster(
      pool_size = c(5, 10),
      pool_number = c(10, 50),
      prevalence = 0.5,
      sensitivity = 1,
      specificity = 1,
      correlation = 0.1,
      form = "beta",
      real_scale = FALSE
    ),
    "Error in integration of likelihoods. Likelihoods do not add to 1 or derivatives of likelihood with respect to parameters do not sum to 0"
  )
})

test_that("bad real_scale input caught in fi_pool_cluster()", {
  expect_error(
    fi_pool_cluster(
      pool_size = 5, pool_number = 10, prevalence = 0.01, correlation = 0.1,
      sensitivity = 1, specificity = 1, real_scale = 0),
    "0 is not TRUE/FALSE")
  expect_error(
    fi_pool_cluster(
      pool_size = 5, pool_number = 10, prevalence = 0.01, correlation = 0.1,
      sensitivity = 1, specificity = 1, real_scale = "yes"),
    "yes is not TRUE/FALSE")
})
