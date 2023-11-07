# Tolerance to address floating point precision errors
test_that("fi_pool() works with expected input ranges", {
  # Tests mainly to ensure same outputs when refactoring fi_pool internals
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
  expect_true(all.equal(
    fi_pool_cluster(
      pool_size = c(10, 20),
      pool_number = c(10, 20),
      prevalence = 0.9,
      sensitivity = 1,
      specificity = 1,
      correlation = 0.5,
      form = "beta",
      real.scale = TRUE
    ),
    matrix(c(5.431112, -1.928655, -1.928655, 1.123625), nrow = 2),
    tolerance = 1e-6
  ))
  expect_true(all.equal(
    fi_pool_cluster(
      pool_size = c(5, 10),
      pool_number = c(10, 15),
      prevalence = 0.9,
      sensitivity = 1,
      specificity = 1,
      correlation = 0.5,
      form = "beta",
      real.scale = TRUE
    ),
    matrix(c(9.028451, -1.718530, -1.718530, 1.130158), nrow = 2),
    tolerance = 1e-6
  ))
  expect_true(all.equal(
    fi_pool_cluster(
      pool_size = c(5, 10),
      pool_number = c(10, 15),
      prevalence = 0.7,
      sensitivity = 0.95,
      specificity = 0.8,
      correlation = 0.9,
      form = "cloglognorm",
      real.scale = FALSE
    ),
    matrix(c(4.692347, -2.270677, -2.270677, 8.494053), nrow = 2),
    tolerance = 1e-6
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
      real.scale = FALSE
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
      real.scale = FALSE
    ),
    "the integral is probably divergent"
  )
})

test_that("fi_pool_cluster() fails when likelihoods or derivatives do not add up", {
  expect_error(
    fi_pool_cluster(
      pool_size = c(5, 10),
      pool_number = c(10, 50),
      prevalence = 0.5,
      sensitivity = 1,
      specificity = 1,
      correlation = 0.1,
      form = "beta",
      real.scale = FALSE
    ),
    "Error in integration of likelihoods. Likelihoods do not add to 1 or derivatives of likelihood with respect to parameters do not sum to 0"
  )
})
