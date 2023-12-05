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

# Input argument checks --------------------------------------------------------
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

test_that("fi_pool_cluster() s and N are numeric/integer vectors of same length", {
  e <- "pool_size and pool_number must be vectors of positive numbers of common length. pool_size can be non-integer, but pool_number must be integer"
  expect_error(fi_pool_cluster(
    pool_size = 10, pool_number = c(10, 20, 30), prevalence = 0.01,
    correlation = 0.05, sensitivity = 0.95, specificity = 0.99), e)
  expect_error(fi_pool_cluster(
    pool_size = "foo", pool_number = "bar", prevalence = 0.01,
    correlation = 0.05, sensitivity = 0.95, specificity = 0.99), e)
  expect_error(fi_pool_cluster(
    pool_size = c(1.1, 1), pool_number = c(1.1, 1), prevalence = 0.01,
    correlation = 0.05, sensitivity = 0.95, specificity = 0.99), e)
  expect_error(fi_pool_cluster(
    pool_size = 0, pool_number = 10, prevalence = 0.01,
    correlation = 0.05, sensitivity = 0.95, specificity = 0.99), e)
})

test_that("fi_pool_cluster() when K == 1 && N == 1 && s == 1", {
  expect_equal(fi_pool_cluster(
    pool_size = 1, pool_number = 1, prevalence = 0.01, correlation = 0.05,
    sensitivity = 0.95, specificity = 0.99),
    matrix(c(46.44747, 0, 0, 0), nrow = 2), tolerance = 1e-5)
})

test_that("fi_pool_cluster() when form == 'discrete'", {
  expect_equal(fi_pool_cluster(
    pool_size = 10, pool_number = 5, prevalence = 0.01, correlation = 0.05,
    sensitivity = 0.95, specificity = 0.99, form = "discrete"),
    matrix(c(3793.83363, -39.6160580, -39.6160580, 0.6903514), nrow = 2),
    tolerance = 1e-6)
})

test_that("fi_pool_cluster() when real_scale", {
  expect_equal(fi_pool_cluster(
    pool_size = 10, pool_number = 5, prevalence = 0.01, correlation = 0.05,
    sensitivity = 0.95, specificity = 0.99, form = "cloglognorm", real_scale = T),
    matrix(c(0.1536556, 0.1461399, 0.1461399, 0.2246681), nrow = 2),
    tolerance = 1e-6)
})

test_that("fi_pool_cluster() when rho 0 or 1", {
  expect_equal(fi_pool_cluster(
    pool_size = 10, pool_number = 5, prevalence = 0.01, correlation = 0,
    sensitivity = 0.95, specificity = 0.99),
    4100.879, tolerance = 1e-3)
  expect_error(fi_pool_cluster(
    pool_size = 10, pool_number = 5, prevalence = 0.01, correlation = 1,
    sensitivity = 0.95, specificity = 0.99),
    "correlation = 1 \\(perfect correlation\\) means that units from the same location are prefectly correlated. Do not use cluster surveys in this case.")
})

test_that("fi_pool_cluster() when s != 1 and form = discrete", {
  expect_error(fi_pool_cluster(
    pool_size = c(1, 2), pool_number = c(5, 10), prevalence = 0.01, correlation = 0.05,
    sensitivity = 0.95, specificity = 0.99, form = "discrete", real_scale = T),
    "unequal pool size not implemented for form = 'discrete'")
})

# Integration fails ------------------------------------------------------------
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
