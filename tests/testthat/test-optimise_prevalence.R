test_that("optimise_sN_prevalence() gives consistent output for basic tests", {
  # Reasonable params
  expect_true(all.equal(
    optimise_sN_prevalence(
      prevalence = 0.01, cost_unit = 5, cost_pool = 10,
      cost_cluster = 100, correlation = 0.05
    ), list(s=5,cost =0.2513798, catch=20, N=4),tolerance=1e-7))
  # Rest not so much
  expect_true(all.equal(
    optimise_sN_prevalence(
      prevalence = 0.9,
      cost_unit = 1,
      cost_pool = 10,
      cost_cluster = 20,
      correlation = 0.9,
      form = "beta",
      sensitivity = 1,
      specificity = 1,
      max_s = 50,
      max_N = 20
    ),
    list(s = 1, cost = 3.591, catch = 2, N = 2)
  ))
  expect_true(all.equal(
    optimise_sN_prevalence(
      prevalence = 0.5,
      cost_unit = 10,
      cost_pool = 100,
      cost_cluster = 2,
      correlation = 0.2,
      form = "beta",
      sensitivity = 1,
      specificity = 1,
      max_s = 10,
      max_N = 10
    ),
    list(s = 1, cost = 33.3, catch = 2, N = 2)
  ))
})

test_that("optimise_sN_prevalence() when correlation == 0", {
  expect_true(all.equal(
    optimise_sN_prevalence(
      prevalence = 0.9,
      cost_unit = 1,
      cost_pool = 10,
      cost_cluster = 20,
      correlation = 0,
      form = "beta",
      sensitivity = 1,
      specificity = 1,
      max_s = 50,
      max_N = 20
    ),
    list(s = 1, cost = 0.99, catch = Inf, N = Inf)
  ))
})

test_that("optimise_sN_prevalence() when opt$N == max_N", {
  expect_warning(
    optimise_sN_prevalence(
      prevalence = 0.01, cost_unit = 5, cost_pool = 10,
      cost_cluster = 100, correlation = 0.05, max_N = 4),
    "Maximum cost effectivness is achieved at or above the maximum number of pools allowed. Consider increasing max_N")
})

# TODO: fix this
#test_that("optimise_sN_prevalence() when opt$s == max_s", {
#  # Not an amazing text, catches all the optimise_s_prevalence() warnings too
#  tryCatch(
#  optimise_sN_prevalence(
#    prevalence = 0.01, cost_unit = 5, cost_pool = 10,
#    cost_cluster = 100, correlation = 0.05, max_s = 5)
#  )
#  expect_warning()
#
#})

test_that("optimise_s_prevalence() gives consistent output for basic tests", {
  # Reasonable parameters
  expect_true(all.equal(
    optimise_s_prevalence(prevalence = 0.01, cost_unit = 5, cost_pool = 10),
    list(s=19, cost=0.05998076, catch=19), tolerance = 1e-7
  ))
  # Not very
  expect_true(all.equal(
    optimise_s_prevalence(
      prevalence = 0.7, cost_unit = 10, cost_pool = 100,
      cost_cluster = NA, correlation = NA, pool_number = 1, form = "beta", sensitivity = 1,
      specificity = 1, max_s = 50, interval = 0
    ),
    list(s = 1, cost = 23.1, catch = 1)
  ))
  expect_true(all.equal(
    optimise_s_prevalence(
      prevalence = 0.2, cost_unit = 1, cost_pool = 200,
      cost_cluster = NA, correlation = NA, pool_number = 10, form = "beta", sensitivity = 0.8,
      specificity = 0.9, max_s = 500, interval = 2
    ),
    list(
      s = 5, cost = 24.43913, catch = 50, s_interval = c(2, 11),
      cost_interval = c(36.73102, 59.40823), catch_interval = c(20, 110)
    ),
    tolerance = 1e-6
  ))
})

test_that("optimise_s_prevalence() throws error when form == 'discrete'", {
  expect_error(optimise_s_prevalence(
    prevalence = 0.01, cost_unit = 5, cost_pool = 10, form = "discrete"),
    'When form = "discrete" the cost of unit information function with
         respect to s often has mulitple minima and therefore the discrete
         distribution is not currently supported for optimisation')
})

test_that("optimise_s_prevalence() when cost_unit == Inf", {
  expect_true(all.equal(
    optimise_s_prevalence(prevalence = 0.01, cost_unit = Inf, cost_pool = 10, interval = 0.1),
    list(s=1, cost=NA, catch=1, s_interval=c(1,19), cost_interval=NA, catch_interval=c(1,19))
    ))
})

test_that("optimise_s_prevalence() when cost_pool == Inf", {
  expect_true(all.equal(
    optimise_s_prevalence(prevalence = 0.1, cost_unit = 5, cost_pool = Inf, interval = 0.1),
    list(s=15, cost=NA, catch=15, s_interval=c(10,21), cost_interval=NA, catch_interval=c(10,21))
  ))
})

test_that("optimise_s_prevalence() hits max_s when determining cost floor/ceiling", {
  expect_warning(
    optimise_s_prevalence(prevalence = 0.01, cost_unit = 5, cost_pool = Inf, max_s = 50),
    'Maximum cost effectivness is achieved at or above the maximum size of pools allowed. Consider increasing max_s')
})

test_that("optimise_s_prevalence() when cost(max_s) < max_cost", {
  expect_warning(
    optimise_s_prevalence(prevalence = 0.1, cost_pool = 1, cost_unit = 0.5, max_s = 6, interval = 0.1),
    'A pool size greater than max_s may fall within the specified range of cost effectiveness. Consider increasing max_s')
})

test_that("optimise_s_prevalence() extremely bad integrand behaviour", {
  # This takes a few seconds to run
  expect_error(
    optimise_s_prevalence(
      prevalence = 0.2, cost_unit = 1, cost_pool = 200,
      cost_cluster = 1, correlation = 0.6, pool_number = 10, form = "beta", sensitivity = 0.8,
      specificity = 0.9, max_s = 500, interval = 2
    ),
    "extremely bad integrand behaviour"
    )})
