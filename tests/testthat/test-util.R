# As util helper functions will be used to catch warnings soon
test_that("interval is >= 0", {
  expect_error(
    optimise_s_prevalence(prevalence = 0.01, cost_unit = 5, cost_pool = 10, interval = -1),
    'interval must be 0 or higher')
  expect_equal(
    optimise_s_prevalence(prevalence = 0.1, cost_unit = 5, cost_pool = 10, interval = 2),
    list(s=5, cost=0.786439, catch =5, s_interval=c(1,25),
         cost_interval=c(1.35, 2.262155), catch_interval=c(1,25)),
    tolerance = 1e-7)
})
