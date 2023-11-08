test_that("optimise_sN_prevalence() gives consistent output for basic tests", {
  expect_true(all.equal(optimise_sN_prevalence(
    prevalence = 0.9,
    cost.unit = 1,
    cost.pool = 10,
    cost.location = 20,
    correlation = 0.9,
    form = 'beta',
    sensitivity = 1,
    specificity = 1,
    max.s = 50,
    max.N = 20),
    list(s=1, cost=3.591, catch=2, N=2)))
  expect_true(all.equal(optimise_sN_prevalence(
    prevalence = 0.5,
    cost.unit = 10,
    cost.pool = 100,
    cost.location = 2,
    correlation = 0.2,
    form = 'beta',
    sensitivity = 1,
    specificity = 1,
    max.s = 10,
    max.N = 10),
    list(s=1, cost=33.3, catch=2, N=2)))
})

test_that("optimise_s_prevalence() gives consistent output for basic tests", {
  expect_true(all.equal(optimise_s_prevalence(
    prevalence = 0.7, cost.unit = 10, cost.pool = 100,
    cost.location = NA, correlation = NA, N = 1, form = 'beta', sensitivity = 1,
    specificity = 1, max.s = 50, interval = 0),
    list(s=1, cost=23.1, catch=1)
  ))
  expect_true(all.equal(optimise_s_prevalence(
    prevalence = 0.2, cost.unit = 1, cost.pool = 200,
    cost.location = NA, correlation = NA, N = 10, form = 'beta', sensitivity = 0.8,
    specificity = 0.9, max.s = 500, interval = 2),
    list(s=5, cost=24.43913, catch=50, s_interval=c(2,11),
         cost_interval=c(36.73102, 59.40823), catch_interval=c(20, 110)
    ), tolerance = 1e-6
  ))
})
