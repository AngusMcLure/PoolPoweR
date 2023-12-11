### design_effect() -----------------------------------------------------------
test_that("design_effect() gives consistent output for basic tests", {
  # This one has reasonable inputs
  expect_equal(
    design_effect(
      pool_size = 5,
      pool_number = 10,
      prevalence = 0.01,
      correlation = 0.05,
      sensitivity = 0.99,
      specificity = 0.95),
    0.7240988, tolerance = 1e-7
  )
  expect_equal(
    design_effect(
      pool_size = 10,
      pool_number = 10,
      prevalence = 0.9,
      sensitivity = 1,
      specificity = 1,
      correlation = 0.9,
      form = "beta"
    ),
    118.9243,
    tolerance = 1e-4
  )
  expect_equal(
    design_effect(
      pool_size = 10,
      pool_number = 10,
      prevalence = 0.2,
      sensitivity = 0.9,
      specificity = 0.8,
      correlation = 0.2,
      form = "beta"
    ),
    26.50055,
    tolerance = 1e-5
  )
})

test_that("design_effect() fails for some very unusual parameters because integral in call to fi_pool_cluster() appears starts to have numeric issues", {
  expect_error(
    design_effect(
      pool_size = 100,
      pool_number = 10,
      prevalence = 0.8,
      sensitivity = 1,
      specificity = 0.7,
      correlation = 0.8,
      form = "cloglognorm"
    ),
    "Error in integration of likelihoods. Likelihoods do not add to 1 or derivatives of likelihood with respect to parameters do not sum to 0"
  )
})


test_that("bad inputs caught in design_effect()", {
  expect_error(design_effect(pool_size = TRUE), "TRUE is a logical")
  expect_error(design_effect(pool_size = 5, pool_number = 10, prevalence = 10), "10 is > 1")
  expect_error(design_effect(pool_size = 5, pool_number = 10, prevalence = 0.01, correlation = 0.1, sensitivity = 1, specificity = 1, form = "binomal"),
               "form must be one of 'beta', 'logitnorm', 'cloglognorm', or 'discrete'.")
})

### design_effect_random() ----------------------------------------------------

