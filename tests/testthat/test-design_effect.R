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
  expect_equal(
    design_effect(
      pool_size = 100,
      pool_number = 10,
      prevalence = 0.8,
      sensitivity = 1,
      specificity = 0.7,
      correlation = 0.8,
      form = "cloglognorm"
    ),
    4422,
    tolerance = 1e-4
  )
})
