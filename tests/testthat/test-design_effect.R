# fixed_design ----
fd <- fixed_design(
  pool_size = 5, pool_number = 10, sensitivity = 0.99, specificity = 0.95
)

test_that("fixed design_effect() gives consistent output for basic tests", {
  # This one has reasonable inputs
  expect_equal(
    design_effect(fd, prevalence = 0.01, correlation = 0.05),
    0.7240988, tolerance = 1e-7
  )
  expect_equal(
    design_effect(fixed_design(10, 10), prevalence = 0.9, correlation = 0.9),
    118.9243, tolerance = 1e-4
  )
  expect_equal(
    design_effect(
      fixed_design(10, 10, 0.9, 0.8), prevalence = 0.2, correlation = 0.2), 
      26.50055, tolerance = 1e-5
  )
})

test_that("design_effect() fails for some very unusual parameters because integral in call to fi_pool_cluster() appears starts to have numeric issues", {
  expect_error(
    design_effect(
      fixed_design(100, 10, 1, 0.7),
      prevalence = 0.8,
      correlation = 0.8,
      form = "cloglognorm"
    ),
    "Error in integration of likelihoods. Likelihoods do not add to 1 or derivatives of likelihood with respect to parameters do not sum to 0"
  )
})

### variable_design() ----
vd_target <- variable_design(nb_catch(5, 7), pool_target_number(10))
vd_max <- variable_design(nb_catch(5, 7), pool_max_size(10))

test_that("variable design_effect()", {
  expect_equal(
    design_effect(vd_target, prevalence = 0.01, correlation = 0.05),
    1.256262, tolerance = 1e-6
  )
  expect_equal(
    design_effect(vd_max, prevalence = 0.01, correlation = 0.05),
    3.726256, tolerance = 1e-6
  )
})

