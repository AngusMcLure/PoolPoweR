# fixtures ----
fd <- fixed_design(
  pool_size = 5, pool_number = 10, sensitivity = 0.99, specificity = 0.95
)

# fixed_design ----
test_that("design_effect() gives consistent output for basic tests", {
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

### design_effect_random() ----

