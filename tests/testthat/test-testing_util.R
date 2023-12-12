test_that("numeric_rel_diff()", {
  expect_true(numeric_rel_diff(0.0234491, 0.0234492, tolerance = 1e-5))
  expect_false(numeric_rel_diff(0.0234491, 0.0234492, tolerance = 1e-6))
})

test_that("list_rel_diff()", {
  act <- list(x=1, y=2, z=0.2513798)
  exp <- list(x=1, y=2, z=0.2514)
  expect_true(list_rel_diff(act, exp, 1e-4))
  expect_false(list_rel_diff(act, exp, 1e-5))
})

test_that("matrix_rel_diff() with rounding", {
  act <- matrix(c(0.3931327, -0.6625920, 0.8810527, -0.8209938), nrow = 2)
  exp <- matrix(c(0.3931, -0.6626, 0.8811, -0.8210), nrow = 2)
  expect_true(matrix_rel_diff(act, exp, tolerance = 1e-4))
  expect_false(matrix_rel_diff(act, exp, tolerance = 1e-5))

  # What happens when you have a big range
  act <- matrix(c(888, 0.0012345, 777, 222), nrow = 2)
  exp <- matrix(c(888, 0.001235, 777, 222), nrow = 2)
  expect_true(matrix_rel_diff(act, exp, tolerance = 1e-5))

  act <- matrix(c(888, 0.0012345, 777, 222), nrow = 2)
  exp <- matrix(c(900, 0.001235, 777, 222), nrow = 2)
  expect_true(matrix_rel_diff(act, exp, tolerance = 1e-2))
  # hmm
})
