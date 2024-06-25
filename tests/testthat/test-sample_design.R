# fixtures ----
fixed_perfect <- fixed_design(
  pool_size = 10, pool_number = NULL, sensitivity = 1, specificity = 1
)

fixed_null <- fixed_design() # sens/spec == 1, pool_size/num == NULL

# fixed_design ----
test_that("fixed_design constructor", {
  expect_equal(class(fixed_perfect), c("fixed_design", "sample_design"))
  expect_equal(fixed_perfect$pool_size, 10)
  expect_equal(fixed_perfect$pool_number, NULL)
  expect_equal(fixed_perfect$sensitivity, 1)
  expect_equal(fixed_perfect$specificity, 1)
})

test_that("fixed_design default", {
  expect_equal(class(fixed_perfect), c("fixed_design", "sample_design"))
  expect_equal(fixed_null$pool_size, NULL)
  expect_equal(fixed_null$pool_number, NULL)
  expect_equal(fixed_null$sensitivity, 1)
  expect_equal(fixed_null$specificity, 1)
}) 

test_that("fixed_design bad inputs caught", {
  expect_error(fixed_design(pool_size = -1))
  expect_error(fixed_design(pool_number = -1))
})
