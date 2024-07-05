# fixed_design ----
## fixtures ----
fixed_perfect <- fixed_design(
  pool_size = 10, pool_number = NULL, sensitivity = 1, specificity = 1
)

fixed_null <- fixed_design() # sens/spec == 1, pool_size/num == NULL

## test ----
test_that("fixed_design constructor", {
  expect_equal(length(fixed_perfect), 5)
  expect_equal(class(fixed_perfect), c("need_N", "fixed_design", "sample_design"))
  expect_equal(fixed_perfect$pool_size, 10)
  expect_equal(fixed_perfect$pool_number, NULL)
  expect_equal(fixed_perfect$sensitivity, 1)
  expect_equal(fixed_perfect$specificity, 1)
})

test_that("fixed_design default", {
  expect_equal(class(fixed_null), c("need_sN", "fixed_design", "sample_design"))
  expect_equal(fixed_null$pool_size, NULL)
  expect_equal(fixed_null$pool_number, NULL)
  expect_equal(fixed_null$sensitivity, 1)
  expect_equal(fixed_null$specificity, 1)
}) 


test_that("fixed_design bad inputs caught", {
  expect_error(fixed_design(pool_size = -1))
  expect_error(fixed_design(pool_number = -1))
})

# variable_design ----
## fixtures ----
var_target <- variable_design(
  catch_dist = nb_catch(5, 10), 
  pool_strat = pool_target_number(20)
)

var_max <- variable_design(
  catch_dist = nb_catch(5, 10), 
  pool_strat = pool_max_size(20)
)

## tests ----
test_that("variable_design constructor (target_size)", {
  expect_equal(class(var_target), c("variable_design", "sample_design"))
  expect_equal(var_target$catch_dist, nb_catch(5, 10))
  expect_equal(var_target$pool_strat, pool_target_number(20))
  expect_equal(var_target$sensitivity, 1)
  expect_equal(var_target$specificity, 1)
})

test_that("variable_design constructor (max_size)", {
  expect_equal(class(var_max), c("variable_design", "sample_design"))
  expect_equal(var_max$pool_strat, pool_max_size(20))
})

test_that("null variable_design", {
  expect_error(variable_design())
})

# Helpers ----
test_that("get_optimise_subclass", {
  expect_equal(get_optimise_subclass(NULL, NULL), "need_sN")
  expect_equal(get_optimise_subclass(1, NULL), "need_N")
  expect_equal(get_optimise_subclass(NULL, 1), "need_s")
})