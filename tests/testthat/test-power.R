test_that(
  "power_pool()", {
  act <- power_pool(pool_size = 10, pool_number = 2, cluster_number = 50, prevalence_null = 0.01, prevalence_alt = 0.02)
  expect_equal(act, 0.7617911, tolerance = 1e-6) 

  act <- power_pool(pool_size = 10, pool_number = 2, cluster_number = 50, prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01)
  expect_equal(act, 0.6213165, tolerance = 1e-6)
  }
)

test_that(
  "sample_size_pool()", {
    act <- sample_size_pool(pool_size = 10, pool_number = 2, prevalence_null = 0.01, prevalence_alt = 0.02)
    exp <- list(clusters = 55, pools = 110, units = 1100)
    expect_equal(act, exp)
    
    act <- sample_size_pool(pool_size = 10, pool_number = 2, prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01)
    exp <- list(clusters = 74, pools = 148, units = 1480)
    expect_equal(act, exp)
  }
)

test_that(
  "power_pool_random()", {
    act <- power_pool_random(nb_catch(20,25), pool_target_number(2), cluster_number = 50, prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01)
    expect_equal(act, 0.613303, tolerance = 1e-6)
  }
)

test_that(
  "sample_size_pool_random()", {
    act <- sample_size_pool_random(nb_catch(20,25), pool_target_number(2), prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01)
    exp <- list(clusters = 75, expected_pools = 150, expected_units = 1500)
    expect_equal(act, exp)
  }
)
