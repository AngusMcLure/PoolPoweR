# fixtures ----
exp1 <- function(act) {
  # Common output across tests
  expect_equal(act$sample_design$cluster_number, 75)
  expect_equal(act$sample_design$exp_total_pools, 150)
  expect_equal(act$sample_design$exp_total_units, 1500)
}

fd <- fixed_design(pool_size = 10, pool_number = 2, cluster_number = 50)
fd_imp <- fixed_design(15, 3, 20, 0.99, 0.98)

# power_pool ----
test_that(
  "power_pool() no corr", {
    act <- pool_power(fd, prevalence_null = 0.01, prevalence_alt = 0.02)
    expect_equal(act$stat_test$power, 0.7617911, tolerance = 1e-6) 
    # tests for totals as they aren't direct inputs
    expect_equal(act$sample_design$total_pools, 100) 
    expect_equal(act$sample_design$total_units, 1000) 
    expect_equal(act$text, "A survey design using a perfect diagnostic test on pooled samples with the above parameters has a statistical power of 0.762") 
  }
)

test_that(
  "power_pool() with corr", {
    act <- pool_power(fd, prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01)
    expect_equal(act$stat_test$power, 0.6213165, tolerance = 1e-6)
    expect_equal(act$sample_design$total_pools, 100) 
    expect_equal(act$sample_design$total_units, 1000) 
    expect_equal(act$text, "A survey design using a perfect diagnostic test on pooled samples with the above parameters has a statistical power of 0.621") 
  }
)

test_that(
  "power_pool() imperfect", {
    act <- pool_power(fd_imp, prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01)
    expect_equal(act$stat_test$power, 0.4350865, tolerance = 1e-6)
    expect_equal(act$sample_design$total_pools, 60) 
    expect_equal(act$sample_design$total_units, 900) 
    expect_equal(act$text, "A survey design using an imperfect diagnostic test on pooled samples with the above parameters has a statistical power of 0.435") 
  }
)

test_that(
  "power_pool() links", {
  act <- pool_power(fd, prevalence_null = 0.01, prevalence_alt = 0.02, link = "identity")
  expect_equal(act$stat_test$power, 0.8448746, tolerance = 1e-6) 
  
  act <- pool_power(fd, prevalence_null = 0.01, prevalence_alt = 0.02, link = "log")
  expect_equal(act$stat_test$power, 0.7598709, tolerance = 1e-6) 
  
  act <- pool_power(fd, prevalence_null = 0.01, prevalence_alt = 0.02, link = "cloglog")
  expect_equal(act$stat_test$power, 0.7608382, tolerance = 1e-6) 
  }
)

# power_pool_random ----
test_that(
  "power_pool_random()", {
    act <- power_pool_random(nb_catch(20,25), pool_target_number(2), cluster_number = 50, prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01)
    expect_equal(act$stat_test$power, 0.613303, tolerance = 1e-6)
  }
)

test_that(
  "power_pool_random() with links", {
    act <- power_pool_random(nb_catch(20,25), pool_target_number(2), cluster_number = 50, prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01, link = "cloglog")
    expect_equal(act$stat_test$power, 0.6117088, tolerance = 1e-6)
    
    act <- power_pool_random(nb_catch(20,25), pool_target_number(2), cluster_number = 50, prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01, link = "log")
    expect_equal(act$stat_test$power, 0.6100939, tolerance = 1e-6)
    
    act <- power_pool_random(nb_catch(20,25), pool_target_number(2), cluster_number = 50, prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01, link = "identity")
    expect_equal(act$stat_test$power, 0.7568458, tolerance = 1e-6)
  }
)

test_that(
  "sample_size_pool()", {
    act <- sample_size_pool(pool_size = 10, pool_number = 2, prevalence_null = 0.01, prevalence_alt = 0.02)
    expect_equal(act$sample_design$cluster_number, 55)
    expect_equal(act$sample_design$total_pools, 110)
    expect_equal(act$sample_design$total_units, 1100)
    expect_equal(
      act$text, 
      "A survey design using a perfect diagnostic test on pooled samples with the above parameters requires a total of 55 clusters, 110 total pools, and 1100 total units."
    )
  }
)
    
test_that(
  "sample_size_pool with corr", {
    # Correlation
    act <- sample_size_pool(pool_size = 10, pool_number = 2, prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01)
    expect_equal(act$sample_design$cluster_number, 74)
    expect_equal(act$sample_design$total_pools, 148)
    expect_equal(act$sample_design$total_units, 1480)
    expect_equal(
      act$text,
      "A survey design using a perfect diagnostic test on pooled samples with the above parameters requires a total of 74 clusters, 148 total pools, and 1480 total units."
    )
  }
)

test_that(
  "sample_size_pool imperfect test", {
    act <- sample_size_pool(pool_size = 10, pool_number = 2, prevalence_null = 0.01, prevalence_alt = 0.02, sensitivity = 0.98, specificity = 0.99)
    expect_equal(act$sample_design$cluster_number, 61)
    expect_equal(act$sample_design$total_pools, 122)
    expect_equal(act$sample_design$total_units, 1220)
    expect_equal(
      act$text, 
      "A survey design using an imperfect diagnostic test on pooled samples with the above parameters requires a total of 61 clusters, 122 total pools, and 1220 total units."
    )
  }
)

test_that(
  "sample_size_pool() links", {
    act <- sample_size_pool(pool_size = 10, pool_number = 2, prevalence_null = 0.01, prevalence_alt = 0.02, link = "identity")
    exp <- list(clusters = 43, pools = 86, units = 860)
    expect_equal(act$sample_design$cluster_number, 43)
    expect_equal(act$sample_design$total_pools, 86)
    expect_equal(act$sample_design$total_units, 860)
    
    act <- sample_size_pool(pool_size = 10, pool_number = 2, prevalence_null = 0.01, prevalence_alt = 0.02, link = "log")
    exp <- list(clusters = 55, pools = 110, units = 1100) # same as link = "logit"
    expect_equal(act$sample_design$cluster_number, 55)
    expect_equal(act$sample_design$total_pools, 110)
    expect_equal(act$sample_design$total_units, 1100)
    
    act <- sample_size_pool(pool_size = 10, pool_number = 2, prevalence_null = 0.01, prevalence_alt = 0.02, link = "cloglog")
    # same as link = "logit"
    expect_equal(act$sample_design$cluster_number, 55)
    expect_equal(act$sample_design$total_pools, 110)
    expect_equal(act$sample_design$total_units, 1100)
  }
)

test_that(
  "sample_size_pool_random()", {
    act <- sample_size_pool_random(nb_catch(20,25), pool_target_number(2), prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01)
    exp1(act)
  }
)

test_that(
  "sample_size_pool_random() links", {
    act <- sample_size_pool_random(nb_catch(20,25), pool_target_number(2), prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01, link = "cloglog")
    exp1(act)
    
    act <- sample_size_pool_random(nb_catch(20,25), pool_target_number(2), prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01, link = "log")
    exp1(act)
    
    act <- sample_size_pool_random(nb_catch(20,25), pool_target_number(2), prevalence_null = 0.01, prevalence_alt = 0.02, correlation = 0.01, link = "identity")
    expect_equal(act$sample_design$cluster_number, 59)
    expect_equal(act$sample_design$exp_total_pools, 118)
    expect_equal(act$sample_design$exp_total_units, 1180)
  }
)