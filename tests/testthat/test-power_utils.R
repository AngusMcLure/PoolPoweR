test_that(
  "power_size_results constructor", {
  # e.g. from power_pool() example
  act <- power_size_results(
    sensitivity = 1, 
    specificity = 1, 
    prev_null = 0.01,
    prev_alt = 0.02,
    correlation = 0,
    sig_level = 0.05, 
    power = 0.76, 
    alternative = "greater",
    pool_size = 10, 
    pool_number = 2, 
    cluster_number = 50,
    total_pools = 100,
    total_units = 1000,
    text = "A survey design using a perfect diagnostic test on pooled samples with the above parameters has a statistical power of 0.762"
  )
  expect_s3_class(act, "power_size_results")
  expect_equal(length(act), 5)
  }
)