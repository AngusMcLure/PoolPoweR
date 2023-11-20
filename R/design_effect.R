design_effect <- function(pool_size,
                          pool_number,
                          prevalence,
                          correlation,
                          sensitivity,
                          specificity,
                          form = "beta") {

  pool_number * pool_size * fi_pool(pool_size = 1, prevalence, sensitivity, specificity) *
    solve(fi_pool_cluster(
      pool_size, pool_number, prevalence,
      correlation, sensitivity, specificity, form)
    )[1, 1]
}
