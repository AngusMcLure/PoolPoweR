optimised_results <- function(sample_design,
                              prevalence,
                              cost_unit,
                              cost_pool,
                              cost_period = NA,
                              cost_cluster,
                              correlation) {
  
  msg <- optimise_msg(
    sample_design$pool_size, 
    sample_design$pool_number,
    correlation
  )
  
  structure(
    list(
      sample_design = sample_design,
      prevalence = prevalence,
      cost_unit = cost_unit,
      cost_pool = cost_pool,
      cost_period = NA,
      cost_cluster = cost_cluster,
      correlation = correlation,
      message = msg
    ), class = c("optimised_results")
    # TODO: Will likely need to add classes for fixed/variable_design
    # maybe even clustered/unclustered based on the correlation value to help
    # with triaging functions in PoolTools (analysis_type()).
  )
  
}


pluralise <- function(obj, text) {
  ifelse(obj > 1, paste0(text, "s"), text)
}

optimise_msg <- function(pool_size, pool_number, correlation) {
  # Generate prose to explain the outcome of the optimisation results.
  # Aim to make this the exact output of what you want to see on PoolTools to
  # minimise processing on the Shiny side.
  pre <- "The optimal design is to sample "
  p_units <- pluralise(pool_size, "unit")
  if (!is.na(correlation) && !is.infinite(pool_number)) {
    paste0(
      pre, pool_size * pool_number, " units per collection site, across ",
      pool_number, " pools with ", pool_size, " ", p_units, " each pool."
    )
  } else {
    paste0(pre, pool_size, " ", p_units, " per pool." )
  }
}