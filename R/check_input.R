# NOTE: Will discontinue - use check_input2.R 

check_input <- function(argument_name, input_value) {
  # Wrapper function to triage arguments, so one function can be used for all
  # inputs instead of remembering which one to use.
  geq0 <- c('pool_size', 'pool_number', 'cost_unit', 'cost_pool', 'interval')
  geq1 <- c('max_s', 'max_N')
  range <- c('prevalence', 'sensitivity', 'specificity')

  if(argument_name %in% geq0) check_geq(argument_name, input_value, min = 0)
  else if(argument_name %in% geq1) check_geq(argument_name, input_value, min = 1)
  else if(argument_name %in% range) check_in_range(argument_name, input_value)
  else if(argument_name == "correlation") check_rho(input_value)
  else if(argument_name == "form") check_form(input_value)
  else if(argument_name == "real_scale") check_scale(input_value)
  else if(argument_name == "cost_cluster") check_cost_cluster(input_value)
}

check_geq <- function(argument_name, input_value, min) {
  # Check that an input value is numeric, and greater than `min`
  # `min` either 0 or 1
  # Stop condition for developers only
  geq0 <- c('pool_size', 'pool_number', 'cost_unit',
            'cost_pool', 'cost_cluster', 'interval')
  geq1 <- c('max_s', 'max_N')
  if(!argument_name %in% c(geq0, geq1)) stop("Needs to be one of the accepted_args")

  if(!is.numeric(input_value) | input_value < min) {
    message(glue::glue("{argument_name} must be a numeric value {min} or greater."))
  }
    if(!is.numeric(input_value)) {
    stop(glue::glue("{input_value} is a {class(input_value)}."))
  }
  if(input_value < min) {
    stop(glue::glue("{input_value} is < {min}"))
  }
}

check_in_range <- function(argument_name, input_value) {
  # Is the input between 0 to 1, inclusive?
  # Stop condition for developers only
  accepted_args <- c('prevalence', 'correlation', 'sensitivity', 'specificity')
  if(!argument_name %in% accepted_args) stop("Needs to be one of the accepted_args")

  if(!is.numeric(input_value)) {
    stop(glue::glue("{input_value} is a {class(input_value)}."))
  }
  if(input_value < 0 | input_value > 1) {
    message(glue::glue("{argument_name} must be a numeric value between 0 and 1, inclusive."))
  }
  if(input_value < 0) stop(glue::glue("{input_value} is < 0"))
  if(input_value > 1) stop(glue::glue("{input_value} is > 1"))
}

check_rho <- function(rho) {
  alert_msg <- "correlation must be a numeric value between 0 and 1, or NA"
  if(!is.numeric(rho) & !is.na(rho)) {
    message(alert_msg)
    stop(glue::glue("{rho} is a {class(rho)}."))
  }
  if(is.numeric(rho)) {
    # `<` and `>` can't deal with NAs
    if(rho < 0 || rho > 1) {
      message(alert_msg)
      if(rho < 0) stop(glue::glue("{rho} is < 0"))
      if(rho > 1) stop(glue::glue("{rho} is > 1"))
    }
  }
}

check_form <- function(form) {
  forms <- c('beta', 'logitnorm', 'cloglognorm', 'discrete')
  if(!form %in% forms) {
    stop(glue::glue("form must be one of 'beta', 'logitnorm', 'cloglognorm', or 'discrete'."))
  }
}

check_scale <- function(real_scale) {
  if(!is.logical(real_scale)) {
    message("real_scale must be either TRUE/FALSE")
    stop(glue::glue("{real_scale} is not TRUE/FALSE"))
  }
}

check_cost_cluster <- function(cost_cluster) {
  # Function really similar to check_rho()
  alert_msg <- "cost_cluster must be a numeric value 0 or greater, or NA"
  if(!is.numeric(cost_cluster) & !is.na(cost_cluster)) {
    message(alert_msg)
    stop(glue::glue("{cost_cluster} is a {class(cost_cluster)}."))
  }
  if(is.numeric(cost_cluster)) {
    # `<` and `>` can't deal with NAs
    if(cost_cluster < 0) {
      message(alert_msg)
      stop(glue::glue("{cost_cluster} is < 0"))
    }
  }
}
