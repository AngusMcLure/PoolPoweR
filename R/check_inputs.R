check_input <- function(argument_name, value) {
  na <- c('cost_cluster', 'correlation')
  forms <- c('beta', 'logitnorm', 'cloglognorm')
}

check_geq <- function(argument_name, input_value) {
  # Check that an input value is numeric, and greater than `min`
  # Stop condition for developers only
  geq0 <- c('pool_size', 'pool_number', 'cost_unit',
            'cost_pool', 'cost_cluster', 'interval')
  geq1 <- c('max_s', 'max_N')
  if(!argument_name %in% c(geq0, geq1)) stop("Needs to be one of the accepted_args")
  # So one function can be used for all >= checks
  min <- 0
  if(argument_name %in% geq1) min <- 1

  if(!is.numeric(input_value) | input_value < min) {
    cli::cli_alert_info("{.field {argument_name}} must be a numeric value {min} or greater.")
  }
  if(!is.numeric(input_value)) {
    cli::cli_alert_danger("{.val {input_value}} is a {class(input_value)}.")
  }
  if(input_value < min) {
    cli::cli_alert_danger("{.val {input_value}} is < {min}")
  }
}

check_in_range <- function(argument_name, input_value) {
  # Is the input between 0 to 1, inclusive?
  # Stop condition for developers only
  accepted_args <- c('prevalence', 'correlation', 'sensitivity', 'specificity')
  if(!argument_name %in% accepted_args) stop("Needs to be one of the accepted_args")

  if(input_value < 0 | input_value > 1) {
    cli::cli_alert_info("{.field {argument_name}} must be a numeric value between 0 and 1, inclusive.")
  }
  if(!is.numeric(input_value)) {
    cli::cli_alert_danger("{.val {input_value}} is a {class(input_value)}.")
  }
  if(input_value < 0) cli::cli_alert_danger("{.val {input_value}} is < 0")
  if(input_value > 1) cli::cli_alert_danger("{.val {input_value}} is > 1")
}
