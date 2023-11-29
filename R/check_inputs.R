check_input <- function(argument_name, value) {
  na <- c('cost_cluster', 'correlation')
  forms <- c('beta', 'logitnorm', 'cloglognorm')
}

check_geq <- function(argument_name, input_value) {
  # Check that an input value is numeric, and greater than `min`
  # Set `min` dynamically depending on the argument requirements
  geq0 <- c('pool_size', 'pool_number', 'cost_unit',
            'cost_pool', 'cost_cluster', 'interval')
  geq1 <- c('max_s', 'max_N')
  min <- 0
  if(argument_name %in% geq1) min <- 1

  if(!is.numeric(input_value) | input_value < min) {
    cli::cli_alert_info("{.field {argument_name}} needs to be a numeric value {min} or greater.")
  }
  if(!is.numeric(input_value)) {
    cli::cli_alert_danger("{.val {input_value}} is a {class(input_value)}.")
  }
  if(input_value < min) {
    cli::cli_alert_danger("{.val {input_value}} is < {min}")
  }
}
