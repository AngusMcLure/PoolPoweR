#' S3 class for power and sample size calculation outputs with pooled tests
#'
#' As many variables are shared across groups of functions (e.g. `power_pool`
#' and `sample_size_pool`), this class and print method ensures that outputs are
#' stored and displayed consistently, summarising calculations results and the
#' inputs for context.
#' 
#' @param design required; sample design
#' @param prev_null,prev_alt,prev either prev, or prev_null and prev_alt required
#' @param correlation required
#' @param sig_level 
#' @param power required
#' @param alternative 
#' @param cluster_number required 
#'
#' @return An object of class \code{power_size_results} containing selected
#'   input parameters and results.
#' @export
#'


power_size_results <- function(design,
                               cluster_number,
                               prev_null = NULL,
                               prev_alt = NULL,
                               prev = NULL,
                               correlation,
                               sig_level = NULL,
                               power,
                               alternative = NULL) {
  
  # Group parameters to different lists for printing
  diag_test <- list(
    title = "DIAGNOSTIC TEST",
    sensitivity = design$sensitivity,
    specificity = design$specificity
  )
  
  prev <- list(
    title = "POPULATION",
    `prevalence (null)` = prev_null,
    `prevalence (alternative)` = prev_alt,
    `prevalence` = prev,
    correlation = correlation
  )
  
  stat <- list(
    title = "STATISTICAL PROPERTIES",
    `significance level` = sig_level,
    power = power,
    alternative = alternative
  )
  
  #
  if(inherits(design,'fixed_design')){
    
    total_pools <- cluster_number * design$pool_number
    total_units <- cluster_number * design$pool_number * design$pool_size
    
    survey_design <- list(
      title = "SURVEY DESIGN",
      design = design,
      `pool size` = design$pool_size,
      `pools per cluster` = design$pool_number,
      `total pools` = total_pools,
      `total units` = total_units,
      `clusters` = cluster_number
    )
    
    text <- paste0(
      "A survey design using ", is_perfect_test_temp(design$sensitivity, design$specificity), 
      " diagnostic test on pooled samples with the above parameters requires a total of ",
      cluster_number, " clusters, ", 
      total_pools, " total pools, and ", 
      total_units, " total units."
    )
  }else if(inherits(design,'variable_design')){
    
    exp_total_units <- round(distrEx::E(design$catch_dist) * cluster_number, 1)
    exp_total_pools <- round(ev(\(catch) sum(design$pool_strat(catch)$pool_number),
                                design$catch_dist) * cluster_number, 1)
    
    # Prepare output
    
    survey_design <- list(
      title = "SURVEY DESIGN",
      design = design,
      `expected catch per site` = design$catch_mean,
      `pooling strategy` = design$pool_strat,
      `clusters` = cluster_number,
      `total expected pools` = exp_total_pools,
      `total expected units` = exp_total_units
    )
    
    text <- paste0(
      "A survey design using ", is_perfect_test_temp(design$sensitivity, design$specificity), 
      " diagnostic test on pooled samples with the above parameters requires a total of ",
      cluster_number, " clusters, ", 
      exp_total_pools, " expected total pools, and ", 
      exp_total_units, " expected total units."
    )
    
  }
  
  results <- structure(
    list(
      diag_test = diag_test,
      prev = prev,
      stat = stat,
      survey_design = survey_design,
      text = text 
    ),
    class = "power_size_results"
  )
  return(results)
}


#' @method print power_size_results
#' @export
print.power_size_results <- function(x, ...) {
  
  text <- x$text
  design <- x$survey_design$design
  
  # Remove 'text' and 'design' from instance for printing
  x <- x[names(x) != 'text']
  xsd <- x[['survey_design']]
  x[['survey_design']] <- xsd[names(xsd) != 'design']
  
  for (list in x) {
    cat(paste("\n", list$title, "\n"))
    for (name in names(list[-1])) { # skip title
      if (name == "power") {
        value <- round(list[[name]], 3)
      } else {
        value = list[[name]]
      }
      if(!is.null(value)){
        cat(paste(format(name, width = 25L, justify = "right"),
                  value, sep = " = "),
            sep = "\n")
      }
    } 
  }
  #cat("\n", text)
  invisible(x)
}

g_switch <- function(link) {
  switch(link,
         logit = stats::qlogis,
         cloglog = cloglog,
         log = log,
         identity = function(x){x}
  )
}

gdivinv_switch <- function(link) {
  switch(link,
         logit = function(x){x * (1-x)},
         cloglog = function(x){-log1p(-x) * (1-x)},
         log = function(x){x},
         identity = function(x){1})
}

is_perfect_test_temp <- function(sensitivity, specificity) {
  # TODO: replace usage with is_perfect_test(sample_design) method
  if (sensitivity == 1 && specificity == 1) {
    return("a perfect")
  }
  return("an imperfect")
}

