#' S3 class for power and sample size calculation outputs with pooled tests
#'
#' As many variables are shared across groups of functions (e.g. `power_pool`
#' and `sample_size_pool`), this class and print method ensures that outputs are
#' stored and displayed consistently, summarising calculations results and the
#' inputs for context.
#' 
#' @param sensitivity .
#' @param specificity .
#' @param prev_null .
#' @param prev_alt .
#' @param correlation .
#' @param sig_level .
#' @param power .
#' @param alternative .
#' @param pool_size .
#' @param pool_number .
#' @param cluster_number . 
#' @param total_pools .
#' @param total_units .
#' @param text chr Explanatory summary text to be printed at the end
#'
#' @return An object of class \code{power_size_results} containing selected
#'   input parameters and results.
#' @export
#'
#' @examples
#' # For power_pool()
#' result <- power_size_results(
#'   sensitivity = 1, specificity = 1, prev_null = 0.01, prev_alt = 0.02,
#'   correlation = 0, sig_level = 0.05, power = 0.76, # rounded in e.g. only
#'   alternative = "greater", pool_size = 10, pool_number = 2, 
#'   cluster_number = 50, total_pools = 100, total_units = 1000,
#'   text = "... has a statistical power of 0.762"
#' )
#' 
#' print(result) # pretty print
#' result$stat_test$power
power_size_results <- function(sensitivity, specificity, prev_null, prev_alt, 
                               correlation, sig_level, power, alternative,
                               pool_size, pool_number, cluster_number, 
                               total_pools, total_units, text) {
  
  # Group parameters to different lists for printing
  diag_test <- list(
    title = "DIAGNOSTIC TEST",
    sensitivity = sensitivity,
    specificity = specificity
  )
  
  prev <- list(
    title = "PREVALENCE",
    prev_null = prev_null,
    prev_alt = prev_alt,
    correlation = correlation
  )
  
  stat_test <- list(
    title = "STATISTICAL TEST",
    sig_level = sig_level,
    power = power,
    alternative = alternative
  )
  
  sample_design <- list(
    title = "SAMPLE DESIGN",
    pool_size = pool_size,
    pool_number = pool_number,
    cluster_number = cluster_number,
    total_pools = total_pools,
    total_units = total_units
  )
  
  results <- structure(
    list(
      diag_test = diag_test,
      prev = prev,
      stat_test = stat_test,
      sample_design = sample_design,
      text = text 
    ),
    class = "power_size_results"
  )
  return(results)
}

#' @method print power_size_results
#' @export
print.power_size_results <- function(x, ...) {
  # Remove 'text' from instance for printing
  text_idx <- which(names(x) == "text")
  text <- x[[text_idx]]
  x <- x[-text_idx]
  
  for (list in x) {
    cat(paste("\n", list$title, "\n"))
    for (name in names(list[-1])) { # skip title
      if (name == "power") {
        value <- round(list[[name]], 3)
      } else {
        value = list[[name]]
      }
      cat(paste(format(name, width = 15L, justify = "right"), value, sep = " = "
), sep = "\n"
      )
    } 
  }
  cat("\n", text)
  invisible(x)
}