#' S3 class for power and sample size calculation outputs with pooled tests
#'
#' As many variables are shared across groups of functions (e.g. `power_pool`
#' and `sample_size_pool`), this class and print method ensures that outputs are
#' stored and displayed consistently, summarising calculations results and the
#' inputs for context.
#'
#' @param ... Input and inferred variables from relevant functions.
#'
#' @return An object of class \code{power_size_results} containing selected
#'   input parameters and results.
#' @export
#'
#' @examples
#' # For power_pool()
#' result <- power_size_results(pool_size = pool_size, pool_number = pool_number, cluster_number = cluster_number,
#'                              prev_null = theta0, prev_alt = theta0, sig_level = sig_level,
#'                              power = power, alternative = alternative, target = "power")
#' print(result)

power_size_results <- function(pool_size, pool_number, cluster_number,
                               prev_null, prev_alt, sig_level, power,
                               alternative, target) {
  results <- structure(
    list(
      pool_size = pool_size,
      pool_number = pool_number,
      cluster_number = cluster_number,
      prev_null = prev_null,
      prev_alt = prev_alt,
      sig_level = sig_level,
      power = power,
      alternative = alternative,
      target = target
    ),
    class = "power_size_results"
  )
  return(results)
}

#' @method print power_size_results
#' @export
print.power_size_results <- function(x,...) {
  # Remove 'target' from instance to decorate line of the
  # inferred variable(s)
  target_idx <- which(names(x) == "target")
  target <- x[[target_idx]]
  x <- x[-target_idx]

  # decorate target line
  name_idx <- which(names(x) == target)
  names(x)[name_idx] <- paste("-->", names(x)[name_idx]) 
   
  # adapted from `stats::power.prop.test`
  cat(paste(format(names(x), width = 15L, justify = "right"),
            format(x), sep = " = "), sep = "\n")
  invisible(x)
}