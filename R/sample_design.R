#' fixed_design s3 constructor
#'
#' Stores parameters related to the sampling design. Aims to reduce having to
#' input each param separately across functions (e.g. power/optimise).
#'
#' @param pool_size numeric/NULL The number of units per pool. Must be a numeric
#'   value greater than 0.
#' @param pool_number numeric/NULL The number of pools per cluster. Numeric
#'   inputs must be an integer greater than or equal to 1.
#' @param sensitivity numeric The probability that the test correctly identifies
#'   a true positive. Must be a numeric value between 0 and 1, inclusive of
#'   both. A value of 1 indicates that the test can perfectly identify all true
#'   positives.
#' @param specificity numeric The probability that the test correctly identifies
#'   a true negative. Must be a numeric value between 0 and 1, inclusive of
#'   both. A value of 1 indicates that the test can perfectly identify all true
#'   negatives.
#'
#' @return An object of class \code{sample_design}
#' @export
#'
#' @examples
#' perfect <- fixed_design(pool_size = 10)
#' 
#' imperfect <- fixed_design(
#'   pool_size = 10, pool_number = NULL, sensitivity = 0.95, specificity = 0.99
#' )
fixed_design <- function(pool_size = NULL,
                         pool_number = NULL,
                         sensitivity = 1,
                         specificity = 1) {

  if (!is.null(pool_size)) {
    check_geq2(pool_size, 0)
  }
  if (!is.null(pool_number)) {
    check_geq2(pool_number, 0)
  }
  # sens and spec cannot be NULL
  check_in_range2(sensitivity)
  check_in_range2(specificity)

  structure(
    list(
      pool_size = pool_size,
      pool_number = pool_number,
      sensitivity = sensitivity,
      specificity = specificity
    ),
    class = c("fixed_design", "sample_design")
  )
}

