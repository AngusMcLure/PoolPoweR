#' S3 sample_design constructors
#'
#' Stores parameters related to the sampling design. Aims to reduce having to
#' input each param separately across functions (e.g. power/optimise). Can
#' either be of class `fixed_design` or `variable_design`.
#'
#' @param pool_size numeric/NULL The number of units per pool. Must be a numeric
#'   value greater than 0. `fixed_design` only.
#' @param pool_number numeric/NULL The number of pools per cluster. Numeric
#'   inputs must be an integer greater than or equal to 1. `fixed_design` only.
#' @param catch_dist An object of class `distribution` (e.g. produced by
#'   `nb_catch()`) defining the distribution of the possible catch. If
#'   `correlation = 0` the catch is for the whole survey. For `correlation > 0`
#'   the catch is per cluster (i.e. cluster size). `variable_design` only.
#' @param pool_strat function Defines a rule for how a number of units will be
#'   divided into pools. Must take a single numeric argument and return a named
#'   list of pool sizes and pool numbers. `pool_max_size()` and
#'   `pool_target_number` provide convenience functions for defining common
#'   pooling strategies. `variable_design` only.
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
#' fd_perfect <- fixed_design(pool_size = 10)
#' 
#' fd_imperfect <- fixed_design(
#'   pool_size = 10, pool_number = NULL, sensitivity = 0.95, specificity = 0.99
#' )
#' 
#' vd_target <- variable_design(
#'   catch_dist = nb_catch(10, 11),
#'   pool_strat = pool_target_number(20)
#' )
#' 
#' vd_max <- variable_design(
#'   catch_dist = nb_catch(10, 11),
#'   pool_strat = pool_max_size(20)
#' )
#' 
#' vd_max_imperfect <- variable_design(
#'   catch_dist = nb_catch(10, 11),
#'   pool_strat = pool_max_size(20),
#'   sensitivity = 0.95,
#'   specificity = 0.98
#' )
fixed_design <- function(pool_size = NULL,
                         pool_number = NULL,
                         sensitivity = 1,
                         specificity = 1) {

  ## Input checks ----
  # allow NULLs for optimise functions to identify which 
  # variable should be optimised
  if (!is.null(pool_size)) {
    check_geq2(pool_size, 0)
  }
  if (!is.null(pool_number)) {
    check_geq2(pool_number, 0)
  }
  # sens and spec cannot be NULL
  check_in_range2(sensitivity)
  check_in_range2(specificity)

  ## Subclasses for optimisation ----
  # To dispatch optimise() based on NULLs 
  opt_class <- get_optimise_subclass(pool_size, pool_number)
  
  ## Parse total parameters ----
  # TODO: Add total_pools here once number of clusters is added
  if (opt_class == "complete_params") {
    total_units = pool_size * pool_number
  } else {
    total_units = NA
  }
  
  ## Output ----
  structure(
    list(
      pool_size = pool_size,
      pool_number = pool_number,
      total_units = total_units,
      sensitivity = sensitivity,
      specificity = specificity
    ),
    class = c(opt_class, "fixed_design", "sample_design")
  )
}

#' @rdname fixed_design
#' @export
variable_design <- function(catch_dist, 
                            pool_strat,
                            sensitivity = 1,
                            specificity = 1) {

  # sens and spec cannot be NULL
  check_in_range2(sensitivity)
  check_in_range2(specificity)

  structure(
    list(
      catch_dist = catch_dist,
      pool_strat = pool_strat,
      sensitivity = sensitivity,
      specificity = specificity
    ),
    class = c("variable_design", "sample_design")
  )
}


get_optimise_subclass <- function(pool_size, pool_number) {
  if (is.null(pool_size) && is.null(pool_number)) {
    return("need_sN")
  } else if (is.null(pool_size)) {
    return("need_s")
  } else if (is.null(pool_number)) {
    return("need_N")
  } else {
    return("complete_params")
  }
}
