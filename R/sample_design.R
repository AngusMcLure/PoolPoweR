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
#' @param total_units numeric/NULL internal use only for cases when this needs
#' to be Inf.
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
                         specificity = 1,
                         total_units = NULL) {

  ## Input checks ----
  # allow NULLs for optimise functions to identify which 
  # variable should be optimised
  if (!is.null(pool_size)) {
    check_geq2(pool_size, 0)
  }
  # Allow NA when optimise_prevalence(correlation = 0)
  if (!is.null(pool_number) && !is.na(pool_number)) {
    check_geq2(pool_number, 0)
  }
  if (!is.null(total_units) && !is.na(total_units)) {
    check_geq2(total_units, 0)
  }
  # sens and spec cannot be NULL
  check_in_range2(sensitivity)
  check_in_range2(specificity)

  ## Subclasses for optimisation ----
  # To dispatch optimise() based on NULLs 
  opt_class <- null_pools(pool_size, pool_number)
  opt_class <- paste0("fixed_design_optimise_", opt_class)

  ## Parse total parameters ----
  if (opt_class == "fixed_design_optimise_complete_params") {
    if (is.null(total_units)) {
      # When pool_size and pool_number are filled
      total_units = pool_size * pool_number
    }
    # Ensure that manually input total_units matches. 
    # TODO: Best to replace total_units arg with ...
    stopifnot(total_units == pool_size * pool_number || total_units == Inf || is.na(total_units))
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
                            pool_strat = NULL,
                            pool_strat_family = NULL,
                            sensitivity = 1,
                            specificity = 1) {

  # sens and spec cannot be NULL
  check_in_range2(sensitivity)
  check_in_range2(specificity)
  
  #must provide either pool_strat or pool_strat_family
  if(is.null(pool_strat) & is.null(pool_strat_family)){
    stop('Must provide either a valid `pool_strat` or `pool_strat_family`')
  }
  
  #
  if(!is.null(pool_strat) & !inherits(pool_strat, 'pool_strat')){
    stop('`pool_strat` must be a object of class `pool_strat`')
  }
  
  #Fill pool_strat_family based on pool_strat if former is not provided
  if(!is.null(pool_strat) & is.null(pool_strat_family)){
    pool_strat_family <- attr(pool_strat, 'family')
  }
  
  #Note there is currently no check that the pool_strat_family matches
  #pool_strat if both are supplied

  structure(
    list(
      catch_dist = catch_dist,
      pool_strat = pool_strat,
      pool_strat_family = pool_strat_family,
      sensitivity = sensitivity,
      specificity = specificity
    ),
    class = c("variable_design", "sample_design")
  )
}


null_pools <- function(pool_size, pool_number) {
  if (is.null(pool_size) && is.null(pool_number)) {
    return("sN")
  } else if (is.null(pool_size)) {
    return("s")
  } else if (is.null(pool_number)) {
    return("N")
  } else {
    return("complete_params")
  }
}

is_perfect_test <- function(x, ...) {
  UseMethod("is_perfect_test")
}

#' @method is_perfect_test sample_design 
#' @export
#' @noRd
is_perfect_test.sample_design <- function(x, ...) {
  if (x$sensitivity == 1 && x$specificity == 1) {
    return(TRUE)
  }
  return(FALSE)
}
