#' Calculate the design effect for pooled testing.
#'
#' This function calculates the design effect (D) for survey designs using pool
#' testing compared to a simple random survey with individual tests of the same
#' number of units. This allows the comparison of the Fisher Information per
#' unit sampled across different pooling and sampling strategies. A design
#' effect `D>1` (`D<1`) indicates that the pooling/sampling strategy reduces
#' (increases) the Fisher information per unit; the total sample size will have
#' to be multiplied by a factor of D to achieve the same degree of precision in
#' estimating prevalence as a simple random survey with individual tests.
#' Supports both cluster and simple random sampling with perfect or imperfect
#' tests.
#'
#' @param pool_size numeric The number of units per pool. Must be a numeric
#'   value greater than or equal to 0.
#' @param pool_number numeric The number of pools per cluster. Must be a numeric
#'   value greater than or equal to 0.
#' @param prevalence numeric The proportion of units that carry the marker of
#'   interest (i.e. true positive). Must be be a numeric value between 0 and 1,
#'   inclusive of both.
#' @param correlation numeric The correlation between test results within a
#'   single cluster (units in different clusters are assumed to be
#'   uncorrelated). Must be a numeric value between 0 and 1, inclusive of both.
#'   A value of 1 indicates that units within clusters are perfectly correlated
#'   (there are no differences units within a single cluster). A value of 0
#'   indicates that units within clusters are no more correlated than units in
#'   different clusters.
#' @param sensitivity numeric The probability that the test correctly identifies
#'   a true positive. Must be a numeric value between 0 and 1, inclusive of
#'   both. A value of 1 indicates that the test can perfectly identify all true
#'   positives.
#' @param specificity numeric The probability that the test correctly identifies
#'   a true negative. Must be a numeric value between 0 and 1, inclusive of
#'   both. A value of 1 indicates that the test can perfectly identify all true
#'   negatives.
#' @param form string The distribution used to model the cluster-level
#'   prevalence and correlation of units within cluster. Select one of "beta",
#'   "logitnorm" or "cloglognorm". See details.
#'
#' @return A numeric value of the design effect `D`.
#' @export
#'
#' @examples
#' design_effect(
#'   pool_size = 5, pool_number = 10, prevalence = 0.01,
#'   correlation = 0.05, sensitivity = 0.99, specificity = 0.95
#'   )
design_effect <- function(pool_size,
                          pool_number,
                          prevalence,
                          correlation,
                          sensitivity,
                          specificity,
                          form = "beta") {

  check_input("pool_size", pool_size)
  check_input("pool_number", pool_number)
  check_input("prevalence", prevalence)
  check_input("correlation", correlation)
  check_input("sensitivity", sensitivity)
  check_input("specificity", specificity)
  check_input("form", form)

  pool_number * pool_size * fi_pool(pool_size = 1, prevalence, sensitivity, specificity) *
    solve(fi_pool_cluster(
      pool_size, pool_number, prevalence,
      correlation, sensitivity, specificity, form)
    )[1, 1]
}
