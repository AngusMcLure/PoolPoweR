#' Calculate the design effect for pooled testing.
#'
#' These functions calculate the design effect (D) for survey designs using pool
#' testing compared to a simple random survey with individual tests of the same
#' number of units. This allows the comparison of the Fisher Information per
#' unit sampled across different pooling and sampling strategies. A design
#' effect `D>1` (`D<1`) indicates that the pooling/sampling strategy reduces
#' (increases) the Fisher information per unit; the total sample size will have
#' to be multiplied by a factor of D to achieve the same degree of precision in
#' estimating prevalence as a simple random survey with individual tests. The
#' functions support cluster and simple random sampling with perfect or
#' imperfect tests, and either fixed sample sizes (`design_effect()`) or random
#' sample sizes (`design_effect_random()`).
#'
#' @param pool_size numeric The number of units per pool. Must be a numeric
#'   value greater than or equal to 0.
#' @param pool_number numeric The number of pools per cluster. Must be a numeric
#'   value greater than or equal to 0.
#' @param catch_dist An object of class `distribution` (e.g. produced by
#'   `nb_catch()`) defining the distribution of the possible catch. If
#'   `correlation = 0` the catch is for the whole survey. For `correlation > 0`
#'   the catch is per cluster (i.e. cluster size).
#' @param pool_strat function Defines a rule for how a number of units will be
#'   divided into pools. Must take a single numeric argument and return a named
#'   list of pool sizes and pool numbers. `pool_max_size()` and
#'   `pool_target_number` provide convenience functions for defining common
#'   pooling strategies.
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
                          form = "beta"){

  pool_number * pool_size * fi_pool(pool_size = 1, prevalence, sensitivity, specificity) *
    solve(fi_pool_cluster(
      pool_size, pool_number, prevalence,
      correlation, sensitivity, specificity, form)
    )[1, 1]
}


#' @rdname design_effect
#' @export
design_effect_random <- function(catch_dist,
                                 pool_strat,
                                 prevalence,
                                 correlation,
                                 sensitivity,
                                 specificity,
                                 form = "beta") {
  
  mean(catch_dist) * fi_pool(pool_size = 1, prevalence, sensitivity, specificity) *
    solve(fi_pool_cluster_random(
      catch_dist, pool_strat, prevalence,
      correlation, sensitivity, specificity, form)
    )[1, 1]
}
