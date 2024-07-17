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
#' imperfect tests, and either fixed sample sizes 
#' (`design_effect(fixed_design, ...)`) or variable sample sizes 
#' (`design_effect(variable_design, ...)`).
#'
#' @param x a sample_design object 
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
#' @param form string The distribution used to model the cluster-level
#'   prevalence and correlation of units within cluster. Select one of "beta",
#'   "logitnorm" or "cloglognorm".
#'
#' @return A numeric value of the design effect `D`.
#' @export
#' 
#' @examples
#' design_effect(fixed_design(10, 2, cluster_number = 1), prevalence = 0.01, correlation = 0.05)
#' 
#' vd <- variable_design(nb_catch(10, 13), pool_target_number(20))
#' design_effect(vd, prevalence = 0.01, correlation = 0.05)
design_effect <- function(x, prevalence, correlation, form) {
  check_in_range2(prevalence)
  check_in_range2(correlation)
  # No input check for form as done in downstream functions/methods
  UseMethod("design_effect")
}

#' @rdname design_effect
#' @method design_effect fixed_design
#' @export
design_effect.fixed_design <- function(x,
                                       prevalence,
                                       correlation,
                                       form = "logitnorm") {
  
  x$pool_number * x$pool_size * fi_pool(pool_size = 1, prevalence, x$sensitivity, x$specificity) *
    solve(fi_pool_cluster(
      x$pool_size, x$pool_number, prevalence,
      correlation, x$sensitivity, x$specificity, form)
    )[1, 1]
}

#' @rdname design_effect
#' @method design_effect variable_design
#' @export
design_effect.variable_design <- function(x,
                                          prevalence,
                                          correlation,
                                          form = "logitnorm") {
  
  mean(x$catch_dist) * fi_pool(pool_size = 1, prevalence, x$sensitivity, x$specificity) *
    solve(fi_pool_cluster_random(
      x$catch_dist, x$pool_strat, prevalence,
      correlation, x$sensitivity, x$specificity, form)
    )[1, 1]
}
