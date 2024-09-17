#' Power and sample size calculations for detecting presence and estimating
#' population prevalence from pooled samples
#'
#' `power_threshold()` calculates the statistical power of a pooled survey
#' design to determine whether population prevalence is different from a
#' threshold.
#'
#' `power_detect()` calculations the statistical power of a pooled survey design
#' to detect a marker in a population, given it is present at a specified
#' prevalence.
#'
#' `sample_size_` family of functions calculate the sample size required for a
#' pooled survey to achieve the specified power for the given objective.
#'
#' @param design sample_design object.
#' @param cluster_number numeric The total number of clusters in a cluster
#'   survey design (should be greater than 1) or 1 surveys where all collection
#'   happens at a single site or collected via simple random sampling from the
#'   target population.
#' @param prevalence_null,prevalence_alt,prevalence numeric The proportion of
#'   units that carry the marker of interest (i.e. true positive). For
#'   `_threshold` functions `prevalence_null` is the threshold to compare to and
#'   `prevalence_alt` is the design prevalence. For `_detect` functions,
#'   `prevalence` is the design prevalence. Must be be a numeric value between 0
#'   and 1, inclusive of both.
#' @param correlation numeric The correlation between test results within a
#'   single cluster (units in different clusters are assumed to be
#'   uncorrelated). Must be a numeric value between 0 and 1, inclusive of both.
#'   A value of 1 indicates that units within clusters are perfectly correlated
#'   (there are no differences units within a single cluster). A value of 0
#'   indicates that units within clusters are no more correlated than units in
#'   different clusters.
#' @param power numeric The desired statistical power of the survey.
#' @param sig_level numeric Significance level for statistical test. Defaults to
#'   0.05. Must be strictly between 0 and 1.
#' @param alternative string The kind of comparison to make. If "greater" or
#'   "less" (default), computes power of tests for one-sided comparisons that
#'   population prevalence (`prevalence_alt`) is greater than or less than the
#'   reference on threshold prevalence (`prevalence_null`) respectively. If
#'   "two.sided" computes power for a two-sided comparison (`power_threshold()`
#'   only).
#' @param form string The distribution used to model the cluster-level
#'   prevalence and correlation of units within cluster. Select one of "beta",
#'   "logitnorm" or "cloglognorm". See details.
#' @param link string Transformation to be applied to prevalence estimates for
#'   the purposes of calculating confidence intervals. Options are 'identity'
#'   (i.e. no transformation), 'logit' (default), 'cloglog' and 'log'.
#' @param max_iter numeric Maximum number of iterations (possible catch sizes)
#'   to consider when calculating expected FI over random catch sizes. Generally
#'   needs to be large enough so that the nearly all catch sizes will be less
#'   than `max_iter` otherwise algorithm will terminate early (with a warning)
#' @param rel_tol numeric Relative tolerance for determining convergence when
#'   calculating expected FI over random catch sizes. Must be positive and
#'   should be much smaller than 1.
#' @param ... Additional parameters.
#'
#' @return The statistical power of the proposed design with regards to
#'   comparing prevalence to a threshold (`power_threshold()`) or a list with
#'   the sample size (number of clusters, pools, and units) required to achieve
#'   desired power (`sample_size_threshold()`)
#' @export
#'
#' @examples
#' # Fixed design examples ----
#' fd <- fixed_design(pool_size = 10, pool_number = 2)
#'
#' ## Unclustered ----
#' power_threshold(fd, cluster_number = 50,
#'            prevalence_null = 0.02, prevalence_alt = 0.1)
#'
#' sample_size_threshold(fd, prevalence_null = 0.02, prevalence_alt = 0.1)
#'
#' ## Clustered ----
#' power_threshold(fd, cluster_number = 50,
#'                 prevalence_null = 0.02, prevalence_alt = 0.01,
#'                 correlation = 0.1)
#'
#' sample_size_threshold(fd,
#'                       prevalence_null = 0.02, prevalence_alt = 0.01,
#'                       correlation = 0.1)
#'
#' # Variable design examples ----
#'
#' vd <- variable_design(nb_catch(20,25), pool_target_number(2))
#' power_threshold(vd, cluster_number = 50,
#'                 prevalence_null = 0.02, prevalence_alt = 0.01,
#'                 correlation = 0.1)
#'
#' sample_size_threshold(vd,
#'                       prevalence_null = 0.02, prevalence_alt = 0.01,
#'                       correlation = 0.1)

# _threshold ----

power_threshold <- function(design, 
                            cluster_number,
                            prevalence_null, 
                            prevalence_alt,
                            correlation,
                            sig_level,
                            alternative,
                            form,
                            link,
                            ...) {
  UseMethod("power_threshold")
}

#' @method power_threshold fixed_design
#' @export
power_threshold.fixed_design <- function(design,
                                         cluster_number,
                                         prevalence_null, 
                                         prevalence_alt,
                                         correlation = 0, 
                                         sig_level = 0.05, 
                                         alternative = 'less',
                                         form = 'logitnorm', 
                                         link = 'logit',
                                         ...) {
  # Input checks
  if (correlation > 0 & cluster_number <= 1) {
    stop('The number of clusters (cluster_number) must be (substantially) greater than 1 if there is non-zero correlation between units in a cluster')
  }
  
  if (correlation > 0 & cluster_number <= 10) {
    warning('Estimated power may be unreliable if number of clusters (cluster_number) is less than 10')
  }
  
  if (correlation == 0 & cluster_number * design$pool_number <= 10) {
    warning('Estimated power may be unreliable if total number of pools (cluster_number * pool_number) is less than 10')
  }
  
  thetaa <- prevalence_alt
  theta0 <- prevalence_null
  
  # Get link functions
  g <- g_switch(link)
  gdivinv <- gdivinv_switch(link)
  
  # Calculate Fisher information
  fia <- cluster_number * gdivinv(thetaa)^2 /
    solve(fi_pool_cluster(pool_size = design$pool_size,
                          pool_number = design$pool_number,
                          prevalence = thetaa,
                          correlation = correlation,
                          sensitivity = design$sensitivity,
                          specificity = design$specificity,
                          form = form))[1,1]
  
  fi0 <- cluster_number * gdivinv(theta0)^2/
    solve(fi_pool_cluster(pool_size = design$pool_size,
                          pool_number = design$pool_number,
                          prevalence = theta0,  #should this be theta0 or thetaa?
                          correlation = correlation,
                          sensitivity = design$sensitivity,
                          specificity = design$specificity,
                          form = form))[1,1]
  
  power <- switch(alternative,
                  less = stats::pnorm(((g(theta0) - g(thetaa))  - stats::qnorm(1-sig_level)/sqrt(fi0)) * sqrt(fia)),
                  greater = stats::pnorm(((g(thetaa) - g(theta0))  - stats::qnorm(1-sig_level)/sqrt(fi0)) * sqrt(fia)),
                  two.sided = stats::pnorm(((g(theta0) - g(thetaa))  - stats::qnorm(1-sig_level/2)/sqrt(fi0)) * sqrt(fia)) +
                    stats::pnorm(((g(thetaa) - g(theta0))  - stats::qnorm(1-sig_level/2)/sqrt(fi0)) * sqrt(fia)),
                  stop('invalid alternative. options are less, greater, and two.sided')
  )
  
  # Prepare output
  power_size_results(
    # sample design
    design = design,
    cluster_number = cluster_number,
    # prevalence
    prev_null = theta0,
    prev_alt = thetaa,
    correlation = correlation,
    # statistical test
    sig_level = sig_level,
    power = power,
    alternative = alternative
  )
  
}

#' @method power_threshold variable_design
#' @export
power_threshold.variable_design <- function(design,
                                            cluster_number,
                                            prevalence_null,
                                            prevalence_alt,
                                            correlation = 0,
                                            sig_level = 0.05,
                                            alternative = 'less',
                                            form = 'logitnorm',
                                            link = 'logit',
                                            max_iter = 10000,
                                            rel_tol = 1e-6,
                                            ...){

  
  if(!inherits(design$pool_strat, 'pool_strat')){
    stop('pool design must include a valid pooling strategy of class `pool_strat`',
    ' (not just a pooling strategy family)')}
  
  exp_pools <- ev(\(catch) sum(design$pool_strat(catch)$pool_number),
                  design$catch_dist) * cluster_number
  
  if(correlation > 0 & cluster_number <= 1){
    stop('The number of clusters (cluster_number) must be (substantially) greater than 1 if there is non-zero correlation between units in a cluster')
  }else if(cluster_number < 1){
    stop('The number of clusters (cluster_number) must be at least 1.')
  }
  
  if(correlation > 0 & cluster_number <= 10){
    warning('Estimated power may be unreliable if number of clusters (cluster_number) is less than 10')
  }
  
  if(correlation == 0 & exp_pools <= 10){
    warning('The expected number of pools from this survey is less than ', ceiling(exp_pools),
            '. Estimated power may be unreliable if total number of pools is less than 10')
  }
  
  thetaa <- prevalence_alt
  theta0 <- prevalence_null
  
  # Get link functions
  g <- g_switch(link)
  gdivinv <- gdivinv_switch(link)
  
  fia <- cluster_number * gdivinv(thetaa)^2 /
    solve(fi_pool_cluster_random(catch_dist = design$catch_dist, 
                                 pool_strat = design$pool_strat, 
                                 prevalence = thetaa,
                                 correlation = correlation,
                                 sensitivity = design$sensitivity,
                                 specificity = design$specificity,
                                 form = form,
                                 max_iter = max_iter,
                                 rel_tol = rel_tol-6))[1,1]
  
  fi0 <- cluster_number * gdivinv(theta0)^2 /
    solve(fi_pool_cluster_random(catch_dist = design$catch_dist, 
                                 pool_strat = design$pool_strat, 
                                 prevalence = theta0,  #should this be theta0 or thetaa?
                                 correlation = correlation,
                                 sensitivity = design$sensitivity,
                                 specificity = design$specificity,
                                 form = form,
                                 max_iter = max_iter,
                                 rel_tol = rel_tol-6))[1,1]
  
  power <- switch(alternative,
                  less    = stats::pnorm(((g(theta0) - g(thetaa))  - stats::qnorm(1-sig_level)/sqrt(fi0)) * sqrt(fia)),
                  greater = stats::pnorm(((g(thetaa) - g(theta0))  - stats::qnorm(1-sig_level)/sqrt(fi0)) * sqrt(fia)),
                  two.sided = stats::pnorm(((g(theta0) - g(thetaa))  - stats::qnorm(1-sig_level/2)/sqrt(fi0)) * sqrt(fia)) +
                    stats::pnorm(((g(thetaa) - g(theta0))  - stats::qnorm(1-sig_level/2)/sqrt(fi0)) * sqrt(fia)),
                  stop('invalid alternative. options are less, greater, and two.sided')
  )

  
  power_size_results( 
    # sample design
    design = design,
    cluster_number = cluster_number,
    # prevalence
    prev_null = theta0,
    prev_alt = thetaa,
    correlation = correlation,
    # statistical test
    sig_level = sig_level,
    power = power,
    alternative = alternative
  )
}

#' @rdname power_threshold
#' @export
sample_size_threshold <- function(design, 
                                  prevalence_null, 
                                  prevalence_alt,
                                  correlation,
                                  power, 
                                  sig_level,
                                  alternative,
                                  form,
                                  link,
                                  ...) {
  UseMethod("sample_size_threshold")
}

#' @method sample_size_threshold fixed_design 
#' @export
sample_size_threshold.fixed_design <- function(design,
                                               prevalence_null,
                                               prevalence_alt,
                                               correlation = 0,
                                               power = 0.8,
                                               sig_level = 0.05,
                                               alternative = "less",
                                               form = "logitnorm",
                                               link = "logit",
                                               ...) {
  thetaa <- prevalence_alt
  theta0 <- prevalence_null
  
  if (!(alternative %in% c("less", "greater"))) {
    stop("currently only supports one-sided tests. Valid options for alternative are less and greater")
  }
  
  if (alternative == "less" && theta0 < thetaa) {
    stop("If alternative == 'less', then prevalence.altnerative must be less than or equal to prevalence_null")
  }
  
  if (alternative == "greater" && theta0 > thetaa) {
    stop("If alternative == 'greater', then prevalence.altnerative must be greater than or equal to prevalence_null")
  }
  
  # Get link functions
  g <- g_switch(link)
  gdivinv <- gdivinv_switch(link)
  
  fia <- gdivinv(thetaa)^2 /
    solve(fi_pool_cluster(pool_size = design$pool_size,
                          pool_number = design$pool_number,
                          prevalence = thetaa,
                          correlation = correlation,
                          sensitivity = design$sensitivity,
                          specificity = design$specificity,
                          form = form))[1,1]
  fi0 <- gdivinv(theta0)^2 /
    solve(fi_pool_cluster(pool_size = design$pool_size,
                          pool_number = design$pool_number,
                          prevalence = theta0,
                          correlation = correlation,
                          sensitivity = design$sensitivity,
                          specificity = design$specificity,
                          form = form))[1,1]
  
  # Note that the below is correct for either kind of one-sided test, but not for two sided tests
  total_clusters_raw <- ((stats::qnorm(power)/sqrt(fia) + stats::qnorm(1 - sig_level)/sqrt(fi0))/(g(theta0) - g(thetaa)))^2
  
  total_clusters <- ceiling(total_clusters_raw)
  
  power_size_results(
    # sample design
    design = design,
    cluster_number = total_clusters,
    # prevalence
    prev_null = theta0,
    prev_alt = thetaa,
    correlation = correlation,
    # statistical test
    sig_level = sig_level,
    power = power,
    alternative = alternative
  )
}



#' @method sample_size_threshold variable_design 
#' @export

sample_size_threshold.variable_design <- function(design,
                                                  prevalence_null,
                                                  prevalence_alt,
                                                  correlation = 0,
                                                  power = 0.8,
                                                  sig_level = 0.05,
                                                  alternative = 'less',
                                                  form = 'logitnorm',
                                                  link = 'logit',
                                                  max_iter = 10000,
                                                  rel_tol = 1e-6,
                                                  ...){
  
  if(!inherits(design$pool_strat, 'pool_strat')){
    stop('pool design must include a valid pooling strategy of class `pool_strat`',
         ' (not just a pooling strategy family)')}
  
  thetaa <- prevalence_alt
  theta0 <- prevalence_null
  
  if(!(alternative %in% c('less', 'greater'))){
    stop('currently only supports one-sided tests. Valid options for alternative are less and greater')
  }
  
  if(alternative == 'less' & theta0 < thetaa){
    stop('If alternative == "less", then prevalence.altnerative must be less than or equal to prevalence_null' )
  }
  
  if(alternative == 'greater' & theta0 > thetaa){
    stop('If alternative == "greater", then prevalence.altnerative must be greater than or equal to prevalence_null' )
  }
  
  # Get link functions
  g <- g_switch(link)
  gdivinv <- gdivinv_switch(link)
  
  fia <- gdivinv(thetaa)^2 /
    solve(fi_pool_cluster_random(catch_dist = design$catch_dist, 
                                 pool_strat = design$pool_strat, 
                                 prevalence = thetaa,
                                 correlation = correlation,
                                 sensitivity = design$sensitivity,
                                 specificity = design$specificity,
                                 form = form,
                                 max_iter = max_iter,
                                 rel_tol = rel_tol-6))[1,1]
  fi0 <- gdivinv(theta0)^2 /
    solve(fi_pool_cluster_random(catch_dist = design$catch_dist, 
                                 pool_strat = design$pool_strat, 
                                 prevalence = theta0,  #should this be theta0 or thetaa?
                                 correlation = correlation,
                                 sensitivity = design$sensitivity,
                                 specificity = design$specificity,
                                 form = form,
                                 max_iter = max_iter,
                                 rel_tol = rel_tol-6))[1,1]
  
  # Note that the below is correct for either kind of one-sided test, but not for two sided tests
  total_cluster_raw <- ((stats::qnorm(power)/sqrt(fia) + stats::qnorm(1 - sig_level)/sqrt(fi0))/(g(theta0) - g(thetaa)))^2
  
  total_clusters <- ceiling(total_cluster_raw)
  
  power_size_results(  
    # sample design
    design = design,
    cluster_number = total_clusters,
    # prevalence
    prev_null = theta0,
    prev_alt = thetaa,
    correlation = correlation,
    # statistical test
    sig_level = sig_level,
    power = power,
    alternative = alternative,
  )
}


# _detect ----

#' @rdname power_threshold
#' @export

power_detect <- function(design,
                         cluster_number,
                         prevalence,
                         correlation = 0,
                         form = 'logitnorm',
                         ...) {
  
  power <- 1 - detection_errors(design, cluster_number, prevalence,
                                correlation, form, ...)$typeII
  
  power_size_results(
    # sample design
    design = design,
    cluster_number = cluster_number,
    # prevalence
    prev = prevalence,
    correlation = correlation,
    # statistical test
    power = power
  )
  
}

#' @rdname power_threshold
#' @export
sample_size_detect <- function(design, 
                               prevalence,
                               correlation,
                               power, 
                               form,
                               ...) {
  UseMethod("sample_size_detect")
}

#' @method sample_size_detect fixed_design
#' @export

sample_size_detect.fixed_design <- function(design,
                                            prevalence,
                                            correlation = 0,
                                            power = 0.8,
                                            form = "logitnorm",
                                            ...) {
  
  # Get detection probability for a single cluster
  errors <- detection_errors(design,1,prevalence,correlation,form)
  
  total_clusters_raw <- log1p(-power)/log(errors$typeII)
  
  #Could be modularised? Everything below here is almost identical (except for
  #how prevalence, alternative, and sig_level are returned) to
  #sample_size_threshold.fixed_design
  
  total_clusters <- ceiling(total_clusters_raw)
  
  power_size_results(  
    # sample design
    design = design,
    cluster_number = total_clusters,
    # prevalence
    prev = prevalence,
    correlation = correlation,
    # statistical test
    power = power
    )
}


#' @method sample_size_detect variable_design 
#' @export

sample_size_detect.variable_design <- function(design,
                                               prevalence,
                                               correlation = 0,
                                               power = 0.8,
                                               form = 'logitnorm',
                                               ...){
  
  if(!inherits(design$pool_strat, 'pool_strat')){
    stop('pool design must include a valid pooling strategy of class `pool_strat`',
         ' (not just a pooling strategy family)')}
  
  # Get detection probability for a single cluster
  errors <- detection_errors(design,1,prevalence,correlation,form)
  
  total_clusters_raw <- log1p(-power)/log(errors$typeII)
  
  total_clusters <- ceiling(total_cluster_raw)
  
  power_size_results(
    # sample design
    design = design,
    cluster_number = total_clusters,
    # prevalence
    prev = prevalence,
    correlation = correlation,
    # statistical test
    power = power
  )
}
