#' Power and sample size calculations for estimating population prevalence from
#' pooled samples
#'
#' `power_pool()` calculates the statistical power of a pooled survey design to
#' determine whether population prevalence is different from a threshold.
#' `sample_size_pool()` calculate the sample size required for a pooled survey
#' to achieve a specified power.
#'
#' @param pool_size numeric The number of units per pool. Must be a numeric
#'   value or vector of values greater than 0.
#' @param pool_number numeric The number of pools per cluster. Must be a integer
#'   value or a vector of integer values greater than or equal to 1.
#' @param cluster_number numeric The total number of clusters in a cluster
#'   survey design (should be greater than 1) or 1 surveys where all collection
#'   happens at a single site or collected via simple random sampling from the
#'   target population.
#' @param catch_dist An object of class `distribution` (e.g. produced by
#'   `nb_catch()`) defining the distribution of the possible catch.
#' @param pool_strat function Defines a rule for how a number of units will be
#'   divided into pools. Must take a single numeric argument and return a named
#'   list of pool sizes and pool numbers. `pool_max_size()` and
#'   `pool_target_number()` provide convenience functions for defining common
#'   pooling strategies.
#' @param prevalence_null,prevalence_alt numeric The proportion of units that
#'   carry the marker of interest (i.e. true positive). `prevalence_null` is the
#'   threshold to compare to and `prevalence_alt` is the design prevalence. Must
#'   be be a numeric value between 0 and 1, inclusive of both.
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
#' @param power numeric The desired statistical power of the survey.
#' @param sig_level numeric Significance level for statistical test. Defaults to
#'   0.05. Must be strictly between 0 and 1.
#' @param alternative string The kind of comparison to make. If "greater"
#'   (default) or "less" computes power of tests for one-sided comparisons that
#'   population prevalence (`prevalence_alt`) is greater than or less than the
#'   reference on threshold prevalence (`prevalence_null`) respectively. If
#'   "two.sided" computes power for a two-sided comparison (`power_pool()`
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
#'
#'
#' @return The statistical power of the proposed design with regards to
#'   comparing prevalence to a threshold (`power_pool()`) or a list with the
#'   sample size (number of clusters, pools, and units) required to achieve
#'   desired power (`sample_size_pool()`)
#' @export
#'
#' @examples
#' power_pool(pool_size = 10, pool_number = 2, cluster_number = 50,
#'            prevalence_null = 0.01, prevalence_alt = 0.02)
#'
#' sample_size_pool(pool_size = 10, pool_number = 2,
#'                  prevalence_null = 0.01, prevalence_alt = 0.02)
#'
#' power_pool(pool_size = 10, pool_number = 2, cluster_number = 50,
#'            prevalence_null = 0.01, prevalence_alt = 0.02,
#'            correlation = 0.01)
#' 
#' sample_size_pool(pool_size = 10, pool_number = 2,
#'                  prevalence_null = 0.01, prevalence_alt = 0.02,
#'                  correlation = 0.01)
#' 
#' power_pool_random(nb_catch(20,25), pool_target_number(2), cluster_number = 50,
#'                   prevalence_null = 0.01, prevalence_alt = 0.02,
#'                   correlation = 0.01)
#' 
#' sample_size_pool_random(nb_catch(20,25), pool_target_number(2),
#'                         prevalence_null = 0.01, prevalence_alt = 0.02,
#'                          correlation = 0.01)



power_pool <- function(pool_size, pool_number, cluster_number,
                       prevalence_null, prevalence_alt,
                       correlation = 0, sensitivity = 1, specificity = 1,
                       sig_level = 0.05, alternative = 'greater',
                       form = 'logitnorm', link = 'logit'){
  if(correlation > 0 & cluster_number <= 1){
    stop('The number of clusters (cluster_number) must be (substantially) greater than 1 if there is non-zero correlation between units in a cluster')
  }
  
  if(correlation > 0 & cluster_number <= 10){
    warning('Estimated power may be unreliable if number of clusters (cluster_number) is less than 10')
  }
  
  if(correlation == 0 & cluster_number * pool_number <= 10){
    warning('Estimated power may be unreliable if total number of pools (cluster_number * pool_number) is less than 10')
  }
  
  thetaa <- prevalence_alt
  theta0 <- prevalence_null
  
  g <- switch(link,
              logit = stats::qlogis,
              cloglog = cloglog,
              log = log,
              identity = function(x){x})
  
  gdivinv <- switch(link,
                    logit = function(x){x * (1-x)},
                    cloglog = function(x){-log1p(-x) * (1-x)},
                    log = function(x){x},
                    identity = function(x){1})
  
  fia <- cluster_number * gdivinv(thetaa)^2 /
    solve(fi_pool_cluster(pool_size = pool_size,
                          pool_number = pool_number,
                          prevalence = thetaa,
                          correlation = correlation,
                          sensitivity = sensitivity,
                          specificity = specificity,
                          form = form))[1,1]
  #print(fia)
  fi0 <- cluster_number * gdivinv(theta0)^2/
    solve(fi_pool_cluster(pool_size = pool_size,
                          pool_number = pool_number,
                          prevalence = theta0,  #should this be theta0 or thetaa?
                          correlation = correlation,
                          sensitivity = sensitivity,
                          specificity = specificity,
                          form = form))[1,1]
  #print(fi0)
  
  power <- switch(alternative,
                  less = stats::pnorm(((g(theta0) - g(thetaa))  - stats::qnorm(1-sig_level)/sqrt(fi0)) * sqrt(fia)),
                  greater = stats::pnorm(((g(thetaa) - g(theta0))  - stats::qnorm(1-sig_level)/sqrt(fi0)) * sqrt(fia)),
                  two.sided = stats::pnorm(((g(theta0) - g(thetaa))  - stats::qnorm(1-sig_level/2)/sqrt(fi0)) * sqrt(fia)) +
                    stats::pnorm(((g(thetaa) - g(theta0))  - stats::qnorm(1-sig_level/2)/sqrt(fi0)) * sqrt(fia)),
                  stop('invalid alternative. options are less, greater, and two.sided')
  )
  
  # Prepare output
  total_pools <- cluster_number * pool_number
  total_units <- total_pools * pool_size
  text <- paste(
    "A survey design using", 
    is_perfect_test(sensitivity, specificity), 
    "diagnostic test on pooled samples with the above parameters has a statistical power of", 
    round(power, 3)
  )
  
  power_size_results(  
    sensitivity = sensitivity,
    specificity = specificity,
    # prevalence
    prev_null = theta0,
    prev_alt = thetaa,
    correlation = correlation,
    # statistical test
    sig_level = sig_level,
    power = power,
    alternative = alternative,
    # sample design
    pool_size = pool_size,
    pool_number = pool_number,
    cluster_number = cluster_number,
    total_pools = total_pools,
    total_units = total_units,
    # parsing
    text = text
  )
  
}

#' @rdname power_pool
#' @export

power_pool_random <- function(catch_dist, pool_strat, cluster_number,
                              prevalence_null, prevalence_alt,
                              correlation = 0, sensitivity = 1, specificity = 1,
                              sig_level = 0.05, alternative = 'greater',
                              form = 'logitnorm', link = 'logit',
                              max_iter = 10000, rel_tol = 1e-6){
  
  exp_pools <- ev(\(catch) sum(pool_strat(catch)$pool_number),
                  catch_dist, max_iter, rel_tol) * cluster_number
  
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
  
  g <- switch(link,
              logit = stats::qlogis,
              cloglog = cloglog,
              log = log,
              identity = function(x){x})
  
  gdivinv <- switch(link,
                    logit = function(x){x * (1-x)},
                    cloglog = function(x){-log1p(-x) * (1-x)},
                    log = function(x){x},
                    identity = function(x){1})
  
  fia <- cluster_number * gdivinv(thetaa)^2 /
    solve(fi_pool_cluster_random(catch_dist = catch_dist, 
                                 pool_strat = pool_strat, 
                                 prevalence = thetaa,
                                 correlation = correlation,
                                 sensitivity = sensitivity,
                                 specificity = specificity,
                                 form = form,
                                 max_iter = max_iter,
                                 rel_tol = rel_tol-6))[1,1]
  
  fi0 <- cluster_number * gdivinv(theta0)^2 /
    solve(fi_pool_cluster_random(catch_dist = catch_dist, 
                                 pool_strat = pool_strat, 
                                 prevalence = theta0,  #should this be theta0 or thetaa?
                                 correlation = correlation,
                                 sensitivity = sensitivity,
                                 specificity = specificity,
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
  
  # Prepare output
  if (sensitivity == 1 && specificity == 1) {
    perf = "a perfect"
  } else {
    perf = "an imperfect"
  }
  text = paste(
    "A survey design using", 
    is_perfect_test(sensitivity, specificity), 
    "diagnostic test on pooled samples with the above parameters has a statistical power of",
    round(power, 3)
  )
  
  power_size_results(  
    sensitivity = sensitivity,
    specificity = specificity,
    # prevalence
    prev_null = theta0,
    prev_alt = thetaa,
    correlation = correlation,
    # statistical test
    sig_level = sig_level,
    power = power,
    alternative = alternative,
    # sample design
    catch_dist = catch_dist,
    pool_strat = as.character(pool_strat),
    cluster_number = cluster_number,
    # parsing
    text = text
  )
}

#' @rdname power_pool
#' @export

sample_size_pool <- function(pool_size, pool_number,
                             prevalence_null, prevalence_alt,
                             correlation = 0, sensitivity = 1, specificity = 1,
                             power = 0.8, sig_level = 0.05,
                             alternative = 'greater',
                             form = 'logitnorm',
                             link = 'logit'){
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
  
  g <- switch(link,
              logit = stats::qlogis,
              cloglog = cloglog,
              log = log,
              identity = function(x){x})
  
  gdivinv <- switch(link,
                    logit = function(x){x * (1-x)},
                    cloglog = function(x){-log1p(-x) * (1-x)},
                    log = function(x){x},
                    identity = function(x){1})
  
  fia <- gdivinv(thetaa)^2 /
    solve(fi_pool_cluster(pool_size = pool_size,
                          pool_number = pool_number,
                          prevalence = thetaa,
                          correlation = correlation,
                          sensitivity = sensitivity,
                          specificity = specificity,
                          form = form))[1,1]
  fi0 <- gdivinv(theta0)^2 /
    solve(fi_pool_cluster(pool_size = pool_size,
                          pool_number = pool_number,
                          prevalence = theta0,
                          correlation = correlation,
                          sensitivity = sensitivity,
                          specificity = specificity,
                          form = form))[1,1]
  
  # Note that the below is correct for either kind of one-sided test, but not for two sided tests
  total_clusters_raw <- ((stats::qnorm(power)/sqrt(fia) + stats::qnorm(1 - sig_level)/sqrt(fi0))/(g(theta0) - g(thetaa)))^2
  
  total_clusters <- ceiling(total_clusters_raw)
  total_pools <- total_clusters * pool_number
  total_units <- total_pools * pool_size
  text <- paste0(
    "A survey design using ", is_perfect_test(sensitivity, specificity), 
    " diagnostic test on pooled samples with the above parameters requires a total of ",
    total_clusters, " clusters, ", 
    total_pools, " total pools, and ", 
    total_units, " total units."
  )
  
  power_size_results(  
    sensitivity = sensitivity,
    specificity = specificity,
    # prevalence
    prev_null = theta0,
    prev_alt = thetaa,
    correlation = correlation,
    # statistical test
    sig_level = sig_level,
    power = power,
    alternative = alternative,
    # sample design
    pool_size = pool_size,
    pool_number = pool_number,
    cluster_number = total_clusters,
    total_pools = total_pools,
    total_units = total_units,
    # parsing
    text = text
  )
}




#' @rdname power_pool
#' @export
#' 

sample_size_pool_random <- function(catch_dist, pool_strat,
                                    prevalence_null, prevalence_alt,
                                    correlation = 0, sensitivity = 1, specificity = 1,
                                    power = 0.8, sig_level = 0.05,
                                    alternative = 'greater',
                                    form = 'logitnorm', link = 'logit',
                                    max_iter = 10000, rel_tol = 1e-6){
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
  
  g <- switch(link,
              logit = stats::qlogis,
              cloglog = cloglog,
              log = log,
              identity = function(x){x})
  
  gdivinv <- switch(link,
                    logit = function(x){x * (1-x)},
                    cloglog = function(x){-log1p(-x) * (1-x)},
                    log = function(x){x},
                    identity = function(x){1})
  
  fia <- gdivinv(thetaa)^2 /
    solve(fi_pool_cluster_random(catch_dist = catch_dist, 
                                 pool_strat = pool_strat, 
                                 prevalence = thetaa,
                                 correlation = correlation,
                                 sensitivity = sensitivity,
                                 specificity = specificity,
                                 form = form,
                                 max_iter = max_iter,
                                 rel_tol = rel_tol-6))[1,1]
  fi0 <- gdivinv(theta0)^2 /
    solve(fi_pool_cluster_random(catch_dist = catch_dist, 
                                 pool_strat = pool_strat, 
                                 prevalence = theta0,  #should this be theta0 or thetaa?
                                 correlation = correlation,
                                 sensitivity = sensitivity,
                                 specificity = specificity,
                                 form = form,
                                 max_iter = max_iter,
                                 rel_tol = rel_tol-6))[1,1]
  
  # Note that the below is correct for either kind of one-sided test, but not for two sided tests
  total_cluster_raw <- ((stats::qnorm(power)/sqrt(fia) + stats::qnorm(1 - sig_level)/sqrt(fi0))/(g(theta0) - g(thetaa)))^2
  
  total_clusters <- ceiling(total_cluster_raw)
  exp_total_units <- round(mean(catch_dist) * total_clusters,1)
  exp_total_pools <- round(ev(\(catch) sum(pool_strat(catch)$pool_number),
                              catch_dist, max_iter, rel_tol) * total_clusters, 1)
  
  return(list(clusters = total_clusters,
              expected_pools = exp_total_pools,
              expected_units = exp_total_units))
  
}