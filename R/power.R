#' Power and sample size calculations for estimating population prevalence from
#' pooled samples
#'
#' `power_pool()` calculates the statistical power of a pooled survey design to
#' determine whether population prevalence is different from a threshold.
#' `sample_size_pool()` calculate the sample size required for a pooled survey
#' to achieve a specified power.
#'
#' @param sample_size numeric The total number of units across the whole sample.
#'   Should be greater than `pool_size * pool_number`
#' @param pool_size numeric The number of units per pool. Must be a numeric
#'   value or vector of values greater than 0.
#' @param pool_number numeric The number of pools per cluster. Must be a integer
#'   value or a vector of integer values greater than or equal to 1.
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
#'
#'
#' @return The statistical power of the proposed design with regards to
#'   comparing prevalence to a threshold (`power_pool()`)
#' @export
#'
#' @examples
#' power_pool(sample_size = 1000, pool_size = 10, pool_number = 2,
#'            prevalence_null = 0.01, prevalence_alt = 0.02)
#'
#' sample_size_pool(pool_size = 10, pool_number = 2,
#'                  prevalence_null = 0.01, prevalence_alt = 0.02)
#'
#' power_pool(sample_size = 1000, pool_size = 10, pool_number = 2,
#'            prevalence_null = 0.01, prevalence_alt = 0.02,
#'            correlation = 0.01)
#'
#' sample_size_pool(pool_size = 10, pool_number = 2,
#'                  prevalence_null = 0.01, prevalence_alt = 0.02,
#'                  correlation = 0.01)


power_pool <- function(sample_size, pool_size, pool_number, prevalence_null, prevalence_alt,
                       correlation = 0, sensitivity = 1, specificity = 1,
                       sig_level = 0.05, alternative = 'greater',
                       form = 'beta', link = 'logit'){
  if(sample_size < pool_size*pool_number){
    stop('The total number of units in sample (sample_size) must be greater than the pool size times pool number (pool_size * pool_number)')
  }
  thetaa <- prevalence_alt
  theta0 <- prevalence_null
  
  # The idea here was that the hypothesis test would be for mu rather than theta, with
  # the idea that this would be equivalent to a test on theta, i.e. if mu0 =
  # g(theta0)<g(thetaa) = mua then theta0 < thetaa. However, this is not true
  # since differences in rho/correlations allows for theta0 > thetaa
  
  # if(real.scale & form %in% c('logitnorm', 'cloglognorm')){ g <- function(x){
  # #calculate mu from theta and rho .var <- correlation * x * (1-x)
  # mu_sigma_linknorm(x,.var, link = switch(form, logitnorm = stats::qlogis,
  # cloglognorm = cloglog), invlink = switch(form, logitnorm = plogis,
  # cloglognorm = cloglog_inv))[1] } gdivinv <- function(x){1}
  # }else{
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
  # }
  
  fia <- sample_size/(pool_size*pool_number) * gdivinv(thetaa)^2 /
    solve(fi_pool_cluster(pool_size = pool_size,
                          pool_number = pool_number,
                          prevalence = thetaa,
                          correlation = correlation,
                          sensitivity = sensitivity,
                          specificity = specificity,
                          form = form))[1,1]
  #print(fia)
  fi0 <- sample_size/(pool_size*pool_number) * gdivinv(theta0)^2 /
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
  power
}

#' @rdname power_pool
#' @export

sample_size_pool <- function(pool_size, pool_number,
                             prevalence_null, prevalence_alt,
                             correlation = 0, sensitivity = 1, specificity = 1,
                             power = 0.8, sig_level = 0.05,
                             alternative = 'greater',
                             form = 'beta',
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
  
  unit_fia <- 1/(pool_size*pool_number) * gdivinv(thetaa)^2 /
    solve(fi_pool_cluster(pool_size = pool_size,
                          pool_number = pool_number,
                          prevalence = thetaa,
                          correlation = correlation,
                          sensitivity = sensitivity,
                          specificity = specificity,
                          form = form))[1,1]
  unit_fi0 <- 1/(pool_size*pool_number) * gdivinv(theta0)^2 /
    solve(fi_pool_cluster(pool_size = pool_size,
                          pool_number = pool_number,
                          prevalence = theta0,
                          correlation = correlation,
                          sensitivity = sensitivity,
                          specificity = specificity,
                          form = form))[1,1]
  
  # Note that the below is correct for either kind of one-sided test, but not for two sided tests
  ((stats::qnorm(power)/sqrt(unit_fia) + stats::qnorm(1 - sig_level)/sqrt(unit_fi0))/(g(theta0) - g(thetaa)))^2
}