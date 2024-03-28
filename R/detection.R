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
#' @param prevalence_null,prevlanece_alt numeric The proportion of units that
#'   carry the marker of interest (i.e. true positive). `prevalence_null` is the
#'   threshold to compare to and `prevlanece_alt` is the design prevalence. Must
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
#' @param sig_level numeric Signifigance level for statistical test. Defaults to
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
#' power_pool(sample_size = 1000, pool_size = 10,
#'            pool_number = 2, prevalence_null = 0.01,
#'            prevalence_alt = 0.02, correlation = 0.01)
#'
#' sample_size_pool(pool_size = 10, pool_number = 2,
#'                 prevalence_null = 0.01, prevalence_alt = 0.02,
#'                 correlation = 0.01)





detection_errors <- function(prevalence, pool_size,N,M,sensitivity, specificity, correlation,
                             periods_per_location, periods_total,
                             catch.mean, catch.dispersion,
                             form = 'beta', link = NULL){
  rho <- correlation
  
  link <- switch(form, logitnorm = qlogis, cloglognorm = cloglog, beta = function(x){x})
  invlink <- switch(form, logitnorm = plogis, cloglognorm = cloglog_inv, beta = function(x){x})
  
  if(form %in% c('logitnorm', 'cloglognorm')){
    pars <- mu_sigma_linknorm(prevalence,prevalence * (1- prevalence) * rho, link, invlink)
    mu <- pars[1]
    sigma <- pars[2]
    density <- function(x){dnorm(x, mean = mu, sd = sigma)}
  }
  if(form == 'beta'){
    Alpha <- prevalence * (rho^-1 -1)
    Beta <- (1-prevalence) * (rho^-1 -1)
    density <- function(x){dbeta(x,Alpha, Beta)}
  }
  
  if(missing(N) & missing(periods_per_location) & missing(periods_total)){
    stop('One of the following must be provided:
             N (the number of groups per location)
             periods_per_location (the number of sampling periods per location)
             periods_total (total the number of sampling periods across all locations)')
  }
  if(missing(N) & missing(periods_per_location)){
    periods_per_location <- periods_total/M
    if(periods_per_location%%1) warning('Inputs imply a fractional number of sampling periods per sampling location')
  }
  
  if(rho == 0){
    if(missing(N)){ #Case where we assume random (negative binomial) catch sizes at each location
      warning('For correlation = 0, a heirarchical/cluster survey design with M locations and p sampling periods per location is approximately equivalent a simple random survey with p*M sampling periods per location')
      const <- catch.dispersion/(catch.mean + catch.dispersion)
      q <- (1 - (1 - sensitivity - specificity) * (1-prevalence)^pool_size - sensitivity) ^ (1/pool_size)
      typeII <- (const/(1 - q * (1 - const)))^(M * periods_per_location * catch.dispersion)
    }else{
      warning('For correlation = 0, a heirarchical/cluster survey design with M locations and N groups per location is equivalent a simple random survey with N*M groups')
      typeII <- (1 - (1 - sensitivity - specificity) * (1-prevalence)^pool_size - sensitivity)^(N*M)
    }
  }else{
    if(missing(N)){ #Case where we assume random (negative binomial) catch sizes at each location
      f <- function(x){
        q <- (1 - (1 - sensitivity - specificity) * (1-invlink(x))^pool_size - sensitivity) ^ (1/pool_size)
        density(x) *
          (const/(1 - q * (1 - const)))^(periods_per_location * catch.dispersion)
      }
      typeI <- 1 - (const/(1 - specificity^(1/pool_size) * (1 - const)))^(periods_per_location * catch.dispersion * M)
      
    }else{ #Case with fixed number of pools per site
      f <- function(x){
        density(x) *
          (1 - (1 - sensitivity - specificity) * (1-invlink(x))^pool_size - sensitivity)^N
      }
    }
    #typeI <-   1- specificity^(N * M)
    typeI <- -expm1(log(specificity)*N * M) # equivalent to the above commented code, but more numerically stable for high specificity and large N*M
  }
  
  lb <- switch(form, logitnorm = -Inf, cloglognorm = -Inf, beta = 0)
  ub <- switch(form, logitnorm = Inf, cloglognorm = Inf, beta = 1)
  
  
  if(form == 'beta'){
    if(missing(N)){stop('Have not implemented negative binomial sample size with form = beta. It has a nice closed form solution in terms of hypergeometric functions for the case with a perfec test. See paper notes')}
    if(sensitivity ==1){
      typeII <- exp((log(specificity) * N + lbeta(Alpha, Beta + pool_size * N) -  lbeta(Alpha, Beta)) * M)
    }else{
      z <- 0:N
      summand <- (1-sensitivity)^N/beta(Alpha, Beta) * choose(N,z) * ((1-sensitivity - specificity)/(sensitivity - 1))^z * beta(Alpha, Beta + z * pool_size)
      typeII <- sum(summand) ^ M
    }
  }else{
    typeII <- stats::integrate(f, lb, ub)$value ^ M
  }
  list(typeI = typeI, typeII = typeII)
}
