#' Power and sample size calculations for estimating population prevalence from
#' pooled samples
#'
#' `detection_errors()` calculates the typeI (false detection probability) and
#' typeII error (false non-detection or 1-power of detection) for pool-tested
#' surveys with a known number and size of pools.
#'
#' @param x a sample_design object
#' @param cluster_number numeric The number of clusters. The same sample design,
#'   x, is assumed to be used at each cluster
#' @param prevalence numeric The proportion of units that carry the marker of
#'   interest (i.e. true positive) that you want wish to assess your survey
#'   against
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
#'
#' @return A list with the typeI (false detection probability) and typeII error
#'   (false non-detection or 1-power of detection)
#' @export
#'
#' @examples
#'
#'
#' detection_errors(fixed_design(10, 2), 20,
#'                  prevalence = 0.01,
#'                  sensitivity = 1, specificity = 0.999)
#' 


detection_errors <- function(x, cluster_number, prevalence,
                             correlation = 0, form = 'logitnorm') {
  check_in_range2(prevalence)
  check_in_range2(correlation)
  # No input check for form as done in downstream functions/methods
  UseMethod("detection_errors")
}


#' @rdname detection_errors
#' @export
 
detection_errors.fixed_design <- function(x, cluster_number,
                                          prevalence, correlation = 0,
                                          form = 'logitnorm'){
  
  pool_size <- x$pool_size
  pool_number <- x$pool_number
  sensitivity <- x$sensitivity
  specificity <- x$specificity
  
  if(correlation == 0){ #Non-cluster surveys
    pool_number <- pool_number * cluster_number
    #typeI <- 1 - specificity^pool_number
    #typeII <- (1 - (1 - sensitivity - specificity) * (1-prevalence)^pool_size - sensitivity)^(pool_number)
    
    #equivalent to the above commented code, but more numerically stable for extreme values
    typeI <- -expm1(log(specificity)*pool_number) 
    typeII <- exp(log_one_minus_phi(prevalence, pool_size,sensitivity, specificity) * pool_number)
    
    return(list(typeI = typeI, typeII = typeII))
  }
  
  if(form %in% c('logitnorm', 'cloglognorm')){
    
    link <- switch(form, logitnorm = stats::qlogis, cloglognorm = cloglog)
    invlink <- switch(form, logitnorm = stats::plogis, cloglognorm = cloglog_inv)
    
    pars <- mu_sigma_linknorm(prevalence, prevalence * (1 - prevalence) * correlation,
                              link, invlink)
    mu <- pars[1]
    sigma <- pars[2]
    
    f <- function(x){
      #stats::dnorm(x, mean = mu, sd = sigma) * (1 - (1 - sensitivity - specificity) * (1-invlink(x))^pool_size - sensitivity)^pool_number
      #equivalent to the above, but more numerically stable
      exp(stats::dnorm(x, mean = mu, sd = sigma, log = TRUE) +
            Vectorize(log_one_minus_phi,'p')(invlink(x), pool_size, sensitivity, specificity) * pool_number)
    }
    
    typeII <- stats::integrate(f, -Inf, Inf)$value ^ cluster_number
    
    
  }else if(form == 'beta'){
    Alpha <-      prevalence  * (correlation^-1 - 1)
    Beta  <- (1 - prevalence) * (correlation^-1 - 1)
    
    if(sensitivity == 1){
      typeII <- exp(
        cluster_number * (log(specificity) * pool_number +
                            lbeta(Alpha, Beta + pool_size * pool_number) - 
                            lbeta(Alpha, Beta)) 
      )
    }else{
      z <- 0:pool_number
      summand <- (1-sensitivity)^pool_number/beta(Alpha, Beta) * choose(pool_number,z) *
        ((1-sensitivity - specificity)/(sensitivity - 1))^z *
        beta(Alpha, Beta + z * pool_size)
      typeII <- sum(summand) ^ cluster_number
    }
    
  }else{
    stop(form, 'is not a valid input for form. Supported values are "beta", "logitnorm", and "cloglognorm"')
  }
  
  #typeI <-   1- specificity^(pool_number * cluster_number)
  typeI <- -expm1(log(specificity)*pool_number * cluster_number) # equivalent to the above commented code, but more numerically stable for high specificity and large pool_number*cluster_number
  
  list(typeI = typeI, typeII = typeII)
}

# An old version which combined all the detection_errors functions into one
# Still need to spin out the other functions

# detection_errors <- function(pool_size, pool_number, cluster_number,
#                                      prevalence, correlation,
#                                      sensitivity = 1, specificity = 1,
#                                      periods_per_location, periods_total,
#                                      catch.mean, catch.dispersion,
#                                      form = 'beta', link = NULL){
#   
#   link <- switch(form, logitnorm = stats::qlogis, cloglognorm = cloglog, beta = function(x){x})
#   invlink <- switch(form, logitnorm = stats::plogis, cloglognorm = cloglog_inv, beta = function(x){x})
#   
#   if(form %in% c('logitnorm', 'cloglognorm')){
#     pars <- mu_sigma_linknorm(prevalence, prevalence * (1 - prevalence) * correlation, link, invlink)
#     mu <- pars[1]
#     sigma <- pars[2]
#     density <- function(x){stats::dnorm(x, mean = mu, sd = sigma)}
#   }
#   if(form == 'beta'){
#     Alpha <- prevalence * (correlation^-1 -1)
#     Beta <- (1-prevalence) * (correlation^-1 -1)
#     density <- function(x){dbeta(x,Alpha, Beta)}
#   }
#   
#   if(missing(pool_number) & missing(periods_per_location) & missing(periods_total)){
#     stop('One of the following must be provided:
#              pool_number (the number of groups per location)
#              periods_per_location (the number of sampling periods per location)
#              periods_total (total the number of sampling periods across all locations)')
#   }
#   if(missing(pool_number) & missing(periods_per_location)){
#     periods_per_location <- periods_total/cluster_number
#     if(periods_per_location%%1) warning('Inputs imply a fractional number of sampling periods per sampling location')
#   }
#   
#   if(correlation == 0){
#     if(missing(pool_number)){ #Case where we assume random (negative binomial) catch sizes at each location
#       warning('For correlation = 0, a heirarchical/cluster survey design with cluster_number locations and p sampling periods per location is approximately equivalent a simple random survey with p*cluster_number sampling periods per location')
#       const <- catch.dispersion/(catch.mean + catch.dispersion)
#       q <- (1 - (1 - sensitivity - specificity) * (1-prevalence)^pool_size - sensitivity) ^ (1/pool_size)
#       typeII <- (const/(1 - q * (1 - const)))^(cluster_number * periods_per_location * catch.dispersion)
#     }else{
#       warning('For correlation = 0, a heirarchical/cluster survey design with cluster_number locations and pool_number groups per location is equivalent a simple random survey with pool_number*cluster_number groups')
#       typeII <- (1 - (1 - sensitivity - specificity) * (1-prevalence)^pool_size - sensitivity)^(pool_number*cluster_number)
#     }
#   }else{
#     if(missing(pool_number)){ #Case where we assume random (negative binomial) catch sizes at each location
#       f <- function(x){
#         q <- (1 - (1 - sensitivity - specificity) * (1-invlink(x))^pool_size - sensitivity) ^ (1/pool_size)
#         density(x) *
#           (const/(1 - q * (1 - const)))^(periods_per_location * catch.dispersion)
#       }
#       typeI <- 1 - (const/(1 - specificity^(1/pool_size) * (1 - const)))^(periods_per_location * catch.dispersion * cluster_number)
#       
#     }else{ #Case with fixed number of pools per site
#       f <- function(x){
#         density(x) *
#           (1 - (1 - sensitivity - specificity) * (1-invlink(x))^pool_size - sensitivity)^pool_number
#       }
#     }
#     #typeI <-   1- specificity^(pool_number * cluster_number)
#     typeI <- -expm1(log(specificity)*pool_number * cluster_number) # equivalent to the above commented code, but more numerically stable for high specificity and large pool_number*cluster_number
#   }
#   
#   lb <- switch(form, logitnorm = -Inf, cloglognorm = -Inf, beta = 0)
#   ub <- switch(form, logitnorm = Inf, cloglognorm = Inf, beta = 1)
#   
#   
#   if(form == 'beta'){
#     if(missing(pool_number)){stop('Have not implemented negative binomial sample size with form = beta. It has a nice closed form solution in terms of hypergeometric functions for the case with a perfec test. See paper notes')}
#     if(sensitivity ==1){
#       typeII <- exp((log(specificity) * pool_number + lbeta(Alpha, Beta + pool_size * pool_number) -  lbeta(Alpha, Beta)) * cluster_number)
#     }else{
#       z <- 0:pool_number
#       summand <- (1-sensitivity)^pool_number/beta(Alpha, Beta) * choose(pool_number,z) * ((1-sensitivity - specificity)/(sensitivity - 1))^z * beta(Alpha, Beta + z * pool_size)
#       typeII <- sum(summand) ^ cluster_number
#     }
#   }else{
#     typeII <- stats::integrate(f, lb, ub)$value ^ cluster_number
#   }
#   list(typeI = typeI, typeII = typeII)
# }
