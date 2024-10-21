#' Tools for defining catch-size distributions
#'
#' `nb_catch()` and `pois_catch()` are convenience functions for defining the
#' distribution of catche sizes (aka cluster sizes) for the negative binomial
#' and Poisson distributions respectively. `fit_catch()` provides a convenient
#' way to estimate a (negative binomial) catch distribution from data
#'
#' @param mean numeric The mean number of units per cluster (catch)
#' @param variance numeric The variance of the number of units per cluster
#'   (catch)
#' @param catch numeric A vector of catch sizes
#'
#' @return An object of class `distr` summarising the distribution of
#'   cluster/catch sizes
#' @rdname catch_distributions
#' @export
#'
#' @examples
#' nb_catch(10,20)
#'
#' pois_catch(10)

nb_catch <- function(mean,variance){
  if(mean>=variance){stop('variance must be greater than mean')}
  distr::Nbinom(size = mean^2/(variance - mean), prob = mean/variance)
}

#' @rdname catch_distributions
#' @export
pois_catch <- function(mean){
  distr::Pois(lambda = mean)
}

#' @rdname catch_distributions
#' @export
fit_catch <- function(catch, distribution = 'nb'){
  if(distribution == 'nb'){
    fit <- fitdistrplus::fitdist(catch, 'nbinom')$estimate
    nb_catch(mean = fit['mu'], variance = fit['mu']^2/fit['size'] + fit['mu'])
  }else if(distribution == 'pois'){
    fit <- fitdistrplus::fitdist(catch, 'pois')$estimate
    pois_catch(fit['lambda'])
  }else{
    stop('only support for two distributions: negative binomial ("nb") and Poisson ("pois")')
  }
}


# ## No longer actually used: but this help you multiply out the product over k of
# ## (1-(1-p)^s_k)^y_k as polynomial in (1-p) which lets you write out likelihood
# ## in terms of sums of beta functions.  The output is the coefficients of the
# ## polynomial starting at c_0
# det_poly <- function(s,y){
#   if(length(s) != length(y) || !all((c(s,y) %% 1) == 0)){
#     stop('s and y must be integer vectors of common length')
#   }
# 
#   K <- length(s)
#   out <- 1
#   for(k in 1:K){
#     out <- out * polynom::polynomial(c(1,rep(0,s[k]-1),-1))^y[k]
#   }
#   coef(out)
# }
