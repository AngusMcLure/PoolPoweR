#' Tools for defining catch-size distributions
#'
#' `nb_catch()` and `pois_catch()` are convenience functions for defining the
#' distribution of catche sizes (aka cluster sizes) for the negative binomial
#' and Poisson distributions respectively
#'
#' @param mean numeric The mean number of units per cluster (catch)
#' @param variance numeric The variance of the number of units per cluster
#'   (catch)
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
  if(mean>variance){stop('variance must be greater than or equal to mean')}
  distributions3::NegativeBinomial(size = mean^2/(variance - mean), p = mean/variance)
}

#' @rdname catch_distributions
#' @export
pois_catch <- function(mean){
  distributions3::Poisson(lambda = mean)
}

## Functions for implementing common pooling strategies



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
