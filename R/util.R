# Basic and helper functions that won't get exposed in the package

# Link and inverse link functions
cloglog <- function(x){
  log(-log1p(- x))
}
cloglog_inv <- function(x){
  - expm1(-exp(x))
}

# Calculating parameters for link-normal style distributions

mu_sigma_linknorm <- function(.mean, .var,link,invlink){
  if(.var >= .mean * (1- .mean)){
    stop('a distribution on [0,1] cannot have a variance this large')
  }
  init <- c(link(.mean),2)
  abs.tol <- 0
  f <- function(x){
    mu <- x[1]
    sigma <- abs(x[2])
    fExp <- function(x) invlink(x) * dnorm(x, mean = mu, sd = sigma)
    ..exp <- integrate(fExp, -Inf, Inf, abs.tol = abs.tol)$value
    fVar <- function(x) (invlink(x) - ..exp)^2 * dnorm(x, mean = mu,sd = sigma)
    ..var <- integrate(fVar, -Inf, Inf, abs.tol = abs.tol)$value
    sum(abs(c(..exp,..var)/c(.mean,.var) - 1))
  }
  out <- optim(init,f,control = list(reltol = .Machine$double.eps ^ 0.7))$par
  out[2] <- abs(out[2])
  out
}

## No longer actually used: but this help you multiply out the product over k of
## (1-(1-p)^s_k)^y_k as polynomial in (1-p) which lets you write out likelihood
## in terms of sums of beta functions.  The output is the coefficients of the
## polynomial starting at c_0
det_poly <- function(s,y){
  if(length(s) != length(y) || !all((c(s,y) %% 1) == 0)){
    stop('s and y must be integer vectors of common length')
  }

  K <- length(s)
  out <- 1
  for(k in 1:K){
    out <- out * polynom::polynomial(c(1,rep(0,s[k]-1),-1))^y[k]
  }
  coef(out)
}


## Functions for generating catch distribution objects
nb_catch <- function(mean,var){
  out <- list(pmf = function(catch){dnbinom(catch, size = mean^2/(var - mean), mu = mean)},
              min = 0,
              max = Inf,
              pars = c(mean = mean, var = var))
  class(out) <- 'catch_distribution'
  out
}

## Functions for implementing common pooling strategies

pool_equal_s <- function(s){
  function(catch){
    if(catch<s){
      return(list(s = catch, N = 1))
    }else if(catch%%s == 0){
      return(list(s = s, N = catch/s))
    }else{
      return(list(s = c(catch%%s, s), N = c(1, catch%/%s)))
    }
  }
}

pool_equal_N <- function(N){
  if(N%%1 !=0) stop('N must be an integer')
  function(catch){
    if(catch<N){
      return(list(s = 1, N = catch))
    }else if(catch%%N == 0){
      return(list(s = catch/N, N = N))
    }else{
      base_s <- catch%/%N
      base_N <- N - (catch %% N)
      return(list(s = c(base_s, base_s + 1), N = c(base_N, N-base_N)))
    }
  }
}

pool_fraction <- function(frac){
  if(frac<=0 | frac>1) stop('frac must be greater than 0 or less than or equal to 1')
  function(catch){
    s <- c()
    count <- 0
    while(catch>0){
      count <- count + 1
      s[count] <- ceiling(frac * catch)
      catch <- catch - s[count]
    }
    s <- table(s)
    N <- unname(s)
    s <- as.numeric(names(s))
    list(s = s, N = N)
  }
}

pool_fraction_maxs <- function(maxs){
  function(catch){
    frac <- maxs/catch
    s <- c()
    count <- 0
    while(catch>0){
      count <- count + 1
      s[count] <- ceiling(frac * catch)
      catch <- catch - s[count]
    }
    s <- table(s)
    N <- unname(s)
    s <- as.numeric(names(s))
    list(s = s, N = N)
  }
}
