

## Functions for generating catch distribution objects

# Perhaps need to work out a proper class system for this? Objects will require a min, max, a pmf function, and a summary of parameters

# Only example so far. Generates function that defines a negative binomial catch distribution based on the mean and variance

nb_catch <- function(mean,var){
  pars <- list(mean = mean, var = var)
  out <- list(pmf = function(catch){dnbinom(catch, size = pars$mean^2/(pars$var - pars$mean), mu = pars$mean)},
              min = 0,
              max = Inf,
              pars = pars)
  class(out) <- 'catch_distribution'
  out
}

## Functions for implementing common pooling strategies


#Generates pooling functions that aim for a fixed maximum pool size (max_s) with any remainder placed in a single smaller pool
pool_max_s <- function(max_s){
  function(catch){
    if(catch<max_s){
      return(list(s = catch, N = 1))
    }else if(catch%%max_s == 0){
      return(list(s = max_s, N = catch/max_s))
    }else{
      return(list(s = c(catch%%max_s, max_s), N = c(1, catch%/%max_s)))
    }
  }
}

#Generates pooling functions that aim for a fixed number (N) of approximately equal sized pools
pool_fixed_N <- function(N){
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

#An example of a more complicated pooling function generator. For the first
#pool, the a fraction (specified) of total catch is placed into a pool (rounding
#up). For the second pool, the same fraction of the remaining units (rounding
#up) is placed into a pool. This is continued until there are no units left. The
#total number of pools grows approximately logarithmically with catch size
#(compare pool_max_s where it grows linearly with catch, and pool_fixed_N where
#it is fixed)
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

#Like pool_fraction, but the fraction is determined by the ratio of the max pool size to catch. This is admitedly a little weird!
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
