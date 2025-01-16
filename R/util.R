# Basic and helper functions that won't get exposed in the package

# Link and inverse link functions
cloglog <- function(x){
  log(-log1p(- x))
}
cloglog_inv <- function(x){
  - expm1(-exp(x))
}

# Helper function giving log of 1 minus the pool positivity probability for a
# given prevalence (p), pool size (s), sensitivity (spec), and specificity
# (spec). Not to be exported
log_one_minus_phi <- function(p, s, sens, spec){
  # log(1 - sens - (1 - spec - sens) * (1 - p)^s)
  if(p %in% 0:1){
    log1p(- (1 - (1 - p)^s) * sens - (1 - p)^s * (1 - spec))
  }else if(sens == 1){
    log1p(-p) * s + log(spec)
  }else{
    log(expm1(log1p(-p) * s) * (sens - 1) + exp(log1p(-p) * s) * spec)
  }
}

# Calculating parameters for link-normal style distributions

mu_sigma_linknorm <- function(.mean, .var,link,invlink){
  max_restarts <- 10
  if(.var >= .mean * (1- .mean)){
    stop('a distribution on [0,1] cannot have a variance this large')
  }
  init <- c(link(.mean),2)
  abs.tol <- 0
  f <- function(x){
    mu <- x[1]
    sigma <- abs(x[2])
    fExp <- function(x) invlink(x+mu) * stats::dnorm(x,sd = sigma)
    ..exp <- stats::integrate(fExp, -10 * sigma, 10 * sigma, abs.tol = abs.tol)$value
    fVar <- function(x) (invlink(x+mu) - ..exp)^2 * stats::dnorm(x,sd = sigma)
    ..var <- stats::integrate(fVar, -10 * sigma, 10 * sigma, abs.tol = abs.tol)$value
    abs(c(..exp,..var)/c(.mean,.var) - 1)
    #c(..exp - .mean, ..var - .var)
  }
  # out <- stats::optim(init,f,method = "SANN", control = list(maxit = 1000))$par
  # out <- stats::optim(out, f, control = list(abstol  = 0, reltol = .Machine$double.eps ^ 0.8))$par
  
  sol <- nleqslv::nleqslv(init, f, control = list(xtol = 0, ftol = 1e-9))
  
  # sol$termcd == 1 means convergence. Anything else means convergence issues.
  
  # sol$termcd == 3 means algorithm has gotten stuck, which can usually fixed by
  # restarting with slightly perturbed initial conditions. In this case we
  # restart max_restart times before giving up
  
  restarts <- 0
  while(restarts < max_restarts){
    if(sol$termcd == 3){
      sol <- nleqslv::nleqslv(sol$x * (runif(2, 0.95, 1.05)),
                              f, control = list(xtol = 0, ftol = 1e-9))
      restarts <- restarts + 1
    }else{
      break
    }
  }

  
  if(sol$termcd != 1){
    stop('No solution for .mean = ', .mean, ' .var = ', .var, ' and provided link. ',
         'Convergence flag is ', sol$termcd, ' and f(x) = ', f(sol$x),'.')
  }

  sol$x
}


ev <- function(fn, distr,
               max_iter = 10000,
               rel_tol = 1e-6){
  
  #Helper function for calculating the expected value of a function with respect
  #to a random variable with support on the positive integers. It does not work
  #for more general distributions, e.g. continuous random variables or functions
  #with support on negative and positive numbers. 
  
  supp <- distr::support(distr)
  
  if(min(supp) < 0) {stop('ev does dot support distributions which can take values on negative numbers')}
  
  #It would be good to also have a check for whether the distribution is
  #continuous. I'm not sure how to do this now that we have moved to distr
  #package
  
  #if(distributions3::is_continuous(distr)){stop('ev does dot support
  #continuous distributions')}
  
  #Initialise sum over distribution
  x <- max(0,min(supp) - 1)
  terminate <- FALSE
  E <- 0
  iter <- 0
  E_incr <- list()
  xs <- c()
  
  
  #Main loop for sum
  while(!terminate){
    x <- x + 1
    mass <- distr::d(distr)(x)
    cumm_mass <- distr::p(distr)(x)
    #this avoids unnecessary calls to fn and prevents the early termination of
    #the algorithm for distributions that may have 0 mass for some n but
    #non-zero mass for some m>n (e.g. if distribution only has mass on multiples
    #of 10)
    if(mass == 0){next} 
    
    # Note that iteration counter comes after check for zero mass: for the
    # purposes of early termination, only counts iteration if mass is non-zero
    iter <- iter + 1 
    
    xs[iter] <- x
    E_incr[[iter]] <- mass * fn(x)
    E <- E +  E_incr[[iter]]
    # Stop if increment changes ALL elements of FI by less than fraction rel_tol
    # OR cumm_mass reaches 1 OR if distribution of catch size has finite support
    # (i.e. if there is a maximum possible catch size)

    rel_incr <- abs(E_incr[[iter]]/E)

    if(all(rel_incr <= rel_tol) | 1 - cumm_mass < .Machine$double.eps*10){
      terminate <- TRUE
    }
    if(iter == max_iter){
      print(rel_incr)
      print(iter)
      terminate <- TRUE
      warning('Reached max_iter without converging. Increase max_iter')
      plot(xs, E_incr)
    }
  }
  return(E)
}

sum_n_rv <- function(rv, n){
  if(!is.numeric(n) || n%%1 || n < 1){
    stop('n must be a positive integer')
  }
  srv <- rv
  if(n > 1){
    for(m in 1:(n-1)){
      srv <- srv + rv
    }
  }
  srv
}

