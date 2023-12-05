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
    fExp <- function(x) invlink(x) * stats::dnorm(x, mean = mu, sd = sigma)
    ..exp <- stats::integrate(fExp, -Inf, Inf, abs.tol = abs.tol)$value
    fVar <- function(x) (invlink(x) - ..exp)^2 * stats::dnorm(x, mean = mu,sd = sigma)
    ..var <- stats::integrate(fVar, -Inf, Inf, abs.tol = abs.tol)$value
    sum(abs(c(..exp,..var)/c(.mean,.var) - 1))
  }
  out <- stats::optim(init,f,control = list(reltol = .Machine$double.eps ^ 0.7))$par
  out[2] <- abs(out[2])
  out
}


ev <- function(fn, distr,
               max_iter = 1000,
               rel_tol = 1e-6){
  
  #Helper function for calculating the expected value of a function with respect
  #to a random variable with support on the positive integers. It does not work
  #for more general distributions, e.g. continuous random variables or functions
  #with support on negative and positive numbers. 
  
  supp <- distributions3::support(distr)
  
  if(supp[['min']] < 0) {stop('ev does dot support distributions which can take values on negative numbers')}
  if(distributions3::is_continuous(distr)){stop('ev does dot support continuous distributions')}
  
  #Initialise sum over distribution
  x <- max(0,supp[['min']] - 1)
  max_x <- supp[['max']]
  terminate <- FALSE
  E <- 0
  iter <- 0
  E_incr <- list()
  xs <- c()
  
  
  #Main loop for sum
  while(!terminate){
    x <- x + 1
    mass <- distributions3::pdf(distr,x)
    cumm_mass <- distributions3::cdf(distr,x)
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
    if(all(rel_incr <= rel_tol) | x == max_x | 1 - cumm_mass < .Machine$double.eps*10){
      terminate <- TRUE
    }
    if(iter == max_iter){
      terminate <- TRUE
      warning('Reached max_iter without converging. Increase max_iter')
      plot(xs, E_incr)
    }
  }
  return(E)
}

