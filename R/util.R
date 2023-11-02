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

