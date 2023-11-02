
power_fi <- function(n, s, N, prevalence.null, prevalence.alternative,
                     sig.level = 0.05, alternative = 'greater',
                     sensitivity = 1, specificity = 1,
                     correlation = 0, form = 'beta',
                     link = 'identity'){
  thetaa <- prevalence.alternative
  theta0 <- prevalence.null

  # The idea here was that the hypothesis test would be for mu rather than theta, with
  # the idea that this would be equivalent to a test on theta, i.e. if mu0 =
  # g(theta0)<g(thetaa) = mua then theta0 < thetaa. However, this is not true
  # since differences in rho/correlations allows for theta0 > thetaa

  # if(real.scale & form %in% c('logitnorm', 'cloglognorm')){ g <- function(x){
  # #calculate mu from theta and rho .var <- correlation * x * (1-x)
  # mu_sigma_linknorm(x,.var, link = switch(form, logitnorm = qlogis,
  # cloglognorm = cloglog), invlink = switch(form, logitnorm = plogis,
  # cloglognorm = cloglog_inv))[1] } gdivinv <- function(x){1}
  # }else{
  g <- switch(link,
              logit = qlogis,
              cloglog = cloglog,
              log = log,
              identity = function(x){x})

  gdivinv <- switch(link,
                    logit = function(x){x * (1-x)},
                    cloglog = function(x){-log1p(-x) * (1-x)},
                    log = function(x){x},
                    identity = function(x){1})
  # }

  fia <- n/(s*N) * gdivinv(thetaa)^2 /
    solve(fi_pool_imperfect_cluster(s = s, N = N, p = thetaa,
                                    sensitivity = sensitivity,
                                    specificity = specificity,
                                    correlation = correlation,
                                    form = form))[1,1]
  #print(fia)
  fi0 <- n/(s*N) * gdivinv(theta0)^2 /
    solve(fi_pool_imperfect_cluster(s = s, N = N, p = theta0,  #should this be theta0 or thetaa?
                                    sensitivity = sensitivity,
                                    specificity = specificity,
                                    correlation = correlation,
                                    form = form))[1,1]
  #print(fi0)

  power <- switch(alternative,
                  less = pnorm(((g(theta0) - g(thetaa))  - qnorm(1-sig.level)/sqrt(fi0)) * sqrt(fia)),
                  greater = pnorm(((g(thetaa) - g(theta0))  - qnorm(1-sig.level)/sqrt(fi0)) * sqrt(fia)),
                  two.sided = pnorm(((g(theta0) - g(thetaa))  - qnorm(1-sig.level/2)/sqrt(fi0)) * sqrt(fia)) +
                    pnorm(((g(thetaa) - g(theta0))  - qnorm(1-sig.level/2)/sqrt(fi0)) * sqrt(fia)),
                  stop('invalid alternative. options are less, greater, and two.sided')
  )
  power
}

sample_size_prevalence <- function(s, N, prevalence.null, prevalence.alternative,
                                   power = 0.8, sig.level = 0.05, alternative = 'greater',
                                   sensitivity = 1, specificity = 1,
                                   correlation = 0, form = 'beta',
                                   link = 'identity'){
  thetaa <- prevalence.alternative
  theta0 <- prevalence.null

  if(!(alternative %in% c('less', 'greater'))){
    stop('currently only supports one-sided tests. Valid options for alternative are less and greater')
  }

  if(alternative == 'less' & theta0 < thetaa){
    stop('If alternative == "less", then prevalence.altnerative must be less than or equal to prevalence.null' )
  }

  if(alternative == 'greater' & theta0 > thetaa){
    stop('If alternative == "greater", then prevalence.altnerative must be greater than or equal to prevalence.null' )
  }

  g <- switch(link,
              logit = qlogis,
              cloglog = cloglog,
              log = log,
              identity = function(x){x})

  gdivinv <- switch(link,
                    logit = function(x){x * (1-x)},
                    cloglog = function(x){-log1p(-x) * (1-x)},
                    log = function(x){x},
                    identity = function(x){1})

  unit_fia <- 1/(s*N) * gdivinv(thetaa)^2 /
    solve(fi_pool_imperfect_cluster(s = s, N = N, p = thetaa,
                                    sensitivity = sensitivity,
                                    specificity = specificity,
                                    correlation = correlation,
                                    form = form))[1,1]
  unit_fi0 <- 1/(s*N) * gdivinv(theta0)^2 /
    solve(fi_pool_imperfect_cluster(s = s, N = N, p = theta0,
                                    sensitivity = sensitivity,
                                    specificity = specificity,
                                    correlation = correlation,
                                    form = form))[1,1]

  # Note that the below is correct for either kind of one-sided test, but not for two sided tests
  ((qnorm(power)/sqrt(unit_fia) + qnorm(1 - sig.level)/sqrt(unit_fi0))/(g(theta0) - g(thetaa)))^2
}


#THIS WHOLE FUNCTION COULD PROBABLY BE SIMPLIFIED TO A SPECIAL CASE OF optimise_s_prevalence
optimise_s_prevalence_power <- function(prevalence.null, prevalence.alternative,
                                        cost.unit, cost.pool, cost.location = 0,
                                        sig.level = 0.05,  power = 0.8,
                                        alternative = 'less', link = 'logit',
                                        correlation = 0, N = 1, form = 'beta',
                                        sensitivity = 1, specificity = 1,
                                        max.s = 200, interval = 0){
  ss <- function(s){
    sample_size_prevalence(s, prevalence.null, prevalence.alternative,
                           sig.level,  power,
                           alternative, link,
                           sensitivity, specificity,
                           correlation = correlation,
                           N = N, form = form)
  }

  cost <- function(s){
    n <- ss(s)
    out <- n * (cost.unit + cost.pool/s + cost.location/(s*N))
    out
  }

  opt <- optimise(cost,c(1,max.s))

  cost_opt_ceiling <- cost(ceiling(opt$minimum))
  cost_opt_floor   <- cost(floor(opt$minimum))

  if(cost_opt_ceiling < cost_opt_floor){
    optimal_s <- ceiling(opt$minimum)
    opt_cost <- cost_opt_ceiling
  }else{
    optimal_s <- floor(opt$minimum)
    opt_cost <- cost_opt_floor
  }

  optimal_n <- ceiling(ss(optimal_s)/optimal_s) * optimal_s

  if(optimal_s == max.s) warning('Maximum cost effectivness is achieved at or above the maximum size of pools allowed. Consider increasing max.s')

  if(interval == 0){
    return(list(optimal_s = optimal_s, optimal_cost = opt_cost, optimal_n = optimal_n))
  }else if(interval > 0){
    max_cost <- opt_cost * (1 + interval)
    if(cost(1) < max_cost){
      lower <- 1
    } else{
      lower <- ceiling(uniroot(function(x){cost(x) - max_cost},c(1,optimal_s))$root)
    }
    if(cost(max.s) < max_cost){
      upper <- max.s
      warning('A pool size greater than max.s may fall within the specified range of cost effectiveness. Consider increasing max.s')
    }else{
      upper <- floor(uniroot(function(x){cost(x) - max_cost},c(optimal_s,max.s))$root)
    }

    cost_interval <- c(cost(lower),cost(upper))
    n_interval <- c(ss(lower), ss(upper))
    out <- list(optimal_s = optimal_s,
                optimal_cost = opt_cost,
                optimal_n = optimal_n,
                s_interval = c(lower, upper),
                cost_interval = cost_interval,
                n_interval = n_interval)
    return(out)
  }else{
    stop('interval must be 0 (for no interval) or greater than zero.')
  }
}



power_fi(1000,10,0.01,0.005, specificity = 0.999, link = 'identity')

theta.null <- 0.01
theta.alt  <- 0.005
s <- 10
N <- 2
rho <- 0.1

plot(function(x){binom::binom.power(x,alternative = 'less',n = 750,p = theta.null,method = 'asym')}, to = 0.01, n = 100000)
plot(function(x){binom::binom.power(x,alternative = 'less',n = 750,p = theta.null,method = 'logit')}, to = 0.01, n = 100000,add = TRUE)

de <- design_effect_cluster_fisher(s,N,theta.null,1,1,rho, form = 'beta')
de
ss.indi.arcsin <- pwr::pwr.p.test(h = pwr::ES.h(theta.alt,theta.null), power = 0.8, sig.level = 0.05, alternative = 'less')
ss.indi.arcsin
ss.pool.asympt <- sample_size_prevalence(s,N,theta.null,theta.alt,0.8,sig.level = 0.05,alternative = 'less',correlation = 0.1)
ss.pool.asympt
ss.indi.arcsin$n * de





### Older plots
plot(function(x){Vectorize(power_fi)(7500,10,2,prevalence.null = 0.01, prevalence.alternative = x, alternative = 'less' ,specificity = 1, correlation = 0.1)},  from = 0.01/100, to = 0.01, n = 100, ylim = c(0,1))
plot(function(x){Vectorize(power_fi)(7500,10,2,prevalence.null = 0.01, prevalence.alternative = x, alternative = 'less' ,specificity = 1, correlation = 0.1)},  from = 0.01/100, to = 0.01, n = 100, add= TRUE)

plot(function(x){Vectorize(power_fi)(750,10,2,prevalence.null = 0.01, prevalence.alternative = x, alternative = 'greater' ,specificity = 1, correlation = 0.1)},  from = 0.01, to = 0.03, n = 100, ylim = c(0,1))
plot(function(x){Vectorize(power_fi)(750,10,2,prevalence.null = 0.01, prevalence.alternative = x, alternative = 'greater' ,specificity = 1, correlation = 0.1)},  from = 0.01, to = 0.03, n = 100, add= TRUE)


plot(function(x){Vectorize(power_fi)(750,1,prevalence.null = 0.01, prevalence.alternative = x, link = 'identity', specificity = 0.9)}, to = 0.01, n = 100000, col = 'red',add = TRUE)
plot(function(x){Vectorize(power_fi)(750,1,prevalence.null = 0.01, prevalence.alternative = x, link = 'logit', specificity = 0.9)}, to = 0.01, n = 100000, col = 'red', add = TRUE)

plot(function(x){Vectorize(power_fi)(750,5,prevalence.null = 0.01, prevalence.alternative = x, link = 'identity', specificity = 0.9)}, to = 0.01, n = 100000, col = 'green',add = TRUE)
plot(function(x){Vectorize(power_fi)(750,5,prevalence.null = 0.01, prevalence.alternative = x, link = 'logit', specificity = 0.9)}, to = 0.01, n = 100000, col = 'green', add = TRUE)

plot(function(x){Vectorize(power_fi)(750,1,prevalence.null = 0.01, prevalence.alternative = x, link = 'identity', specificity = 1)}, to = 0.01, n = 100000, col = 'blue',add = TRUE)
plot(function(x){Vectorize(power_fi)(750,1,prevalence.null = 0.01, prevalence.alternative = x, link = 'logit', specificity = 1)}, to = 0.01, n = 100000, col = 'blue', add = TRUE)


plot(function(x){Vectorize(power_fi)(750,x,prevalence.null = 0.01, prevalence.alternative = 0.005, link = 'identity',specificity = 0.9)}, from = 1, to = 100, n = 100, ylim = c(0,1))
plot(function(x){Vectorize(power_fi)(750,x,prevalence.null = 0.01, prevalence.alternative = 0.005, link = 'logit',specificity = 0.9)}, from =1, to = 100, n = 100, add= TRUE)




plot(function(x){power_fi(x,10,prevalence.null = 0.01, prevalence.alternative = 0.005, link = 'logit', sensitivity = 1, specificity = 1, alternative = 'less')}, from = 1, to = 3300, n = 1000, ylim = c(0, 1))
plot(function(x){power_fi(x,10,prevalence.null = 0.01, prevalence.alternative = 0.02, link = 'logit', sensitivity = 1, specificity = 1, alternative = 'greater')}, from = 1, to = 3300, n = 1000, col = 'red', add = TRUE)
plot(function(x){rep(0.8, length(x))}, add = T, from = 1, to = 3300)
sample_size_prevalence(1,0.01,0.005)
sample_size_prevalence(1,0.01,0.02, alternative = 'greater')

plot(function(x){n <- Vectorize(sample_size_prevalence)(x,2,prevalence.null = 0.01, prevalence.alternative = 0.005, correlation = 0.1, link = 'logit', sensitivity = 1, specificity = 1, alternative = 'less')}, from = 1, to = 40, n = 40)
plot(function(x){n <- Vectorize(sample_size_prevalence)(x,2,prevalence.null = 0.01, prevalence.alternative = 0.015, correlation = 0.1, link = 'logit', sensitivity = 1, specificity = 1, alternative = 'greater')}, from = 1, to = 40, n = 20, col = 'red', add = TRUE)


Vectorize(optimise_s_prevalence_power)(0.01,seq(0.001,0.009,0.001),1,1.5, specificity = 0.99,correlation = 0, N = 1,alternative = 'less', interval = 0.05, link= 'logit')

