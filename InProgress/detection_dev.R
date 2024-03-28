#' Probability of no detections or false detections in a heirarchical survey and group sampling and an imperfect test
#'
#' @export
#' @param theta The prevalence in the whole population
#' @param s The number of units per group/pool
#' @param N The number of groups/pools tested per location
#' @param M The number of locations
#' @param sensitivity The sensitivity of the *group* test -- assumed the same for all group sizes
#' @param specificity The specificity of the *group* test -- assumed the same for all group sizes
#' @return The probability that the survey will not have any positive tests
#'
#' @seealso
#' @example
#' @references


detection_errors <- function(theta,s,N,M,sensitivity, specificity, correlation,
                             periods_per_location, periods_total,
                             catch.mean, catch.dispersion,
                             form = 'beta', link = NULL){
  rho <- correlation

  link <- switch(form, logitnorm = qlogis, cloglognorm = cloglog, beta = function(x){x})
  invlink <- switch(form, logitnorm = plogis, cloglognorm = cloglog_inv, beta = function(x){x})

  if(form %in% c('logitnorm', 'cloglognorm')){
    pars <- mu_sigma_linknorm(theta,theta * (1- theta) * rho, link, invlink)
    mu <- pars[1]
    sigma <- pars[2]
    density <- function(x){dnorm(x, mean = mu, sd = sigma)}
  }
  if(form == 'beta'){
    Alpha <- theta * (rho^-1 -1)
    Beta <- (1-theta) * (rho^-1 -1)
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
      q <- (1 - (1 - sensitivity - specificity) * (1-theta)^s - sensitivity) ^ (1/s)
      typeII <- (const/(1 - q * (1 - const)))^(M * periods_per_location * catch.dispersion)
    }else{
      warning('For correlation = 0, a heirarchical/cluster survey design with M locations and N groups per location is equivalent a simple random survey with N*M groups')
      typeII <- (1 - (1 - sensitivity - specificity) * (1-theta)^s - sensitivity)^(N*M)
    }
  }else{
    if(missing(N)){ #Case where we assume random (negative binomial) catch sizes at each location
      f <- function(x){
        q <- (1 - (1 - sensitivity - specificity) * (1-invlink(x))^s - sensitivity) ^ (1/s)
        density(x) *
          (const/(1 - q * (1 - const)))^(periods_per_location * catch.dispersion)
      }
      typeI <- 1 - (const/(1 - specificity^(1/s) * (1 - const)))^(periods_per_location * catch.dispersion * M)

    }else{ #Case with fixed number of pools per site
      f <- function(x){
        density(x) *
          (1 - (1 - sensitivity - specificity) * (1-invlink(x))^s - sensitivity)^N
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
      typeII <- exp((log(specificity) * N + lbeta(Alpha, Beta + s * N) -  lbeta(Alpha, Beta)) * M)
    }else{
      z <- 0:N
      summand <- (1-sensitivity)^N/beta(Alpha, Beta) * choose(N,z) * ((1-sensitivity - specificity)/(sensitivity - 1))^z * beta(Alpha, Beta + z * s)
      typeII <- sum(summand) ^ M
    }
  }else{
    typeII <- stats::integrate(f, lb, ub)$value ^ M
  }
  list(typeI = typeI, typeII = typeII)
}

detection_errors_three_levels <- function(theta,s,N,M1,M2,sensitivity, specificity,
                                          periods_per_location, catch.mean, catch.dispersion,
                                          sigma1, sigma2, distribution = 'logitnormal'){
  if(missing(N) & missing(periods_per_location)){
    stop('N (the number of groups per location) or periods_per_location (the number of sampling periods_per_location) must be provided')
  }
  if(missing(N)){ #Case where we assume random (negative binomial) catch sizes at each location
    mu <- mu_logitnorm(theta,sqrt(sigma1^2 + sigma2^2))
    const <- catch.dispersion/(catch.mean + catch.dispersion)
    f1 <- function(x,m){
      q <- (1 - (1 - sensitivity - specificity) * (1-plogis(x))^s - sensitivity) ^ (1/s)
      dnorm(x, mean = m, sd = sigma1) *
        (const/(1 - q * (1 - const)))^(periods_per_location * catch.dispersion)
    }
    typeI <- 1 - (const/(1 - specificity^(1/s) * (1 - const)))^(periods_per_location * catch.dispersion * M1 * M2)
  }else{ #case where sample size is fixed ahead of time
    if(sigma1 == 0 & sigma2 == 0){
      warning('For sigma1 = sigma2 = 0, a nested heirarchical/cluster survey is equivalent to a simple random survey with N*M1*M2 groups')
      typeII <- (1 - (1 - sensitivity - specificity) * (1-theta)^s - sensitivity)^(N*M1*M2)
    }else{
      if(sigma1 == 0){
        warning('For sigma1 = 0, a nested heirarchical/cluster survey is equivalent to a single-level heirarchical/cluster survey, with N*M1 groups tested at each of M2 locations')
        return(detection_errors(theta,s,N*M1,M2,sensitivity, specificity,sigma2, distribution))
      }
      if(sigma2 == 0){
        warning('For sigma2 = 0, a nested heirarchical/cluster survey is equivalent to a single-level heirarchical/cluster survey, with N groups tested at each of M1*M2 locations')
        return(detection_errors(theta,s,N,M1*M2,sensitivity, specificity,sigma1, distribution))
      }
      mu <- mu_logitnorm(theta,sqrt(sigma1^2 + sigma2^2))
      f1 <- function(x,m){
        dnorm(x, mean = m, sd = sigma1) *
          (1 - (1 - sensitivity - specificity) * (1-plogis(x))^s - sensitivity)^N
      }
    }
    #typeI <-   1- specificity^(N * M1 * M2)
    typeI <- -expm1(log(specificity)*N * M1 * M2) # equivalent to the above commented code, but more numerically stable for high specificity and large N*M
  }

  f2 <- function(x){
    out <- x
    xisextreme <- x == -Inf | x == Inf
    out[xisextreme] <- 0
    out[!xisextreme] <- dnorm(x[!xisextreme], mean = 0, sd = sigma2) *
      sapply(x[!xisextreme], function(.x){
        stats::integrate(f1, -Inf, Inf, m = .x + mu)$value ^ M1
      })
    out
  }
  typeII <- stats::integrate(f2, -Inf, Inf)$value ^ M2
  list(typeI = typeI, typeII = typeII)
}

detection_error_plots <- function(theta,s,N,M,sensitivity, specificity, correlation,
                                  periods_per_location, periods_total,
                                  catch.mean, catch.dispersion,
                                  form){
  argg <- as.list(environment())
  argg <- argg[!sapply(argg,is.name)]
  long_arggs <- map(argg,length) %>% unlist %>% sort(decreasing = TRUE) %>% `[`(.,.>1) %>% names

  argg <- do.call(expand_grid,argg)
  errors <- data.frame()
  for(.n in 1:nrow(argg)){
    errors[.n,c('typeI','typeII')] <- do.call(detection_errors, argg[.n,])
  }

  aes_order <- c('x', 'color','shape')
  aesthetics <- long_arggs %>% `names<-`(.,aes_order[1:length(.)])
  aesthetics <- c(list(y = 'error'), aesthetics)

  argg <- argg %>% mutate(across(!aesthetics[['x']], .fns = as.factor))

  cbind(argg,errors) %>%
    pivot_longer(cols = c(typeI, typeII),names_to = 'type',values_to = 'error') %>%
    ggplot(do.call(aes_string,aesthetics)) +
    geom_point() +
    geom_line() +
    facet_wrap(~type,scales = "free_y") +
    geom_hline(yintercept = 0.05)

}


detection_data <- expand_grid(theta = c(0.01, 0.001), M = c(5,25,50),
                              sensitivity = 1, specificity = 1,
                              correlation = c(0.1,0.3), form = 'beta',
                              catch = 100 * (1:10),
                              cost.unit = 1, cost.pool = 4, cost.location = 10,
                              max.pool = 25) %>%
  mutate(s = min(catch/M,25),
         N = catch/(M * s),
         `total cost` = cost.unit * catch + cost.pool * N * M + cost.location * M) %>%
  rowwise() %>%
  mutate(typeIIerror = detection_errors(theta = theta,s = s,N = N, M = M,
                                        sensitivity = sensitivity,
                                        specificity = specificity,
                                        correlation = correlation, form = form)$typeII,
         `detection probability` = 1 - typeIIerror)

dp <- detection_data %>%
  mutate(ICC  = factor(correlation),
         M = factor(M),
         Prevalence = paste(theta * 100, '%')) %>%
  ggplot(aes(x = catch,y = `detection probability`, color  = M)) +
  geom_line() +
  facet_grid(Prevalence~ICC, labeller = label_both) +
  labs(x = 'Total catch at all sites',
       y = 'Power\n(detection probability)',
       color = 'Number of sampling sites') +
  lims(x = c(0, NA)) +
  theme(strip.placement = 'outside', legend.position = 'bottom',
        text = element_text(size = 12)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,NA))
dp
ggsave('./Figures/DetectionPower.png', dp,
       width = 4, height = 4, unit = 'in',dpi = 'retina')

dp1 <- detection_data %>%
  #pivot_longer(c(`detection probability`, `total cost`), names_to = 'measure', values_to = 'value') %>%
  mutate(correlation  = factor(correlation),
         `Total Catch` = factor(catch),
         `Number of Sites` = factor(M)
  ) %>%
  ggplot(aes(x = theta,y = `detection probability`, color  = correlation)) +
  geom_line() +
  facet_grid(`Total Catch`~`Number of Sites`,labeller = label_both) +
  labs(x = 'Prevalence',
       y = 'Power\n(detection probability)',
       color = 'correlation') +
  lims(y = c(0,NA)) +
  theme(strip.placement = 'outside', legend.position = 'bottom') +
  scale_x_continuous(labels = scales::percent,limits = c(0, NA))
dp1




ICC1 <- 0.1
sigma1 <- sqrt(pi^2/(1/ICC1 - 1)/3); sigma1
ICC2 <- 0.1
sigma2 <- sqrt(pi^2/(1/ICC2 - 1)/3); sigma2
theta <- 0.01
s <- 25
M1 <- 15
M2 <- 1
periods_per_location <- 1
catch.mean <- 100
catch.cv <- 0.5
catch.var <- (catch.mean * catch.cv)^2
catch.disp <- catch.mean^2/(catch.var - catch.mean); catch.disp; qnbinom(c(0.025, 0.5,0.975), size = catch.disp, mu = catch.mean)
N <- catch.mean * periods_per_location / s; N
MeanTests <- N * M1 * M2; MeanTests
replicate <- 1 #the number of times each positive test is checked -- only taken as positive if all are positive
sensitivity <- 1
specificity <- 1

detection_error_plots(theta = 0.01,s = seq(5, 30, by = 5),
                      N = 1, M = M1*M2,
                      sensitivity = 1,
                      specificity = 1,
                      correlation = c(0.1,0.2,0.5), form = 'beta')

detection_error_plots(theta = 0.02,s = seq(25, 25, by = 5),
                      periods_total = 12,
                      M = c(1,2,3,4,6,8,12),
                      sensitivity = 1, specificity = 0.999,
                      catch.disp = catch.disp, catch.mean = catch.mean,
                      sigma = seq(1,3,0.2))

detection_errors(theta = theta,s = s,N = N, M = M1*M2,
                 sensitivity = sensitivity^replicate,
                 specificity = 1-(1-specificity)^replicate,
                 correlation = ICC1 + ICC2)


detection_errors(theta = theta,s = s, M = M1*M2,
                 sensitivity = sensitivity^replicate,
                 specificity = 1-(1-specificity)^replicate,
                 periods_per_location = periods_per_location, catch.mean = catch.mean, catch.dispersion = catch.disp,
                 sigma = sqrt(sigma1^2+sigma2^2))


detection_errors_three_levels(theta = theta,s = s, N = N, M1 = M1, M2 = M2,
                              sensitivity = sensitivity^replicate,
                              specificity = 1-(1-specificity)^replicate,
                              sigma1 = sigma1, sigma2 = sigma2)

detection_errors_three_levels(theta = theta,s = s, M1 = M1, M2 = M2,
                              sensitivity = sensitivity^replicate,
                              specificity = 1-(1-specificity)^replicate,
                              periods_per_location = periods_per_location, catch.mean = catch.mean, catch.dispersion = catch.disp,
                              sigma1 = sigma1, sigma2 = sigma2)



mu_logitnorm <- function(mean, sigma){
  if(mean == 0.5){
    return(0)
  }else if(mean < 0.5){
    upper <- qlogis(mean) / (sigma + 1)
    lower <- qlogis(mean) * (sigma + 1)
  }else if(mean > 0.5){
    lower <- qlogis(mean) / (sigma + 1)
    upper <- qlogis(mean) * (sigma + 1)
  }
  f <- function(x){logitnorm::momentsLogitnorm(x,sigma)['mean'] - mean}
  sol <- uniroot(f,c(lower, upper),check.conv = TRUE, tol = .Machine$double.eps^0.25)
  #print(sol$estim.prec)
  return(sol$root)
}




