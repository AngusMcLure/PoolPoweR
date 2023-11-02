### MOST OF THESE DON'T MAKE SENSE OR ARE AT BEST POORLY DEFINED SINCE WHENEVER SENSITIVITY IS PERFECT IT IS ALWAYS BETTER TO HAVE AS LARGE A POOL AS POSSIBLE (i.e. ONLY ONE POOL PER SITE)
### These first two functions might be OK, but the rest probably don't work at all (and might be subject to bitrot and imcompatible with the latest detection_error functions. In particular the latest detection-error functions are in terms of correlation not sigma)

optimise_detection <- function(theta, max.typeI, max.typeII,
                               sensitivity, specificity,
                               correlation, form = 'beta',
                               cost.unit, cost.pool, cost.location, init = NULL, init_M = 30,
                               LB = c(1,1,1), UB = c(Inf,Inf,Inf), int.step = 10, int.search.only = F){


  #if we assume that typeI error is 1 - specificity^(N*M) then this can be
  #translated into a simple constraint on N and M
  maxtests <- log(1 - max.typeI)/log(specificity)

  #Since a hierarchical survey is less likely to find a positive than a simple random sample of the same size, we can use this to calculate a minimum number of units that need to be tested
  minunits <- log(max.typeII)/log((1-theta) * specificity + theta * (1 - sensitivity))

  if(sensitivity == 1){

    cost <- function(par){
      s <- par[1]; N <- 1; M <- par[2]
      cost.unit * s*N*M + cost.pool * N*M + cost.location * M
    }

    errors <- function(par){
      unlist(detection_errors(theta,par[1],1,par[2],1,specificity,correlation, form = 'beta'))
    }
  }else{
    cost <- function(par){
      s <- par[1]; N <- par[2]; M <- par[3]
      cost.unit * s*N*M + cost.pool * N*M + cost.location * M
    }

    errors <- function(par){
      unlist(detection_errors(theta,par[1],par[2],par[3],sensitivity,specificity,correlation, form = 'beta'))
    }
  }



  #initial guess
  if(is.null(init)){
    init_M <- max(LB[3] + 1, init_M)
    init_N <- if(sensitivity == 1){1}else{max(LB[2] + 1, floor(maxtests/init_M))}
    init_n <- max(LB[1] + 1, ceiling(minunits/(init_N * init_M)))
    if(sensitivity == 1){
      init <- c(init_n, init_M)
      LB <- LB[c(1,3)]
      UB <- UB[c(1,3)]
    }else{
      init <- c(init_n,init_N, init_M)
    }
  }
  if(int.search.only){
    opt <- list(pars = init)
  }else{
    real_opt <- Rsolnp::solnp(init, cost, ineqfun = errors,
                              ineqLB = c(0,0),
                              ineqUB = c(max.typeI, max.typeII),
                              LB = LB,
                              UB = UB,
                              control = list(delta = 1.0e-5,
                                             tol = 1.0e-3)
    )
    opt <- real_opt
  }
  n <- 0
  converged <- FALSE
  while(!converged){
    n <- n + 1
    print(paste("Local integer search:", n, ". Pars: ",paste(round(opt$pars),collapse = " ")))
    opt_new <- opt_local_int(opt$pars,cost,ineqfun = errors,
                             ineqLB = c(0,0),
                             ineqUB = c(max.typeI, max.typeII),
                             LB = LB,
                             UB = UB,
                             steps = int.step)
    if(all(opt_new$pars == opt$pars)) converged = TRUE
    opt <- opt_new
  }

  names(opt) <- c("design", "cost")
  opt$errors <- errors(opt$design)
  if(sensitivity ==1){
    opt$design <- c(opt$design[1],1,opt$design[2])
  }
  names(opt$design) <- c("s","N","M")
  names(opt$cost) <- NULL
  names(opt$errors) <- c("typeI","typeII")
  opt
}

opt_local_int <- function(pars,fun,ineqfun, ineqLB, ineqUB,LB, UB, steps,max_step_doublings = 2){
  satisfied <- FALSE
  step_doublings <- 0
  while(!satisfied){
    search <- do.call(expand.grid,purrr::map(round(pars),~{(.x-steps):(.x+steps)}))
    opt_val <- Inf
    opt_par <- NULL
    for(n in 1:nrow(search)){
      test_par <- unlist(search[n,])
      val <- fun(test_par)
      if(val < opt_val & all(test_par >= LB) & all(test_par <= UB)){
        val_fineq <- ineqfun(test_par)
        if(all(ineqLB <= val_fineq)  & all(val_fineq <= ineqUB)){
          opt_val <- val
          opt_par <- test_par
        }
      }
    }
    if(is.infinite(opt_val) || is.null(opt_par)){
      if(step_doublings == max_step_doublings){
        stop('Could not find any solutions to inequalities')
      }
      steps <- steps * 2
      step_doublings <- step_doublings + 1
    }else{
      satisfied <- TRUE
    }
  }

  list(pars = opt_par, value = opt_val)
}

optimise_detection_fixed <- function(s = NULL, N = NULL,
                                     theta, max.typeI, max.typeII,
                                     sensitivity, specificity,
                                     correlation, form = 'beta',
                                     cost.unit, cost.pool, cost.location,
                                     LB = c(1,1), UB = c(Inf,Inf),
                                     int.step = 10, int.search.only = F, verbose = F){
  #if we assume that typeI error is 1 - specificity^(N*M) then this can be
  #translated into a simple constraint on N and M
  maxtests <- if(specificity <1){log1p(- max.typeI)/log(specificity)}else{Inf}

  #Since a hierarchical survey is less likely to find a positive than a simple random sample of the same size, we can use this to calculate a minimum number of units that need to be tested
  minunits <- log(max.typeII)/log((1-theta) * specificity + theta * (1 - sensitivity))

  if(!is.null(s) && minunits/s > maxtests){
    stop('It is not possible to achieve both the specified type I error and type II error with the test. You must either:
         * increase the design prevlance
         * use larger poolsize (s)
         * use a test with higher sensitivity and/or specificity
         * or relax (increase) your maximum error requirements.
         minunits = ', minunits,'
         maxtests = ', maxtests)
  }

  cost <- function(pars){
    catchpersite <- pars[1];
    M <- pars[2];
    if(!is.null(s)){
      N <- catchpersite/s
    }
    cost.unit * catchpersite * M + cost.pool * N*M + cost.location * M
  }

  errors <- function(pars){
    catchpersite <- pars[1];
    M <- pars[2];
    if(!is.null(s)){
      N <- catchpersite/s
    }else if(!is.null(N)){
      s <- catchpersite/N
    }
    er <- unlist(detection_errors(theta = theta,
                                  s = s, N = N, M = M,
                                  sensitivity = sensitivity,
                                  specificity = specificity,
                                  correlation = correlation,
                                  form = form))
    er
  }

  init <- c(10, minunits/10)


  if(int.search.only){
    opt <- list(pars = init)
  }else{
    opt <- Rsolnp::solnp(init, cost, ineqfun = errors,
                              ineqLB = c(0,0), ineqUB = c(max.typeI, max.typeII),
                              LB = LB, UB = UB,
                              control = list(delta = 1.0e-5,tol = 1.0e-3)
    )
  }
  n <- 0
  converged <- FALSE
  while(!converged){
    n <- n + 1
    if(verbose){
      print(paste("Local integer search:", n, ". Pars: ",paste(round(opt$pars),collapse = " ")))
    }
    opt_new <- opt_local_int(opt$pars,cost,ineqfun = errors,
                             ineqLB = c(0,0),
                             ineqUB = c(max.typeI, max.typeII),
                             LB = LB, UB = UB,
                             steps = int.step)
    if(all(opt_new$pars == opt$pars)) converged = TRUE
    opt <- opt_new
  }

  names(opt) <- c("design", "cost")
  opt$errors <- errors(opt$design)
  optcatchpersite <- opt$design[1];
  optM <- opt$design[2];
  if(!is.null(s)){
    optN <- optcatchpersite/(s)
    opts <- s
  }else if(!is.null(N)){
    opts <- optcatchpersite/(N)
    optN <- N
  }
  opt$design <- c(optcatchpersite,optM,optN, opts)
  names(opt$design) <- c("catch per site","M","N","s")
  names(opt$cost) <- NULL
  names(opt$errors) <- c("typeI","typeII")
  opt
}

optimise_detection_fixed(N = 1, theta = 0.01, max.typeI =  0.3,
                         max.typeII = 0.05, sensitivity = 1, specificity = 1,
                         cost.unit = 1, cost.pool = 4, cost.location = 10,
                         correlation = 0.1, int.step = 100, int.search.only = T)

optimise_detection_data <- expand_grid(theta = seq(0.0005, 0.01, by = 0.0001), N = 1,
                                       sensitivity = 1, specificity = 1,
                                       correlation = c(0.1,0.3), form = 'beta',
                                       cost.unit = 1, cost.pool = 4, cost.location = c(4,40)) %>%
  rowwise() %>%
  mutate(opt = list(optimise_detection_fixed(N = 1, theta = theta,
                                             max.typeI =  0.2, max.typeII = 0.2,
                                             cost.unit = cost.unit, cost.pool = cost.pool, cost.location = cost.location,
                                             sensitivity = sensitivity, specificity = specificity,
                                             correlation = correlation, form = form,
                                             int.step = 30, int.search.only = T)),
         optcost = opt$cost,
         opts = opt$design['s'], optM = opt$design['M'],
         optcatchpersite = opt$design['catch per site'],
         opttotalcatch = optcatchpersite*optM)

odp <- optimise_detection_data %>%
  mutate(opttotalcatch = optcatchpersite*optM) %>%
  pivot_longer(c(optcost,optM,optcatchpersite,opttotalcatch), names_to = 'measure', values_to = 'value') %>%
  mutate(correlation = factor(correlation),
         measure = recode(measure,
                          optcost = 'Optimal\ncost',
                          optM = 'Optimal\nnumber of sites',
                          optcatchpersite = 'Optimal\ncatch per site',
                          opttotalcatch = 'Optimal\ntotal catch') %>%
           fct_relevel('Optimal\nnumber of sites')) %>%
  subset(measure %in% c('Optimal\nnumber of sites', 'Optimal\ncatch per site') &
           theta >= 0.001) %>%
  ggplot(aes(x = theta, y = value, color = correlation)) +
  facet_grid(measure~cost.pool+cost.location,switch = 'y',scales = 'free',
             labeller = label_bquote(rows = .(measure),
                                     cols = c[u] * {phantom() == phantom()} * 1 * ';' *
                                       c[p] * {phantom() == phantom()} * .(cost.pool) * ';' *
                                       c[s] * {phantom() == phantom()} * .(cost.location))) +
  geom_smooth(se = FALSE) +
  #geom_point() +
  lims(y = c(0,NA)) +
  scale_x_continuous(labels = scales::percent, limits = c(0, NA),
                     breaks = c(0,0.002,0.004,0.006,0.008,0.01)) +
  scale_y_log10()+
  labs(x = 'Prevalence',
       y = '',
       color = 'ICC') +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        legend.position = 'bottom', text = element_text(size = 12)) +
  guides(x = guide_axis(angle = -90))
odp

ggsave('./Figures/OptimiseDetectionPower.png', odp,
       width = 4, height = 4, unit = 'in',dpi = 'retina')


cost.unit <- 1
cost.pool <- 1
cost.location <- 1
s <- 10; NM <- 35; M <- 5; detection_errors(theta = 0.01, group_size = s, group_num = NM/M,location_num = M, sensitivity = 0.9, specificity = 0.99, sigma = 1);     cost.unit * s*NM + cost.pool * NM + cost.location * M





optimise_detection_fixedsites <- function(theta, max.typeI, max.typeII, sensitivity, specificity, sigma,
                                          cost.unit, cost.pool, M){
  #if we assume that typeI error is 1 - specificity^(N*M) then this can be
  #translated into a simple constraint on N and M
  maxtests <- if(specificity <1){log1p(- max.typeI)/log(specificity)}else{Inf}

  if(maxtests < M){
    stop('It is not possible to stay below the specified type I error with the given test specificity without pooling specimens collected across multiple sites and this tool does not consider this case.')
  }

  #Since a hierarchical survey is less likely to find a positive than a simple random sample of the same size, we can use this to calculate a minimum number of units that need to be tested
  minunits <- log(max.typeII)/log((1-theta) * specificity + theta * (1 - sensitivity))

  cost <- function(pars){
    s <- pars[1];
    N <- pars[2];
    cost.unit * s*N*M + cost.pool * N*M
  }

  errors <- function(pars){
    s <- pars[1];
    N <- pars[2];
    er <- unlist(detection_errors(p = theta,
                                  group_size = s,
                                  group_num = N,
                                  location_num = M,
                                  sensitivity = sensitivity,
                                  specificity = specificity,
                                  sigma = sigma))
    er
  }

  init <- if(maxtests != Inf){c(minunits/maxtests, maxtests/M)}else{c(10, 10)}

  out <- Rsolnp::solnp(init, cost, ineqfun = errors,
                       ineqUB = c(max.typeI, max.typeII),
                       ineqLB = c(0,0),
                       LB = c(1,5))
  out
}


optimise_detection_fixedsites(theta = 0.001,
                              max.typeI =  0.3,
                              max.typeII = 0.05,
                              sensitivity = 0.9,
                              specificity = 0.999,
                              cost.unit = 1,
                              cost.pool = 1,
                              sigma = 0,
                              M = 10)






theta <- 0.1
sensitivity <- 0.99
specificity <- 1
sigma <- 10
cost.unit <- 1
cost.pool <- 5
cost.location <- 100
ICC <-  sigma^2/(sigma^2 + pi^2/3); ICC


a <- optimise_detection(theta = theta,
                        max.typeI = 0.05, max.typeII = 0.01,
                        sensitivity = sensitivity, specificity = specificity,
                        sigma = sigma,
                        cost.unit = cost.unit, cost.pool = cost.pool, cost.location = cost.location,
                        int.step = 30, int.search.only = F,
                        LB = c(0,10,30),
                        UB = c(Inf, 10, 30))
a
detection_errors(theta,
                 s = 10,sigma = sigma, sensitivity = sensitivity, specificity = specificity)

optimise_detection <- function(theta, max.typeI, max.typeII, sensitivity, specificity, sigma,
                               cost.unit, cost.pool, cost.location, init = NULL,
                               LB = c(1,1,1), UB = c(Inf,Inf,Inf), int.step = 3, int.search.only = F){


  #if we assume that typeI error is 1 - specificity^(N*M) then this can be
  #translated into a simple constraint on N and M
  maxtests <- if(specificity == 1){Inf}else{log(1 - max.typeI)/log(specificity)}

  #Since a hierarchical survey is less likely to find a positive than a simple random sample of the same size, we can use this to calculate a minimum number of units that need to be tested
  minunits <- log(max.typeII)/log((1-theta) * specificity + theta * (1 - sensitivity))

  if(minunits/UB[1] > maxtests){
    print(maxtests)
    print(minunits)
    stop('It is not possible to achieve the desired typeI and typeII error rates even with a simple random sample and using a heirarchical/cluster sampling scheme will perform even more poorly')
  }

  cost <- function(par){
    s <- par[1]; N <- par[2]; M <- par[3]
    cost.unit * s*N*M + cost.pool * N*M + cost.location * M
  }

  cost_totals <- function(par){
    cost.unit * par[1] + cost.pool * par[2] + cost.location * par[3]
  }

  errors <- function(par){
    unlist(detection_errors(theta,par[1],par[2],par[3],sensitivity,specificity,sigma))
  }

  real_constraints <- function(par){
    s <- par[1]/par[2]
    N <- par[2]/par[3]
    M <- par[3]

    c(unlist(detection_errors(theta,s,N,M,sensitivity,specificity,sigma)),
      s,N,M)
  }

  if(is.null(init)){
    ICC <- sigma^2/(sigma^2 + pi^2/3)
    init.M <- minunits * ICC + 1
    init.units <- (1 - ICC)/ (1/minunits -  ICC / init.M)
    init <- c(init.units, min(maxtests-1, init.units), init.M)
  }
  if(int.search.only){
    opt <- list(pars = c(init[1]/init[2],init[2]/init[3],init[3]))
  }else{
    real_opt <- Rsolnp::solnp(init, cost_totals, ineqfun = real_constraints,
                              ineqLB = c(0,0,LB),
                              ineqUB = c(max.typeI, max.typeII,UB),
                              LB = c(0,0,0),
                              UB = c(Inf,Inf,Inf),
                              control = list(delta = 1.0e-2,
                                             tol = 1.0e-3)
    )
    opt <- real_opt
    opt <- list(pars = c(opt$pars[1]/opt$pars[2],opt$pars[2]/opt$pars[3],opt$pars[3]))
  }
  n <- 0
  converged <- FALSE
  while(!converged){
    n <- n + 1
    print(paste("Local integer search:", s, ". Pars: ",paste(round(opt$pars),collapse = " ")))
    opt_new <- opt_local_int(opt$pars,cost,ineqfun = errors,
                             ineqLB = c(0,0),
                             ineqUB = c(max.typeI, max.typeII),
                             LB = LB,
                             UB = UB,
                             steps = int.step)
    if(all(opt_new$pars == opt$pars)) converged = TRUE
    opt <- opt_new
  }
  names(opt) <- c("design", "cost")
  names(opt$design) <- c("s","N","M")
  names(opt$cost) <- NULL
  opt$errors <- errors(opt$design)
  names(opt$errors) <- c("typeI","typeII")
  opt
}


# if(maxtests>minunits){
#   warn_message <- 'It is not possible to achieve the specified type I and and type II error with the given design prevalence, test specificity, and prevalence.
#   You will need to either:
#   * increase the allowable type I rate,
#   * increase the allowable type II error rates
#   * increase the design prevalence'
#   if(sensitivity != 1){
#     warn_message <- paste0(warn_message, '
#                            * use a test with better sensitivity')
#   }
#   if(specificity != 1){
#     warn_message <- paste0(warn_message, '
#                            * use a test with better specificity')
#   }
#   warn_message <- paste0(warn_message,'
#                          Changing the costs will not help.')
#   warning(warn_message)}
