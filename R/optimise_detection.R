### NEED TO CHECK TO SEE WHETHER n HERE MEANS TOTAL SAMPLE SIZE OR POOL SIZE AND CORRECT AS NECESSARY


optimise_detection_fixedn <- function(n,design.prevalence, max.typeI, max.typeII, sensitivity, specificity, sigma,
                                          cost.unit, cost.pool, cost.location){
  #if we assume that typeI error is 1 - specificity^(N*M) then this can be
  #translated into a simple constraint on N and M
  maxtests <- if(specificity <1){log1p(- max.typeI)/log(specificity)}else{Inf}

  #Since a hierarchical survey is less likely to find a positive than a simple random sample of the same size, we can use this to calculate a minimum number of units that need to be tested
  minunits <- log(max.typeII)/log((1-design.prevalence) * specificity + design.prevalence * (1 - sensitivity))

  if(minunits/n > maxtests){
    stop('It is not possible to achieve both the specified type I error and type II error with the test. You must either:
         * increase the design prevlance
         * use larger poolsize (n)
         * use a test with higher sensitivity and/or specificity
         * or relax (increase) your maximum error requirements.
         minunits = ', minunits,'
         maxtests = ', maxtests)
  }


  cost <- function(pars){
    NM <- pars[1];
    M <- 1/pars[2];
    cost.unit * n*NM + cost.pool * NM + cost.location * M
  }

  errors <- function(pars){
    NM <- pars[1];
    M <- 1/pars[2];
    er <- unlist(detection_errors(p = design.prevalence,
                                  group_size = n,
                                  group_num = NM/M,
                                  location_num = M,
                                  sensitivity = sensitivity,
                                  specificity = specificity,
                                  sigma = sigma))
    #print(pars)
    ##print(er)
    er
  }

  init <- c(minunits/n * 1.5, 1/10)
  print(init)
  print(errors(init))

  out <- Rsolnp::solnp(init, cost, ineqfun = errors,
                       ineqUB = c(max.typeI, max.typeII),
                       ineqLB = c(-Inf,-Inf),
                       LB = c(1,0))
  #print(purrr::map((round(out$pars)-10):(round(out$pars)+10), ~c(cost(.x), errors(.x))))
  out
}

optimise_detection_fixedn(n = 10,
                          design.prevalence = 0.01,
                          max.typeI =  0.3,
                          max.typeII = 0.05,
                          sensitivity = 0.9,
                          specificity = 0.99,
                          cost.unit = 1,
                          cost.pool = 1,
                          cost.location = 1,
                          sigma = 1)

cost.unit <- 1
cost.pool <- 1
cost.location <- 1
n <- 10; NM <- 35; M <- 5; detection_errors(p = 0.01, group_size = n, group_num = NM/M,location_num = M, sensitivity = 0.9, specificity = 0.99, sigma = 1);     cost.unit * n*NM + cost.pool * NM + cost.location * M





optimise_detection_fixedsites <- function(design.prevalence, max.typeI, max.typeII, sensitivity, specificity, sigma,
                                          cost.unit, cost.pool, M){
  #if we assume that typeI error is 1 - specificity^(N*M) then this can be
  #translated into a simple constraint on N and M
  maxtests <- if(specificity <1){log1p(- max.typeI)/log(specificity)}else{Inf}

  if(maxtests < M){
    stop('It is not possible to stay below the specified type I error with the given test specificity without pooling specimens collected across multiple sites and this tool does not consider this case.')
  }

  #Since a hierarchical survey is less likely to find a positive than a simple random sample of the same size, we can use this to calculate a minimum number of units that need to be tested
  minunits <- log(max.typeII)/log((1-design.prevalence) * specificity + design.prevalence * (1 - sensitivity))

  cost <- function(pars){
    n <- pars[1];
    N <- pars[2];
    cost.unit * n*N*M + cost.pool * N*M
  }

  errors <- function(pars){
    n <- pars[1];
    N <- pars[2];
    er <- unlist(detection_errors(p = design.prevalence,
                                  group_size = n,
                                  group_num = N,
                                  location_num = M,
                                  sensitivity = sensitivity,
                                  specificity = specificity,
                                  sigma = sigma))
    #print(pars)
    ##print(er)
    er
  }

  init <- if(maxtests != Inf){c(minunits/maxtests, maxtests/M)}else{c(10, 10)}

  out <- Rsolnp::solnp(init, cost, ineqfun = errors,
                       ineqUB = c(max.typeI, max.typeII),
                       ineqLB = c(0,0),
                       LB = c(1,5))
  #print(purrr::map((round(out$pars)-10):(round(out$pars)+10), ~c(cost(.x), errors(.x))))
  out
}


optimise_detection_fixedsites(design.prevalence = 0.001,
                              max.typeI =  0.3,
                              max.typeII = 0.05,
                              sensitivity = 0.9,
                              specificity = 0.999,
                              cost.unit = 1,
                              cost.pool = 1,
                              sigma = 0,
                              M = 10)




optimise_detection <- function(design.prevalence, max.typeI, max.typeII, sensitivity, specificity, sigma,
                               cost.unit, cost.pool, cost.location, init = NULL, init_M = 30,
                               LB = c(1,1,1), UB = c(Inf,Inf,Inf), int.step = 3, int.search.only = F){


  #if we assume that typeI error is 1 - specificity^(N*M) then this can be
  #translated into a simple constraint on N and M
  maxtests <- log(1 - max.typeI)/log(specificity)

  #Since a hierarchical survey is less likely to find a positive than a simple random sample of the same size, we can use this to calculate a minimum number of units that need to be tested
  minunits <- log(max.typeII)/log((1-design.prevalence) * specificity + design.prevalence * (1 - sensitivity))

  cost <- function(par){
    n <- par[1]; N <- par[2]; M <- par[3]
    cost.unit * n*N*M + cost.pool * N*M + cost.location * M
  }

  errors <- function(par){
    unlist(detection_errors(design.prevalence,par[1],par[2],par[3],sensitivity,specificity,sigma))
  }

  #initial guess
  if(is.null(init)){
    init_M <- max(LB[3] + 1, init_M)
    init_N <- max(LB[2] + 1, floor(maxtests/init_M))
    init_n <- max(LB[1] + 1, ceiling(minunits/(init_N * init_M)) * (sigma + 1) * 2)
    init <- c(init_n,init_N, init_M)
  }
  print(init)
  if(int.search.only){
    opt <- list(pars = init)
  }else{
    real_opt <- Rsolnp::solnp(init,cost, ineqfun = errors,
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
  names(opt$design) <- c("n","N","M")
  names(opt$cost) <- NULL
  opt$errors <- errors(opt$design)
  names(opt$errors) <- c("typeI","typeII")
  opt
}

opt_local_int <- function(pars,fun,ineqfun, ineqLB, ineqUB,LB, UB, steps){
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
  list(pars = opt_par, value = opt_val)
}

design.prevalence <- 0.1
sensitivity <- 0.99
specificity <- 1
sigma <- 10
cost.unit <- 1
cost.pool <- 5
cost.location <- 100
ICC <-  sigma^2/(sigma^2 + pi^2/3); ICC


a <- optimise_detection(design.prevalence = design.prevalence,
                        max.typeI = 0.05, max.typeII = 0.01,
                        sensitivity = sensitivity, specificity = specificity,
                        sigma = sigma,
                        cost.unit = cost.unit, cost.pool = cost.pool, cost.location = cost.location,
                        int.step = 30, int.search.only = F,
                        LB = c(0,10,30),
                        UB = c(Inf, 10, 30))
a
detection_errors(design.prevalence,
                 n = 10,sigma = sigma, sensitivity = sensitivity, specificity = specificity)

optimise_detection <- function(design.prevalence, max.typeI, max.typeII, sensitivity, specificity, sigma,
                               cost.unit, cost.pool, cost.location, init = NULL,
                               LB = c(1,1,1), UB = c(Inf,Inf,Inf), int.step = 3, int.search.only = F){


  #if we assume that typeI error is 1 - specificity^(N*M) then this can be
  #translated into a simple constraint on N and M
  maxtests <- if(specificity == 1){Inf}else{log(1 - max.typeI)/log(specificity)}

  #Since a hierarchical survey is less likely to find a positive than a simple random sample of the same size, we can use this to calculate a minimum number of units that need to be tested
  minunits <- log(max.typeII)/log((1-design.prevalence) * specificity + design.prevalence * (1 - sensitivity))

  if(minunits/UB[1] > maxtests){
    print(maxtests)
    print(minunits)
    stop('It is not possible to achieve the desired typeI and typeII error rates even with a simple random sample and using a heirarchical/cluster sampling scheme will perform even more poorly')
  }

  cost <- function(par){
    n <- par[1]; N <- par[2]; M <- par[3]
    cost.unit * n*N*M + cost.pool * N*M + cost.location * M
  }

  cost_totals <- function(par){
    cost.unit * par[1] + cost.pool * par[2] + cost.location * par[3]
  }

  errors <- function(par){
    unlist(detection_errors(design.prevalence,par[1],par[2],par[3],sensitivity,specificity,sigma))
  }

  real_constraints <- function(par){
    n <- par[1]/par[2]
    N <- par[2]/par[3]
    M <- par[3]

    c(unlist(detection_errors(design.prevalence,n,N,M,sensitivity,specificity,sigma)),
      n,N,M)
  }

  if(is.null(init)){
    ICC <- sigma^2/(sigma^2 + pi^2/3)
    init.M <- minunits * ICC + 1
    init.units <- (1 - ICC)/ (1/minunits -  ICC / init.M)
    init <- c(init.units, min(maxtests-1, init.units), init.M)
  }

  print(init)
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
  names(opt$design) <- c("n","N","M")
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
