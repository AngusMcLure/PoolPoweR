design_effect_cluster_fisher <- function(s,N,prevalence,sensitivity,specificity, correlation, form = 'beta'){
   N * s * fi_pool_imperfect(1,prevalence,sensitivity,specificity) *
      solve(fi_pool_imperfect_cluster(s,N,prevalence,sensitivity,specificity, correlation,form))[1,1]
}


#optimising the pool size for estimating prevalence

unit_fi_cost <- function(s,prevalence,sensitivity,specificity, cost.unit, cost.test){
  (cost.unit * s + cost.test)/fi_pool_imperfect(s,prevalence,sensitivity,specificity)
}

unit_fi_cost_clustered <- function(s,N,prevalence,correlation,sensitivity,specificity, cost.unit, cost.test, cost.location,form = 'beta'){
  fi <- fi_pool_imperfect_cluster(s,N,prevalence,sensitivity,specificity,correlation, form = form)
  cost <- cost.unit * s * N + cost.test * N + cost.location
  #print(c(cost = cost, information = fi))
  cost * solve(fi)[1,1]
}

# unit_fi_cost_clustered_alt <- function(s,N,prevalence,correlation,sensitivity,specificity,
#                                       cost.unit, cost.test, cost.location, cost.collect, catch.collect, form = 'beta'){
#
#   fi <- fi_pool_imperfect_cluster(s,N,prevalence,sensitivity,specificity,correlation, form = form)
#
#   ##   expected number of sampling periods required to get N pools of size s.
#   ##   tmight not be an integer!!
#   #t <- N * s/catch.collect
#   #cost <- cost.unit * s * N +  cost.test * N + cost.location + cost.collect * t
#
#   cost <- (cost.collect/catch.collect + cost.unit) * s * N + cost.test * N + cost.location + cost.collect * t
#   cost * solve(fi)[1,1]
# }


optimise_s_prevalence <- function(prevalence, cost.unit, cost.test,
                                  cost.location = 0, correlation = 0,
                                  N = 1, form = 'beta',
                                  sensitivity = 1,specificity = 1,
                                  max.s = 50, interval = 0){
  invalid.cost <- FALSE #trigger for when costs are infinite to ensure that there's no cost output in these cases
  #print(c(theta = prevalence, sens = sensitivity, spec = specificity, unit = cost.unit, test = cost.test, location = cost.location , rho = correlation, N = N, form = form, max.s = max.s))
  theta <- prevalence
  ## The case cost.test == Inf (or equivalently cost.unit = 0) is helpful because
  ## the s in this case (s_opt) is the largest you would ever want
  ## to make a pool. If you can only do N tests and have more than N*s_opt
  ## vectors you should just do N tests of size s_opt and not test the remaining
  ## vectors. In practice this upper limit is usually really high (especially if
  ## prevalence is low and/or specificity is not perfect), so you can assume in
  ## most designs that you want to test all vectors
  if(cost.test == Inf){
    cost.unit <- 0
    cost.test <- 1
    invalid.cost <- TRUE
  }

  ## This case is also very helpful. The solution is trivially s_opt = 1 for
  ## perfect specificity, however it can get very large for low prevalence and
  ## low specificity. This gives you the case where you have decided to not
  ## collect any more vectors and are now just trying to decide how to the pool
  ## them. This gives you the lower bound for the pool size, i.e. even if you
  ## have budget to test them all individually you're actually better off
  ## testing them in pools of this size
  if(cost.unit == Inf){
    cost.test <- 0
    cost.unit <- 1
    invalid.cost <- TRUE

  }
  if(correlation == 0){
    cost <- function(x){
      ufc <- unit_fi_cost(x,theta,sensitivity,specificity, cost.unit, cost.test)
      ufc}
  }else{
    cost <- function(x){
      ufc <- unit_fi_cost_clustered(s = x, N = N, prevalence = theta,
                                    correlation = correlation,
                                    sensitivity = sensitivity,
                                    specificity = specificity,
                                    cost.unit = cost.unit, cost.test = cost.test,
                                    cost.location = cost.location, form = form)
      ufc}
  }

  opt <- optimise(cost,c(1,max.s))
  cost_opt_ceiling <- cost(ceiling(opt$minimum))
  cost_opt_floor   <- cost(floor(opt$minimum))

  if(cost_opt_ceiling < cost_opt_floor){
    s <- ceiling(opt$minimum)
    opt_cost <- cost_opt_ceiling
  }else{
    s <- floor(opt$minimum)
    opt_cost <- cost_opt_floor
  }

  if(s == max.s) warning('Maximum cost effectivness is achieved at or above the maximum size of pools allowed. Consider increasing max.s')

  if(interval == 0){
    if(invalid.cost) opt_cost <- NA
    return(list(s = s, cost = opt_cost, catch = s * N))
  }else if(interval > 0){
    max_cost <- opt_cost * (1 + interval)
    if(cost(1) < max_cost){
      lower <- 1
    } else{
      lower <- ceiling(uniroot(function(x){cost(x) - max_cost},c(1,s))$root)
    }
    if(cost(max.s) < max_cost){
      upper <- max.s
      warning('A pool size greater than max.s may fall within the specified range of cost effectiveness. Consider increasing max.s')
    }else{
      upper <- floor(uniroot(function(x){cost(x) - max_cost},c(s,max.s))$root)
    }
    if(invalid.cost){
      opt_cost <- NA
      cost_interval <- NA}
    else{
      cost_interval <- c(cost(lower),cost(upper))
    }
    out <- list(s = s,
                cost = opt_cost,
                catch = s * N,
                s_interval = c(lower, upper),
                cost_interval = cost_interval,
                catch_interval = N * c(lower, upper))
    return(out)
  }else{
    stop('interval must be 0 (for no interval) or greater than zero.')
  }
}

optimise_sN_prevalence <- function(prevalence, cost.unit, cost.test,
                                   cost.location, correlation, form = 'beta',
                                   sensitivity = 1, specificity = 1,
                                   max.s = 50, max.N = 20){
  #print(c(theta = prevalence, sens = sensitivity, spec = specificity, unit = cost.unit, test = cost.test, location = cost.location , rho = correlation, N = N, form = form, max.s = max.s))
  theta <- prevalence

  if(correlation == 0){
    warning('If there is no correlation between units at locations, then this means sampling at a single random location is a identical to sampling from the whole population. This is assumption unlikely to be true in most settings, but would mean that sampling at a single site is the most cost-effective strategy')
    opt <- c(optimise_s_prevalence(prevalence, cost.unit, cost.test,
                                   cost.location, correlation,
                                   N = 1, form = 'beta',
                                   sensitivity,specificity,
                                   max.s = max.s),
             list(N = Inf))

  }else{
    Nopt <- 1
    opt <- list(cost = Inf)
    while(Nopt < max.N){
      opt_new <- optimise_s_prevalence(prevalence,
                                       cost.unit, cost.test, cost.location,
                                       correlation, Nopt+1, form,
                                       sensitivity,specificity, max.s)
      #print(opt_new)
      if(opt_new$cost > opt$cost){
        break
      }else{
        Nopt <- Nopt + 1
        opt <- opt_new
      }
    }
    opt$N <- Nopt
  }
  opt
}

optimise_N_prevalence <- function(prevalence, cost.unit, cost.test,
                                   cost.location, cost.collect, catch.collect,
                                   correlation, P, form = 'beta',
                                   sensitivity = 1, specificity = 1,
                                   max.s = 50, max.N = 20){
  # For given set up including number of sampling periods and catch per period (i.e. given samples per location)
  # calculate the number of pools that is optimal. Resulting pool size will usually not be an integer
  n <- P * catch.collect
  N.opt <- max(if(correlation == 0){1}else{2}, ceiling(n/max.s)) - 1
  cost.opt <- Inf
  while(N.opt<n){
    cost.new <- unit_fi_cost_clustered(n/(N.opt+1),N.opt+1,prevalence,correlation,
                                       sensitivity,specificity,
                                       cost.unit + cost.collect/catch.collect,
                                       cost.test,cost.location,form)
    if(cost.new > cost.opt){
      break
    }
    else{
      cost.opt <- cost.new
      N.opt <- N.opt + 1
    }
  }
  out <- list(N = N.opt, cost = cost.opt, s = n/N.opt)
  out
}

optimise_NP_prevalence <- function(prevalence, cost.unit, cost.test,
                                   cost.location, cost.collect, catch.collect,
                                   correlation, form = 'beta',
                                   sensitivity = 1, specificity = 1,
                                   max.s = 50, max.P = 20){
  #print(c(theta = prevalence, sens = sensitivity, spec = specificity, unit = cost.unit, test = cost.test, location = cost.location , rho = correlation, N = N, form = form, max.s = max.s))
  theta <- prevalence

  if(correlation == 0){
    stop('If there is no correlation between units at locations, then this means sampling at a single random location is a identical to sampling from the whole population. This is assumption unlikely to be true in most settings, but would mean that sampling at a single site is the most cost-effective strategy')
    # opt <- c(optimise_s_prevalence(prevalence, cost.unit + cost.collect/catch.collect,
    #                                cost.test, cost.location, correlation,
    #                                N = 1, form = 'beta',
    #                                sensitivity,specificity,
    #                                max.s = max.s),
    #          list(P = Inf))

  }else{
    P.opt <- 0
    opt <- list(cost = Inf)
    while(P.opt < max.P){
      opt.new <- optimise_N_prevalence(prevalence, cost.unit, cost.test,
                                       cost.location, cost.collect, catch.collect,
                                       correlation, P.opt+1, form = 'beta',
                                       sensitivity, specificity,
                                       max.s, max.N)

      if(opt.new$cost > opt$cost){
        break
      }else{
        P.opt <- P.opt + 1
        opt <- opt.new
      }
    }
    opt$P <- P.opt
  }
  opt$catch <- opt$P * catch.collect
  opt$s <- opt$P * catch.collect/opt$N
  opt
}








