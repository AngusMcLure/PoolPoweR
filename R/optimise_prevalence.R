#' Optimising the pool size and number for estimating prevalence.
#'
#' These functions determine cost-effective pooling strategies for estimating
#' the prevalence of a marker in a population. Both functions attempt to choose
#' survey designs that maximise the Fisher Information for given cost or effort.
#' `optimise_s_prevalence()` calculates the optimal single pool size that
#' balances the cost and accuracy given the marker prevalence, test sensitivity,
#' and specificity, and works for simple random surveys or cluster surveys.
#' `optimise_sN_prevalence()` also attempts to identify the optimal number of
#' pools per cluster (cluster-surveys only). `optimise_random_prevalence()`
#' works for designs where collection proceeds for a number of collection
#' periods, but the number of units collected (and therefore the pool sizes,
#' pool numbers, and total cost) is random. Given the mean and variance of the
#' catch sizes per collection period and a family of pooling strategies, it
#' calculates the optimal combination of number of sampling periods and pooling
#' strategy.
#'
#' @param pool_number numeric The number of pools per cluster. Must be a numeric
#'   value greater than or equal to 0.
#' @param catch_mean,catch_variance numeric The mean and variance of the number
#'   of units collected per collection period. Both must be greater than 0 and
#'   `catch_variance` must be greater than or equal to `catch_mean`. Number of
#'   units caught per period is assumed to follow a negative binomial
#'   distribution (or Poisson distribution if `catch_mean = catch_variance`)
#' @param pool_strat_family function A function that defines a family of rules
#'   for how a number of units will be divided into pools, e.g.
#'   `pool_max_size()` and `pool_target_number()`. The function must take
#'   positive integer valued parameters, and return a function that defines a
#'   pooling strategy
#' @param prevalence numeric The proportion of units that carry the marker of
#'   interest (i.e. true positive). Must be be a numeric value between 0 and 1,
#'   inclusive of both.
#' @param correlation numeric The correlation between test results within a
#'   single cluster (units in different clusters are assumed to be
#'   uncorrelated). Must be a numeric value between 0 and 1, inclusive of both.
#'   A value of 1 indicates that units within clusters are perfectly correlated
#'   (there are no differences units within a single cluster). A value of 0
#'   indicates that units within clusters are no more correlated than units in
#'   different clusters.
#' @param sensitivity numeric The probability that the test correctly identifies
#'   a true positive. Must be a numeric value between 0 and 1, inclusive of
#'   both. A value of 1 indicates that the test can perfectly identify all true
#'   positives.
#' @param specificity numeric The probability that the test correctly identifies
#'   a true negative. Must be a numeric value between 0 and 1, inclusive of
#'   both. A value of 1 indicates that the test can perfectly identify all true
#'   negatives.
#' @param cost_unit numeric The cost to process a single unit. Must be a numeric
#'   value greater than or equal to 0.
#' @param cost_pool numeric The cost to process a single pool. Must be a numeric
#'   value greater than or equal to 0.
#' @param cost_period numeric The cost per collection period (per collection per
#'   cluster if `correlation` is not `NA`). Most be a numeric value greater than
#'   or equal to 0.
#' @param cost_cluster numeric/NA The cost to process a cluster. Must be a
#'   numeric value greater than or equal to 0. For `optimise_s_prevalence()`,
#'   this can be ignored if correlation is NA.
#' @param max_s numeric The maximum number of units per pool (pool size).
#' @param max_period numeric The maximum number of collection periods (per
#'   cluster if `correlation` is not `NA`)
#' @param form string The distribution used to model the cluster-level
#'   prevalence and correlation of units within cluster. Select one of "beta",
#'   "logitnorm" or "cloglognorm". See details.
#' @param interval Range of near-optimal designs to consider. If interval == 0
#'   (the default) only returns optimal design. If interval > 0, function
#'   identifies range of designs with cost less than the optimal cost * (1 +
#'   interval).
#' @param verbose logical Should function indicate progress by printing each
#'   parameter set calculated to screen? If FALSE (default) only prints optimal
#'   parameters identified in each iteration
#'
#' @returns `optimise_s_prevalence` returns a list with the optimal pool size
#'   `s`, the unit cost of Fisher information, and range of near-optimal designs
#'   `catch`. `optimise_sN_prevalence()` returns the same list as
#'   `optimise_s_prevalence` with an additional optimal pool number `N`.
#'   `optimise_random_prevlance` returns a list including the optimal number of
#'   sampling periods; the unit cost of Fisher information; the mean, variance,
#'   and distribution of the number of units caught across all sampling
#'   periods(per cluster where applicable); and the optimal pooling strategy.
#' @export
#'
#' @examples
#' optimise_s_prevalence(prevalence = 0.01, cost_unit = 5, cost_pool = 10)
#' 
#' optimise_prevalence(
#'   fixed_design(), 
#'   prevalence = 0.01, cost_unit = 5, cost_pool = 10,
#'   cost_cluster = 100, correlation = 0.05, form = "beta"
#' )
optimise_s_prevalence <- function(pool_number = 1,
                                  prevalence,
                                  cost_unit,
                                  cost_pool,
                                  cost_cluster = NA,
                                  correlation = NA,
                                  sensitivity = 1,
                                  specificity = 1,
                                  max_s = 50,
                                  form = "logitnorm",
                                  interval = 0) {
  N <- pool_number

  # Input checks ----
  check_geq2(cost_unit, 0) 
  check_geq2(cost_pool, 0) 
  check_geq2(interval, 0) 
  check_geq2(max_s, 1, allow_na = TRUE)
  check_geq2(cost_cluster, 0, allow_na = TRUE)
  check_in_range2(prevalence)
  check_in_range2(correlation, allow_na = TRUE)

  if (form == "discrete") {
    stop('When form = "discrete" the cost of unit information function with
         respect to s often has mulitple minima and therefore the discrete
         distribution is not currently supported for optimisation')
  }

  # Cost handling ----
  invalid_cost <- FALSE # ensure that there's no cost output when costs = Inf

  ## The case cost_pool == Inf (or equivalently cost_unit = 0) is helpful because
  ## the s in this case (s_opt) is the largest you would ever want
  ## to make a pool. If you can only do N tests and have more than N*s_opt
  ## vectors you should just do N tests of size s_opt and not test the remaining
  ## vectors. In practice this upper limit is usually really high (especially if
  ## prevalence is low and/or specificity is not perfect), so you can assume in
  ## most designs that you want to test all vectors
  if (cost_pool == Inf) {
    cost_unit <- 0
    cost_pool <- 1
    invalid_cost <- TRUE
  }

  ## This case is also very helpful. The solution is trivially s_opt = 1 for
  ## perfect specificity, however it can get very large for low prevalence and
  ## low specificity. This gives you the case where you have decided to not
  ## collect any more vectors and are now just trying to decide how to the pool
  ## them. This gives you the lower bound for the pool size, i.e. even if you
  ## have budget to test them all individually you're actually better off
  ## testing them in pools of this size
  if (cost_unit == Inf) {
    cost_pool <- 0
    cost_unit <- 1
    invalid_cost <- TRUE
  }

  # Assign clustered/unclustered function ----
  if (is.na(correlation)) {
    cost <- function(s) {
      cost_fi(s, prevalence, sensitivity, specificity, cost_unit, cost_pool)
    }
  } else {
    cost <- function(s) {
      cost_fi_cluster(
        pool_size = s, 
        pool_number = N, 
        prevalence = prevalence,
        correlation = correlation,
        sensitivity = sensitivity,
        specificity = specificity,
        cost_unit = cost_unit, 
        cost_pool = cost_pool,
        cost_cluster = cost_cluster, 
        form = form
      )
    }
  }

  # Optimise cost and interval ----
  opt <- stats::optimise(cost, c(1, max_s))
  cost_opt_ceiling <- cost(ceiling(opt$minimum))
  cost_opt_floor <- cost(floor(opt$minimum))

  if (cost_opt_ceiling < cost_opt_floor) {
    s <- ceiling(opt$minimum)
    opt_cost <- cost_opt_ceiling
  } else {
    s <- floor(opt$minimum)
    opt_cost <- cost_opt_floor
  }

  if (s == max_s) {
      warning("Maximum cost effectivness is achieved at or above the maximum size of pools allowed. Consider increasing max_s")
  }

  if (interval == 0) {
    if (invalid_cost) opt_cost <- NA
    return(list(s = s, cost = opt_cost, catch = s * N))
  } else if (interval > 0) {
    max_cost <- opt_cost * (1 + interval)
    if (cost(1) < max_cost) {
      lower <- 1
    } else {
      lower <- ceiling(stats::uniroot(function(x) {
        cost(x) - max_cost
      }, c(1, s))$root)
    }
    if (cost(max_s) < max_cost) {
      upper <- max_s
      warning("A pool size greater than max_s may fall within the specified range of cost effectiveness. Consider increasing max_s")
    } else {
      upper <- floor(stats::uniroot(function(x) {
        cost(x) - max_cost
      }, c(s, max_s))$root)
    }
    if (invalid_cost) {
      opt_cost <- NA
      cost_interval <- NA
    } else {
      cost_interval <- c(cost(lower), cost(upper))
    }
    out <- list(
      s = s,
      cost = opt_cost,
      catch = s * N,
      s_interval = c(lower, upper),
      cost_interval = cost_interval,
      catch_interval = N * c(lower, upper)
    )
  }
  out
}

#' Optimise prevalence
#' 
#' Placeholder documentation to track params for the refactor.
#' 
#' @inheritParams design_effect
#' @param cost_unit numeric The cost to process a single unit. Must be a numeric
#'   value greater than or equal to 0.
#' @param cost_pool numeric The cost to process a single pool. Must be a numeric
#'   value greater than or equal to 0.
#' @param cost_cluster numeric/NA The cost to process a cluster. Must be a
#'   numeric value greater than or equal to 0. For `optimise_s_prevalence()`,
#'   this can be ignored if correlation is NA.
#' @param max_s numeric The maximum number of units per pool (pool size).
#' @param ... Additional arguments passed to methods.
#' @export
optimise_prevalence <- function(x, ...) {
  UseMethod("optimise_prevalence")
}

#' @rdname optimise_prevalence
#' 
#' @inheritParams optimise_prevalence
#' @param max_N numeric The maximum number of pools per cluster (pool number).
#' @method optimise_prevalence fixed_sN
#' @export
optimise_prevalence.fixed_sN <- function(x, 
                                         prevalence,
                                         cost_unit,
                                         cost_pool,
                                         cost_cluster,
                                         correlation,
                                         max_s = 50,
                                         max_N = 20,
                                         form = "logitnorm",
                                         ...) {

  # max_N is the only argument that optimise_s_prevalence doesn't use
  # The rest of the input_checks are in optimise_s_prevalence
  check_geq2(max_N, 1)

  # Simple random sampling ----
  ## for simple random sampling there is no optimal N and for no correlation
  ## cluster survey the optimal approach is to infinite pools at one site
  ## (i.e. a simple random survey)
  if (is.na(correlation) || correlation == 0) { 
    opt <- optimise_s_prevalence(
      pool_number = 1, prevalence, cost_unit, cost_pool, cost_cluster,
      correlation = NA, x$sensitivity, x$specificity, max_s, form
    )

    na_or_inf <- ifelse(is.na(correlation), NA, Inf)
    
    out <- fixed_design(
      pool_size = opt$s,
      pool_number = na_or_inf,
      total_units = na_or_inf,
      sensitivity = x$sensitivity,
      specificity = x$specificity
    )
    
    return(out)
  }

  # Iterate to find optimal N ----
  N <- 1
  opt <- list(cost = Inf)
  while (N < max_N) {
    opt_new <- optimise_s_prevalence(
      pool_number = N + 1, prevalence, cost_unit, cost_pool, cost_cluster,
      correlation, x$sensitivity, x$specificity, max_s, form
    )
    if (opt_new$cost > opt$cost) {
      break
    } else {
      N <- N + 1
      opt <- opt_new
    }
  }

  # Warnings ----
  opt$N <- N
  if (opt$N == max_N) {
    warning("Maximum cost effectivness is achieved at or above the maximum number of pools allowed. Consider increasing max_N")
  }
  if (opt$s == max_s) {
    warning("Maximum cost effectivness is achieved at or above the maximum size of pools allowed. Consider increasing max_s")
  }

  # Output optimal results ----
  fixed_design(
    pool_size = opt$s,
    pool_number = opt$N,
    sensitivity = x$sensitivity,
    specificity = x$specificity
  )
}

#' @rdname optimise_s_prevalence
#' @export
optimise_random_prevalence <- function(catch_mean, catch_variance,
                                       pool_strat_family,
                                       prevalence,
                                       cost_unit, cost_pool, cost_period, cost_cluster = NA,
                                       correlation = NA,
                                       sensitivity = 1, specificity = 1,
                                       max_period = 10, form = "logitnorm",
                                       verbose = FALSE) {
  
  strat_par_names <- names(formals(pool_strat_family))
  
  if(is.na(correlation)){ #If correlation is NA, assumes no correlation. However, sampling everything from the one site is optimal (i.e. periods = Inf). Consequently, we assume period = 1 
    
    pars <- as.list(rep(1, length(strat_par_names)))
    names(pars) <- strat_par_names
    catch <- nb_catch(catch_mean, catch_variance)
    f <- function(prm){
      cost_fi_random(catch,
                     do.call(pool_strat_family, prm),
                     prevalence, sensitivity,specificity,
                     cost_unit, cost_pool, cost_period)
    }
    
    optpars <- optim_local_int(pars,f,verbose = verbose, max_iter = 100)
    
    
    opt <- list(periods = NA,
                cost = optpars$val,
                catch = list(mean = catch_mean,
                             variance = catch_variance,
                             distribution = catch),
                pool_strat = do.call(pool_strat_family, optpars$par),
                pool_strat_pars = optpars$par)
    
    
  }else{
    
    pars <- as.list(rep(1, 1 + length(strat_par_names)))
    names(pars) <- c(".periods", strat_par_names)
    
    f <- function(prm){
      cost_fi_cluster_random(nb_catch(catch_mean * prm$.periods, catch_variance * prm$.periods),
                             do.call(pool_strat_family, prm[strat_par_names]),
                             prevalence, correlation,
                             sensitivity,specificity,
                             cost_unit, cost_pool,cost_cluster + prm$.periods * cost_period,
                             form = form)
    }
    
    optpars <- optim_local_int(pars,f, verbose = verbose)
    
    
    opt <- list(periods = optpars$par$.periods,
                cost = optpars$val,
                catch = list(mean = optpars$par$.periods * catch_mean,
                             variance = optpars$par$.periods * catch_variance,
                             distribution = nb_catch(optpars$par$.periods * catch_mean,
                                                     optpars$par$.periods * catch_variance)),
                pool_strat = do.call(pool_strat_family, optpars$par[-1]),
                pool_strat_pars = optpars$par[-1])
    
  }

  opt
}


optim_local_int <- function(par, fn,
                            width = rep(1, length(par)),
                            par_lower = rep(1, length(par)),
                            par_upper = rep(Inf, length(par)),
                            max_iter = 20, verbose = FALSE){
  
  opt_par <- lapply(par, round)
  searched_par <- as.data.frame(opt_par)
  opt_val <- fn(opt_par)
  
  for(iter in 1:max_iter){
    print(paste0('Iteration ', iter, '. optimal parameters = ', paste(opt_par, collapse = ', '), ', with function value = ', opt_val))
    local_par_ranges <- purrr::pmap(list(op = opt_par, w = width, pl = par_lower, pu = par_upper), 
                                   \(op, w, pl, pu) max((op - w), pl):min((op + w), pu))
    local_par <- do.call(expand.grid,local_par_ranges)
    search <- dplyr::setdiff(local_par, searched_par)
    searched_par <- rbind(searched_par, search)
    
    val <- c()
    for(i in 1:nrow(search)){
      if(verbose){print(paste0('Calculating for paramter set #',i,': ', paste(unlist(search[i,]), collapse = ", ") ))}
      val[i] <- fn(as.list(search[i,,drop = FALSE]))
    }
    
    min_val <- min(val)
    if(opt_val <= min_val){
      break
    }else{
      opt_val <- min_val
      opt_par <- unlist(search[which.min(val),,drop=FALSE])
    }
  }
  if(iter == max_iter){warning('Local integer search reached max_iter and may have terminated early. Consider increasing max_iter')}
  return(list(val = opt_val, par = as.list(opt_par)))
}

cost_fi <- function(pool_size, prevalence, sensitivity, specificity, cost_unit, cost_pool) {
  (sum(cost_unit * pool_size) + cost_pool) /
    fi_pool(pool_size, prevalence, sensitivity, specificity)
}

cost_fi_cluster <- function(pool_size, pool_number,
                            prevalence, correlation,
                            sensitivity, specificity,
                            cost_unit, cost_pool, cost_cluster,
                            form = "logitnorm") {
  s <- pool_size
  N <- pool_number

  fi <- fi_pool_cluster(s, N, prevalence, correlation, sensitivity, specificity, form)
  cost <- sum(cost_unit * s * N) + sum(cost_pool * N) + cost_cluster
  # print(c(cost = cost, information = fi))
  cost * solve(fi)[1, 1]
}

cost_fi_random <- function(catch_dist, pool_strat, prevalence,
                           sensitivity, specificity, cost_unit, cost_pool, cost_fixed){
  
  fn <- function(catch){
    pooling <- pool_strat(catch)
    cost_unit * catch + cost_pool * sum(pooling$pool_number) + cost_fixed
  }
  
  fi <- fi_pool_random(catch_dist,pool_strat,prevalence,sensitivity, specificity)
  mean_cost <- ev(fn, catch_dist)
  mean_cost/fi
}

cost_fi_cluster_random <- function(catch_dist, pool_strat,
                                   prevalence, correlation,
                                   sensitivity, specificity,
                                   cost_unit, cost_pool, cost_cluster,
                                   form = "logitnorm"){

  fn <- function(catch){
    pooling <- pool_strat(catch)
    cost_unit * catch + cost_pool * sum(pooling$pool_number) + cost_cluster
  }
  
  fi <- fi_pool_cluster_random(catch_dist,pool_strat,
                               prevalence, correlation,
                               sensitivity,specificity,
                               form)
  mean_cost <- ev(fn, catch_dist)
  mean_cost * solve(fi)[1,1]
}
