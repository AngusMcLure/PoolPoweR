#' Optimising the pool size and number for estimating prevalence.
#'
#' These functions determine cost-effective pooling strategies for estimating
#' the prevalence of a marker in a population. Both functions attempt to choose
#' survey designs that maximise the Fisher Information for given cost or effort.
#' `optimise_s_prevalence()` calculates the optimal single pool size that
#' balances the cost and accuracy given the marker prevalence, test sensitivity,
#' and specificity, and works for simple random surveys or cluster surveys.
#' `optimise_sN_prevalence` also attempts to identify the optimal number of
#' pools per cluster (cluster-surveys only).
#'
#' @param pool_number numeric The number of pools per cluster. Must be a numeric
#'   value greater than or equal to 0.
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
#' @param cost_cluster numeric/NA The cost to process a cluster. Must be a
#'   numeric value greater than or equal to 0. For `optimise_s_prevalence()`,
#'   this can be ignored if correlation is NA.
#' @param max_s numeric The maximum number of units per pool (pool size).
#' @param max_N numeric The maximum number of pools per cluster (pool number).
#' @param form string The distribution used to model the cluster-level
#'   prevalence and correlation of units within cluster. Select one of "beta",
#'   "logitnorm" or "cloglognorm". See details.
#' @param interval Range of near-optimal designs to consider. If interval == 0
#'   (the default) only returns optimal design. If interval > 0, function
#'   identifies range of designs with cost less than the optimal cost * (1 +
#'   interval).
#'
#' @returns * `optimise_s_prevalence` returns a list with the optimal pool size
#'   `s`, cost, and range of near-optimal designs `catch`. *
#'   `optimise_sN_prevalence()` returns the same list as `optimise_s_prevalence`
#'   with an additional optimal pool number `N`.
#' @rdname optimise_prevalence
#' @export
#'
#' @examples
#' optimise_s_prevalence(prevalence = 0.01, cost_unit = 5, cost_pool = 10)
#'
#' optimise_sN_prevalence(
#'   prevalence = 0.01, cost_unit = 5, cost_pool = 10,
#'   cost_cluster = 100, correlation = 0.05
#'   )
optimise_s_prevalence <- function(pool_number = 1,
                                  prevalence,
                                  cost_unit,
                                  cost_pool,
                                  cost_cluster = NA,
                                  correlation = NA,
                                  sensitivity = 1,
                                  specificity = 1,
                                  max_s = 50,
                                  form = "beta",
                                  interval = 0) {
  N <- pool_number

  if (form == "discrete") {
    stop('When form = "discrete" the cost of unit information function with
         respect to s often has mulitple minima and therefore the discrete
         distribution is not currently supported for optimisation')
  }
  if (interval < 0) stop("interval must be 0 or higher")

  invalid_cost <- FALSE # trigger for when costs are infinite to ensure that there's no cost output in these cases
  # print(c(theta = prevalence, sens = sensitivity, spec = specificity, unit = cost_unit, test = cost_pool, location =cost_cluster , rho = correlation, N = N, form = form, max_s = max_s))

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
  if (is.na(correlation)) {
    cost <- function(s) {
      ufc <- cost_fi(s, prevalence, sensitivity, specificity, cost_unit, cost_pool)
      ufc
    }
  } else {
    cost <- function(s) {
      ufc <- cost_fi_cluster(
        pool_size = s, pool_number = N, prevalence = prevalence,
        correlation = correlation,
        sensitivity = sensitivity,
        specificity = specificity,
        cost_unit = cost_unit, cost_pool = cost_pool,
        cost_cluster = cost_cluster, form = form
      )
      ufc
    }
  }

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

  if (s == max_s) warning("Maximum cost effectivness is achieved at or above the maximum size of pools allowed. Consider increasing max_s")

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

#' @rdname optimise_prevalence
#' @export
optimise_sN_prevalence <- function(prevalence,
                                   cost_unit,
                                   cost_pool,
                                   cost_cluster,
                                   correlation,
                                   sensitivity = 1,
                                   specificity = 1,
                                   max_s = 50,
                                   max_N = 20,
                                   form = "beta") {

  if (is.na(correlation) || correlation == 0) { # for simple random sampling there is no optimal N and for no correlation cluster survey the optimal approach is to infinite pools at one site (i.e. a simple random survey)
    opt <- optimise_s_prevalence(
      pool_number = 1, prevalence, cost_unit, cost_pool, cost_cluster,
      correlation = NA, sensitivity, specificity, max_s, form
    )

    opt$N <- ifelse(is.na(correlation), NA, Inf)
    opt$catch <- opt$N
  } else {
    Nopt <- 1
    opt <- list(cost = Inf)
    while (Nopt < max_N) {
      opt_new <- optimise_s_prevalence(
        pool_number = Nopt + 1, prevalence, cost_unit, cost_pool, cost_cluster,
        correlation, sensitivity, specificity, max_s, form
      )
      # print(opt_new)
      if (opt_new$cost > opt$cost) {
        break
      } else {
        Nopt <- Nopt + 1
        opt <- opt_new
      }
    }
    opt$N <- Nopt
    if (opt$N == max_N) {
      warning("Maximum cost effectivness is achieved at or above the maximum number of pools allowed. Consider increasing max_N")
    }
    if (opt$s == max_s) {
      warning("Maximum cost effectivness is achieved at or above the maximum size of pools allowed. Consider increasing max_s")
    }
  }

  opt
}

#' @rdname optimise_prevalence
#' @export
optimise_random_prevalence <- function(catch_mean, catch_variance,
                                       pool_strat_family,
                                       prevalence,
                                       cost_unit, cost_pool, cost_period, cost_cluster = NA,
                                       correlation = NA,
                                       sensitivity = 1, specificity = 1,
                                       max_period = 10, form = "beta",
                                       verbose = FALSE) {
  
  strat_par_names <- names(formals(pool_strat_family))
  pars <- as.list(rep(1, 1+ length(strat_par_names)))
  names(pars) <- c(".periods", strat_par_names)
  
  if(is.na(correlation)){
    f <- function(prm){
      cost_fi_random(nb_catch(catch_mean * prm$.periods, catch_variance * prm$.periods),
                     do.call(pool_strat_family, prm[strat_par_names]),
                     prevalence, sensitivity,specificity,
                     cost_unit, cost_pool, prm$.periods * cost_period)
    }
    
    opt <- optim_local_int(pars,f,verbose = verbose)
  }else{
    f <- function(prm){
      
      cost_fi_cluster_random(nb_catch(catch_mean * prm$.periods, catch_variance * prm$.periods),
                             do.call(pool_strat_family, prm[strat_par_names]),
                             prevalence, correlation,
                             sensitivity,specificity,
                             cost_unit, cost_pool,cost_cluster + prm$.periods * cost_period,
                             form = form)
    }
    
    opt <- optim_local_int(pars,f, verbose = verbose)
    
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
      val[i] <- fn(as.list(search[i,]))
    }
    
    min_val <- min(val)
    if(opt_val < min_val){
      break
    }else{
      opt_val <- min_val
      opt_par <- unlist(search[which.min(val),])
    }
  }
  return(list(val = opt_val, par = opt_par))
}

cost_fi <- function(pool_size, prevalence, sensitivity, specificity, cost_unit, cost_pool) {
  (sum(cost_unit * pool_size) + cost_pool) /
    fi_pool(pool_size, prevalence, sensitivity, specificity)
}

cost_fi_cluster <- function(pool_size, pool_number,
                            prevalence, correlation,
                            sensitivity, specificity,
                            cost_unit, cost_pool, cost_cluster,
                            form = "beta") {
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
                                   form = "beta"){

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
