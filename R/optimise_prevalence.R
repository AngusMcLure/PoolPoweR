design_effect <- function(pool_size,
                          pool_number,
                          prevalence,
                          correlation,
                          sensitivity,
                          specificity,
                          form = "beta") {
  pool_number * pool_size * fi_pool(1, prevalence, sensitivity, specificity) *
    solve(fi_pool_cluster(
      pool_size,
      pool_number,
      prevalence,
      correlation,
      sensitivity,
      specificity,
      form
    ))[1, 1]
}

# optimising the pool size for estimating prevalence

cost_fi <- function(
    pool_size, prevalence, sensitivity, specificity, cost_unit, cost_pool) {
  (sum(cost_unit * pool_size) + cost_pool) /
    fi_pool(pool_size, prevalence, sensitivity, specificity)
}

cost_fi_cluster <- function(pool_size,
                            pool_number,
                            prevalence,
                            correlation,
                            sensitivity,
                            specificity,
                            cost_unit,
                            cost_pool,
                            cost_cluster,
                            form = "beta") {
  fi <- fi_pool_cluster(
    pool_size, pool_number, prevalence, correlation,
    sensitivity, specificity, form = form
  )
  cost <- sum(cost_unit * pool_size * pool_number) +
    sum(cost_pool * pool_number) + cost_cluster
  # print(c(cost = cost, information = fi))
  cost * solve(fi)[1, 1]
}

optimise_s_prevalence <- function(prevalence,
                                  cost_unit,
                                  cost_pool,
                                  cost_cluster = NA,
                                  correlation = NA,
                                  pool_number = 1,
                                  form = "beta",
                                  sensitivity = 1,
                                  specificity = 1,
                                  max_s = 50,
                                  interval = 0) {
  N <- pool_number

  if (form == "discrete") {
    stop('When form = "discrete" the cost of unit information function with
         respect to s often has mulitple minima and therefore the discrete
         distribution is not currently supported for optimisation')
  }
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

  opt <- optimise(cost, c(1, max_s))
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
      lower <- ceiling(uniroot(function(x) {
        cost(x) - max_cost
      }, c(1, s))$root)
    }
    if (cost(max_s) < max_cost) {
      upper <- max_s
      warning("A pool size greater than max_s may fall within the specified range of cost effectiveness. Consider increasing max_s")
    } else {
      upper <- floor(uniroot(function(x) {
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
  } else {
    stop("interval must be between 0 and 1.")
  }
  out
}

optimise_sN_prevalence <- function(prevalence, cost_unit, cost_pool,
                                   cost_cluster, correlation, form = "beta",
                                   sensitivity = 1, specificity = 1,
                                   max_s = 50, max_N = 20) {
  # print(c(theta = prevalence, sens = sensitivity, spec = specificity, unit = cost_unit, test = cost_pool, location =cost_cluster , rho = correlation, N = N, form = form, max.s = max.s))

  if (is.na(correlation) || correlation == 0) { # for simple random sampling there is no optimal N and for no correlation cluster survey the optimal approach is to infinite pools at one site (i.e. a simple random survey)
    opt <- optimise_s_prevalence(prevalence, cost_unit, cost_pool,
      cost_cluster,
      correlation = NA,
      pool_number = 1, form = form,
      sensitivity, specificity,
      max_s = max_s
    )

    opt$N <- ifelse(is.na(correlation), NA, Inf)
    opt$catch <- opt$N
  } else {
    Nopt <- 1
    opt <- list(cost = Inf)
    while (Nopt < max_N) {
      opt_new <- optimise_s_prevalence(
        prevalence,
        cost_unit, cost_pool, cost_cluster,
        correlation, Nopt + 1, form,
        sensitivity, specificity, max_s
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
      warning("Maximum cost effectivness is achieved at or above the maximum number of pools allowed. Consider increasing max.N")
    }
    if (opt$s == max_s) {
      warning("Maximum cost effectivness is achieved at or above the maximum size of pools allowed. Consider increasing max.s")
    }
  }

  opt
}
