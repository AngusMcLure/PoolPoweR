# These functions assume that traps are set out of for an integer number, P,
# periods and the insects are collected at the end of each period. Insects are
# pooled across the period.

# As currently implemented:
#    * The number of insects caught per unit time is fixed and known
#    * The pooling of insects takes the the total number and divides into equal sized pools (giving non-integer pool sizes)
#    * Assumes that the global minima for the cost per unit information function
#      is the closest minima to N = 1 and P = 1. I.e. FOr fixed P, systematically
#      works its way with increasing N until it finds the first local minima,
#      and for optimised P, repeats this for increasing P until it finds another
#      local minima
#
# All of these could perhaps be improved:
#    * Use fi_pool_imperfect_cluster_unequal to allow for random catch sizes
#    * Use fi_pool_imperfect_cluster_unequal or just pool strategy objects to split up catches between pools
#    * Use a local integer search after finding that first local minima?




#Experimental function that allows for random catches per site (following
#arbitrary distribution defined with catch.dist using a custom class), and
#arbitrary pooling strategies to split these catches between the pools (defined
#in pool.strat which is a function which takes an integer, and returns a set of
#pool sizes (s) and pool numbers (N)). An example catch distribution generating
#function can be found in utils (nb_catch). A few simple example pool.strat
#function generating functions can be found in utils also

fi_pool_imperfect_cluster_unequal <- function(catch.dist, pool.strat, prevalence, sensitivity, specificity,
                                              correlation,form = 'beta', real.scale = FALSE, max.iter = 200){
  
  catch <- max(catch.dist$min-1, 0)
  terminate <- FALSE
  FI <- matrix(0, 2,2)
  rel.tol <- 1e-4
  iter <- 0
  FI.incr <- array(dim = c(2,2,max.iter))
  catches <- c()
  while(!terminate){
    catch <- catch + 1
    mass <- catch.dist$pmf(catch)
    if(mass == 0) next #this avoids unnecessary calls to fi_pool_imperfect_cluster and prevents the early termination of the algorithm for distributions that may have 0 mass for some n but non-zero mass for m>n (e.g. if distribution only has mass on multiples of 10)
    iter <- iter + 1 #only counts iteration if mass is non-zero
    pooling <- pool.strat(catch)
    catches[iter] <- catch
    FI.incr[,,iter] <- mass *
      fi_pool_imperfect_cluster(pooling$s, pooling$N, prevalence,
                                sensitivity,specificity,correlation,
                                form,real.scale)
    FI <- FI +  FI.incr[,,iter]
    # Stop if increment changes ALL elements of FI by less than fraction rel.tol or
    # all elements by less than abs.tol
    rel.incr <- abs(FI.incr[,,iter]/FI)
    #print(rel.incr)
    if(all(rel.incr <= rel.tol) | catch == catch.dist$max){
      terminate <- TRUE
    }
    if(iter == max.iter){
      terminate <- TRUE
      warning('Algorithm reached max.iter without converging')
    }
  }
  catches <- catches[1:iter]
  FI.incr <- FI.incr[,,1:iter]
  plot(catches, FI.incr[1,1,])
  plot(catches, FI.incr[1,2,])
  plot(catches, FI.incr[2,2,])
  
  FI
}





optimise_N_prevalence <- function(prevalence, cost.unit, cost.pool,
                                  cost.location, cost.collect, catch.collect,
                                  correlation, P, form = 'beta',
                                  sensitivity = 1, specificity = 1,
                                  max.s = 50, max.N = 20){
  # For given set up including number of sampling periods (P) and catch per period (i.e. given samples per location)
  # calculate the number of pools that is optimal. Resulting pool size will usually not be an integer
  n <- P * catch.collect
  N.opt <- max(if(correlation == 0){1}else{2}, ceiling(n/max.s)) - 1
  cost.opt <- Inf
  while(N.opt<n){
    cost.new <- cost_fi_cluster(n/(N.opt+1),N.opt+1,prevalence,correlation,
                                sensitivity,specificity,
                                cost.unit + cost.collect/catch.collect,
                                cost.pool,cost.location,form)
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

optimise_NP_prevalence <- function(prevalence, cost.unit, cost.pool,
                                   cost.location, cost.collect, catch.collect,
                                   correlation, form = 'beta',
                                   sensitivity = 1, specificity = 1,
                                   max.s = 50, max.P = 20){
  #print(c(theta = prevalence, sens = sensitivity, spec = specificity, unit = cost.unit, test = cost.pool, location = cost.location , rho = correlation, N = N, form = form, max.s = max.s))
  theta <- prevalence
  
  if(correlation == 0){
    stop('If there is no correlation between units at locations, then this means sampling at a single random location is a identical to sampling from the whole population. This is assumption unlikely to be true in most settings, but would mean that sampling at a single site is the most cost-effective strategy')
  }else{
    P.opt <- 0
    opt <- list(cost = Inf)
    while(P.opt < max.P){
      opt.new <- optimise_N_prevalence(prevalence, cost.unit, cost.pool,
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


