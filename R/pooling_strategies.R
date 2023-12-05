#' Tools for defining pooling strategies
#'
#' `pool_max_size()` and `pool_taget_number()` are functions for defining the how a
#' given number of units should be divided into pools. Using `pool_max_size()`
#' pools have fixed maximum pool size (`max_size`) with any remainder placed in
#' a single smaller pool. Using `pool_target_number()` pools have a target number (target)
#' of approximately equal sized pools, with fewer pools used only if 

#'
#' @param max_size integer The maximum number of units per pool
#' @param target_number integer The variance of the number of units per cluster
#'   (catch)
#'
#' @return A function of class `pool_strat` with a single argument, `catch`, for
#'   the catch size that returns a list of the pool sizes and numbers that the
#'   catch should be split into
#' @rdname pooling_strategies
#' @export
#'
#' @examples
#' pool_max_size(10)
#'
#' pool_target_number(4)


pool_max_size <- function(max_size){
  strat <- function(catch){
    if(catch<max_size){
      return(list(pool_size = catch, pool_number = 1))
    }else if(catch%%max_size == 0){
      return(list(pool_size = max_size, pool_number = catch/max_size))
    }else{
      return(list(pool_size = c(catch%%max_size, max_size), pool_number = c(1, catch%/%max_size)))
    }
  }
  attr(strat, 'call') <- match.call()
  attr(strat, 'description') <- paste('that places units in pools of size', max_size, 'with any remainder placed in a single smaller pool.')
  class(strat) <- 'pool_strat'
  strat
}

#' @rdname pooling_strategies
#' @export
pool_target_number <- function(target_number){
  if(target_number%%1 !=0) stop('target_number must be an integer')
  strat <- function(catch){
    if(catch<target_number){
      return(list(pool_size = 1, pool_number = catch))
    }else if(catch%%target_number == 0){
      return(list(pool_size = catch/target_number, pool_number = target_number))
    }else{
      base_size <- catch%/%target_number
      base_number <- target_number - (catch %% target_number)
      return(list(pool_size = c(base_size, base_size + 1), pool_number = c(base_number, target_number-base_number)))
    }
  }
  attr(strat, 'call') <- match.call()
  attr(strat, 'description') <- paste(if(target_number == 1){'places all units in 1 pool,'}else{paste('aims to distribute units into', target_number, 'equally sized pools,')},'with no maximum pool size')
  class(strat) <- 'pool_strat'
  strat
}



  
# #An example of a more complicated pooling function generator. For the first
# #pool, the a fraction (specified) of total catch is placed into a pool (rounding
# #up). For the second pool, the same fraction of the remaining units (rounding
# #up) is placed into a pool. This is continued until there are no units left. The
# #total number of pools grows approximately logarithmically with catch size
# #(compare pool_max_size where it grows linearly with catch, and pool_fixed_N where
# #it is fixed)
# pool_fraction <- function(frac){
#   if(frac<=0 | frac>1) stop('frac must be greater than 0 or less than or equal to 1')
#   function(catch){
#     pool_size <- c()
#     count <- 0
#     while(catch>0){
#       count <- count + 1
#       pool_size[count] <- ceiling(frac * catch)
#       catch <- catch - pool_size[count]
#     }
#     pool_size <- table(pool_size)
#     target_number <- unname(pool_size)
#     pool_size <- as.numeric(names(pool_size))
#     list(pool_size = pool_size, target_number = target_number)
#   }
# }
# 
# #Like pool_fraction, but the fraction is determined by the ratio of the max pool size to catch. This is admitedly a little weird!
# pool_fraction_maxs <- function(maxs){
#   function(catch){
#     frac <- maxs/catch
#     pool_size <- c()
#     count <- 0
#     while(catch>0){
#       count <- count + 1
#       pool_size[count] <- ceiling(frac * catch)
#       catch <- catch - pool_size[count]
#     }
#     pool_size <- table(pool_size)
#     target_number <- unname(pool_size)
#     pool_size <- as.numeric(names(pool_size))
#     list(pool_size = pool_size, pool_number = target_number)
#   }
# }
