#' Constructor for pool_strat objects
#'
#' @param strat function that takes an integer number of units and returns the
#'   number and size of pools to split the units across
#' @param call call used to generate the pool_strat function. Used for format
#'   and as.character method (pretty printing)
#' @param description text description of pool strategy to be used in print
#'   function
#' @rdname pool_strat_methods

pool_strat <- function(strat, call, description){
  if(!inherits(strat, 'function')) {stop('strat must be a function')}
  ps <- strat
  attr(ps, 'call') <- call
  attr(ps, 'description') <- description
  class(ps) <- c('pool_strat', 'function')
  return(ps)
}

#' @export
format.pool_strat <- function(x,...){
  return(format(attr(x,'call'),...))
}

#' @export
as.character.pool_strat <- function(x,...){
  return(format(attr(x,'call'),...))
}

#' @export
print.pool_strat <- function(x,...){
  if(is.null(attr(x, 'description'))){
    print.function(x)
  }else{
    cat(paste0('A pooling strategy that ', attr(x, 'description')))
  }
  invisible(x)
}


