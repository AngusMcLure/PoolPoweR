#' Constructor for pool_strat and pool_strat_family objects
#'
#' @param strat function that takes an integer number of units and returns the
#'   number and size of pools to split the units across
#' @param family family that the pool strategy function comes from
#' @param call call used to generate the pool_strat function. Used for format
#'   and as.character method (pretty printing)
#' @param description text description of pool strategy to be used in print
#'   function
#' @param name name of pool strategy family to be used in print function
#' @rdname pool_strat_methods

pool_strat <- function(strat, family, call, description){
  if(!inherits(strat, 'function')) {stop('strat must be a function')}
  ps <- strat
  attr(ps, 'call') <- call
  attr(ps, 'family') <- family
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
    cat('A pooling strategy that ', attr(x, 'description'), '\n')
  }
  invisible(x)
}

#' @export
pool_strat_family <- function(strat_family,name){
  psf <- strat_family
  attr(psf, 'name') <- name
  class(psf) <- c('pool_strat_family', 'function')
  psf
}

#' @export
print.pool_strat_family <- function(x,...){
  cat(attr(x, 'name'), '\n')
  invisible(x)
}

#' @export
format.pool_strat_family <- function(x,...){
  attr(x, 'name')
}

#' @export
as.character.pool_strat_family <- function(x,...){
  attr(x, 'name')
}



