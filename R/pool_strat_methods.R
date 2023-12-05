#' Methods for pool_strat obejcts
#'
#' @param x An object of class `pool_strat` (typically produced by one of the
#'   built-in functions for generating them such as `pool_max_size()` and
#'   `pool_taget_number()`
#' @rdname pool_strat_methods
#' @export



print.pool_strat <- function(x,...){
  if(is.null(attr(x, 'description'))){
    print.function(x)
  }else{
    print(paste0('A pooling strategy that ', attr(x, 'description')))
  }
}