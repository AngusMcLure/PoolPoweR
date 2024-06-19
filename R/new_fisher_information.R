#' sample_design constructor
#'
#' @param pool_size
#' @param pool_number
#' @param sensitivity
#' @param specificity
#'
#' @return An object of class \code{sample_design}
#' @export
#'
#' @examples
#' sd <- sample_design(
#'   pool_size = 10, pool_number = NULL, sensitivity = 1, specificity = 1
#' )
#' fi_pool(sd)
sample_design <- function(pool_size,
                          pool_number = NULL,
                          sensitivity = 1,
                          specificity = 1) {
  # TODO: input checks
  structure(
    list(
      pool_size = pool_size,
      pool_size = pool_number,
      sensitivity = sensitivity,
      specificity = specificity
    ),
    class = "sample_design"
  )
}

#' @export
fi_pool <- function(x, ...) {
  UseMethod("fi_pool")
}

#' @method fi_pool sample_design
#' @export
fi_pool.sample_design <- function(x, prevalence) {
  s <- x$pool_size
  theta <- prevalence
  varphi <- x$sensitivity
  psi <- x$specificity

  q <- 1 - theta
  s^2 * (1 - psi - varphi)^2 /
    (q^(2 - 2 * s) * (varphi + q^s * (1 - psi - varphi)) *
       (1 - varphi - q^s * (1 - psi - varphi)))
}
