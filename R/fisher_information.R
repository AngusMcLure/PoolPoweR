## Fisher information for population prevalence (and intra-cluster correlation
## where appropriate)
fi_pool <- function(pool_size, prevalence, sensitivity, specificity) {
  s <- pool_size
  theta <- prevalence
  varphi <- sensitivity
  psi <- specificity

  q <- 1 - theta
  numerator <- s^2 * (1 - psi - varphi)^2
  denominator <- (q^(2 - 2 * s) * (varphi + q^s * (1 - psi - varphi)) *
    (1 - varphi - q^s * (1 - psi - varphi)))

  numerator / denominator
}

fi_pool_cluster <- function(pool_size,
                            pool_number,
                            prevalence,
                            correlation,
                            sensitivity,
                            specificity,
                            form = "beta",
                            real.scale = FALSE) {
  s <- pool_size
  N <- pool_number
  K <- length(N)
  theta <- prevalence
  rho <- correlation
  varphi <- sensitivity
  psi <- specificity

  if (length(s) != length(N) || !all((N %% 1) == 0) || !all(c(N, s) > 0)) {
    stop("s and N must be vectors of positive numbers of common length. s can be non-integer, but N must be integer")
  }

  if (rho == 0) {
    # warning('correlation = 0 (no correlation) would mean that random sampling from a single location is equivalent to simple random sampling from the whole population (i.e. any one location is representative of the whole population). If this is the case, do not use cluster survey. However, this is highly unlikely to be true. Instead choose a small correlation')
    return(sum(fi_pool(s, theta, varphi, psi) * N))
  }
  if (rho == 1) {
    stop("correlation = 1 (perfect correlation) means that units from the same location are prefectly correlated. Do not use cluster surveys in this case.")
    # return(fi_pool(1,theta,varphi,psi))
  }
  if (K == 1 && N == 1 && s == 1) {
    return(matrix(c(fi_pool(s, theta, varphi, psi) * N, 0, 0, 0), nrow = 2))
  }
  phi <- function(theta) {
    varphi + (1 - psi - varphi) * (1 - theta)^s
  }

  if (form %in% c("discrete", "beta")) {
    # In these cases Fisher information matrix is calculated directly in terms
    # of theta and rho and checks on sums of likelihoods and likelihood
    # derivatives can be performed accordingly
    if (form == "discrete") {
      method <- "summation"
      if (length(s) != 1) stop("unequal pool size not implemented for kind = 'discrete'")
      # This models the prevalence at each site as a discrete distribution:
      # values:                0,   theta,           1
      # weights: rho * (1-theta), (1-rho), rho * theta
      # This is a pretty weird distribution and probably not a good way to model
      # site-level prevalence. However, it demonstrates that optimal design depends
      # on the distribution. Namely, under this model the most cost-effective design
      # will usually be to lots of individual tests at a few sites

      y <- 0:N
      phi. <- phi(theta)

      lik <- choose(N, y) * ((1 - rho) * phi.^y * (1 - phi.)^(N - y) +
        rho * theta * varphi^y * (1 - varphi)^(N - y) +
        rho * (1 - theta) * psi^(N - y) * (1 - psi)^y)
      lik_theta <- choose(N, y) * ((1 - rho) * (1 - psi - varphi) * s * (1 - theta)^(s - 1) * phi.^(y - 1) * (1 - phi.)^(N - y - 1) * (N * phi. - y) +
        rho * varphi^y * (1 - varphi)^(N - y) -
        rho * psi^(N - y) * (1 - psi)^y)
      lik_rho <- choose(N, y) * (-phi.^y * (1 - phi.)^(N - y) +
        theta * varphi^y * (1 - varphi)^(N - y) +
        (1 - theta) * psi^(N - y) * (1 - psi)^y)
    } else if (form == "beta") {
      rho. <- rho^-1 - 1
      Alpha <- theta * rho.
      Beta <- (1 - theta) * rho.
      if (varphi == 1 & psi == 1 & length(s) == 1 && (s == 1 | N < 20)) { # Cases where we have solutions in terms of beta functions etc. However this method becomes numerically unstable for large N
        y <- 0:N
        method <- "summation"
        if (s == 1) { # un-pooled case with a perfect test has simple exact solution
          # choose(N,y) * beta(Alpha + y, Beta + N - y)/beta(Alpha, Beta)
          lik <- exp(lchoose(N, y) + lbeta(Alpha + y, Beta + N - y) - lbeta(Alpha, Beta))

          lik_theta <- lik * rho. * (digamma(Alpha + y) - digamma(Alpha) +
            digamma(Beta) - digamma(Beta + N - y))

          lik_rho <- lik / rho^2 * (theta * (digamma(Alpha) - digamma(Alpha + y)) +
            (1 - theta) * (digamma(Beta) - digamma(Beta + N - y)) +
            digamma(Alpha + Beta + N) - digamma(Alpha + Beta))
        } else {
          lik <- sapply(y, function(y) {
            x <- 0:y
            choose(N, y) / beta(Alpha, Beta) * sum((-1)^x * choose(y, x) * beta(Alpha, Beta + s * (N - y + x)))
          })

          lik_theta <- sapply(y, function(y) {
            x <- 0:y
            (Alpha + Beta) * choose(N, y) / beta(Alpha, Beta) *
              sum((-1)^x * choose(y, x) * beta(Alpha, Beta + s * (N - y + x)) *
                (digamma(Beta) - digamma(Beta + s * (N - y + x))))
          })

          lik_rho <- sapply(y, function(y) {
            x <- 0:y
            (Alpha + Beta + 1)^2 * choose(N, y) / beta(Alpha, Beta) *
              sum((-1)^x * choose(y, x) * beta(Alpha, Beta + s * (N - y + x)) *
                ((digamma(Beta) - digamma(Beta + s * (N - y + x))) * (1 - theta) +
                  digamma(Alpha + Beta + s * (N - y + x)) - digamma(Alpha + Beta)))
          })
        }
      } else { # most general case (imperfect test, mix of pool sizes, or large N)
        # requires numerical integration
        method <- "integration"

        # All possible outcomes for N1,N2,... pools of sizes s1,s2,...
        ys <- do.call(expand.grid, purrr::map(N, ~ (0:.x)))

        # probability of observing y positive pools out of N is prod(choose(N,y)) times
        # the integral of integrand. However for alpha < 1 or beta <1 numerical
        # integration can be problematic so we subtract off the poles at p=0 and
        # p=1 with a multiple of a beta density. The resulting integrand can be
        # integrated numerically without the stability issues, however has to be
        # corrected by the integral of the beta density

        # integrand_theta and integrand_rho are the partial derivatives of
        # integrand with respect to theta and rho. It also has poles at p = 0 and
        # p = 1. We remove the poles with a similar trick as above, subtracting
        # off a multiple of a beta density and multiple of beta density times
        # log(p) or log(1-p). The latter still has a closed form expression for
        # the integral but it involves the generalised hypergeometric distribution
        # (see mathematica)

        if (Alpha >= 1 & Beta >= 1) {
          integrand <- function(p, y) { # note that this functions need to be vectorised for p, y, s, N. Vectorisation for p is done via loops, others through vectorised basic operations
            out <- dbeta(p, Alpha, Beta)
            for (j in 1:length(p)) {
              pj <- p[j]
              out[j] <- out[j] * prod(phi(pj)^y * (1 - phi(pj))^(N - y))
            }
            out
          }
          lik_correction <- function(y) {
            0
          }

          integrand_theta <- function(p, y) {
            rho. * integrand(p, y) * (log(p) - log1p(-p) - digamma(Alpha) + digamma(Beta))
          }
          lik_theta_correction <- function(y) {
            0
          }

          integrand_rho <- function(p, y) {
            rho^(-2) * integrand(p, y) * (theta * (digamma(Alpha) - log(p)) + (1 - theta) * (digamma(Beta) - log1p(-p)) - digamma(Alpha + Beta))
          }
          lik_rho_correction <- function(y) {
            0
          }
        } else if (Alpha < 1 || Beta < 1) {
          integrand <- function(p, y) { # note that these functions need to be vectorised for p to be passed to integrate
            out <- dbeta(p, Alpha, Beta)
            c0 <- prod((1 - psi)^y * psi^(N - y))
            c1 <- prod(varphi^y * (1 - varphi)^(N - y))
            for (j in 1:length(p)) {
              pj <- p[j]
              out[j] <- out[j] * (prod(phi(pj)^y * (1 - phi(pj))^(N - y)) - ifelse(pj < 0.5, c0, c1))
            }
            out
          }
          lik_correction <- function(y) {
            pbeta(0.5, Alpha, Beta) * prod(psi^(N - y) * (1 - psi)^y) +
              pbeta(0.5, Beta, Alpha) * prod((1 - varphi)^(N - y) * (varphi)^y)
          }

          integrand_theta <- function(p, y) { # note that this functions need to be vectorised for p
            out <- rho. * dbeta(p, Alpha, Beta)
            c0 <- prod((1 - psi)^y * psi^(N - y))
            c1 <- prod(varphi^y * (1 - varphi)^(N - y))
            for (j in 1:length(p)) {
              pj <- p[j]
              out[j] <- out[j] *
                (-ifelse(pj < 0.5,
                  c0 * (log(pj) - digamma(Alpha) + digamma(Beta)),
                  c1 * (-log1p(-pj) - digamma(Alpha) + digamma(Beta))
                ) +
                  prod(phi(pj)^y * (1 - phi(pj))^(N - y)) * (log(pj) - log1p(-pj) - digamma(Alpha) + digamma(Beta)))
            }
            out
          }
          lik_theta_correction <- function(y) {
            (pbeta(0.5, Alpha, Beta) * (log(0.5) + digamma(Beta) - digamma(Alpha)) -
              0.5^Alpha / (Alpha^2 * beta(Alpha, Beta)) *
                hypergeo::genhypergeo(c(Alpha, Alpha, 1 - Beta), c(Alpha + 1, Alpha + 1), 0.5)) *
              prod((1 - psi)^y * psi^(N - y)) * rho. +
              (pbeta(0.5, Beta, Alpha) * (-log(0.5) + digamma(Beta) - digamma(Alpha)) +
                0.5^Beta / Beta^2 / beta(Beta, Alpha) * hypergeo::genhypergeo(c(Beta, Beta, 1 - Alpha), c(Beta + 1, Beta + 1), 0.5)) *
                prod(varphi^y * (1 - varphi)^(N - y)) * rho.
          }

          integrand_rho <- function(p, y) { # note that this functions need to be vectorised for p
            out <- rho^(-2) * dbeta(p, Alpha, Beta)
            c0 <- prod((1 - psi)^y * psi^(N - y))
            c1 <- prod(varphi^y * (1 - varphi)^(N - y))
            for (j in 1:length(p)) {
              pj <- p[j]
              out[j] <- out[j] *
                (-ifelse(pj < 0.5,
                  c0 * (theta * (digamma(Alpha) - log(pj)) + (1 - theta) * digamma(Beta) - digamma(Alpha + Beta)),
                  c1 * (theta * digamma(Alpha) + (1 - theta) * (digamma(Beta) - log1p(-pj)) - digamma(Alpha + Beta))
                ) +
                  prod(phi(pj)^y * (1 - phi(pj))^(N - y)) * (theta * (digamma(Alpha) - log(pj)) + (1 - theta) * (digamma(Beta) - log1p(-pj)) - digamma(Alpha + Beta)))
            }
            out
          }
          lik_rho_correction <- function(y) {
            rho^(-2) * prod((1 - psi)^y * psi^(N - y)) *
              ((theta * (digamma(Alpha) - log(0.5)) + (1 - theta) * digamma(Beta) - digamma(Alpha + Beta)) * pbeta(0.5, Alpha, Beta) +
                theta * 0.5^Alpha / Alpha^2 / beta(Alpha, Beta) * hypergeo::genhypergeo(c(Alpha, Alpha, 1 - Beta), c(Alpha + 1, Alpha + 1), 0.5, series = TRUE)) +

              rho^(-2) * prod(varphi^y * (1 - varphi)^(N - y)) *
                ((theta * digamma(Alpha) + (1 - theta) * (digamma(Beta) - log(0.5)) - digamma(Alpha + Beta)) * pbeta(0.5, Beta, Alpha) +
                  (1 - theta) * 0.5^Beta / Beta^2 / beta(Beta, Alpha) * hypergeo::genhypergeo(c(Beta, Beta, 1 - Alpha), c(Beta + 1, Beta + 1), 0.5, series = TRUE))
          }
        }

        tol <- .Machine$double.eps^0.8
        lik <- apply(ys, 1, function(y) {
          plot(function(x) {
            integrand(x, y)
          }, main = paste("integrand", y), to = 0.1, n = 10000)
          (integrate(integrand, 0, 1, y = y, rel.tol = tol, abs.tol = tol)$value + lik_correction(y)) *
            prod(choose(N, y))
        })

        lik_theta <- apply(ys, 1, function(y) {
          plot(function(x) {
            integrand_theta(x, y)
          }, main = paste("integrand_theta", y), to = 0.1, n = 10000)
          (integrate(integrand_theta, 0, 1, y = y, rel.tol = tol, abs.tol = tol)$value + lik_theta_correction(y)) *
            prod(choose(N, y))
        })

        lik_rho <- apply(ys, 1, function(y) {
          plot(function(x) {
            integrand_rho(x, y)
          }, main = paste("integrand_rho", y, Alpha, Beta), to = 0.1, n = 10000)
          (integrate(integrand_rho, 0, 1, y = y, rel.tol = tol, abs.tol = tol)$value + lik_rho_correction(y)) *
            prod(choose(N, y))
        })
      }
    }
    tol <- 1e-5
    if (abs(sum(lik) - 1) > tol ||
      abs(sum(lik_theta)) > tol ||
      abs(sum(lik_rho)) > tol) {
      print(lik)
      print(sum(lik) - 1)
      print(lik_theta)
      print(sum(lik_theta))
      print(lik_rho)
      print(sum(lik_rho))
      print(s)
      print(N)
      print(theta)
      print(varphi)
      print(psi)
      print(correlation)
      print(form)
      stop("Error in ", method, " of likelihoods. Likelihoods do not add to 1 or derivatives of likelihood with respect to parameters do not sum to 0")
    }

    return(matrix(
      c(
        sum(lik_theta^2 / lik),
        rep(sum(lik_theta * lik_rho / lik), 2),
        sum(lik_rho^2 / lik)
      ),
      nrow = 2
    ))
  } else if (form %in% c("logitnorm", "cloglognorm")) {
    link <- switch(form,
      logitnorm = qlogis,
      cloglognorm = cloglog
    )
    invlink <- switch(form,
      logitnorm = plogis,
      cloglognorm = cloglog_inv
    )

    # calculate parameters on the real scale (mu, sigma of normal distribution)
    pars <- mu_sigma_linknorm(theta, theta * (1 - theta) * rho, link, invlink)
    mu <- pars[1]
    sigma <- pars[2]

    # All possible outcomes for N1,N2,... pools of sizes s1,s2,...
    ys <- do.call(expand.grid, purrr::map(N, ~ (0:.x)))

    # calculate Fisher information matrix for alternative parameters on the real scale
    integrand <- function(z, y) { # note that this functions need to be vectorised for p
      out <- dnorm(z, mean = mu, sd = sigma)
      p <- invlink(z)
      for (j in 1:length(p)) {
        pj <- p[j]
        out[j] <- out[j] * prod(phi(pj)^y * (1 - phi(pj))^(N - y))
      }
      out
    }

    integrand_mu <- function(z, y) {
      integrand(z, y) * (z - mu) / sigma
    }

    integrand_sigma <- function(z, y) {
      integrand(z, y) * ((z - mu)^2 - sigma^2) / sigma^3
    }

    tol <- .Machine$double.eps^0.8
    lik <- apply(ys, 1, function(y) {
      # plot(function(x){integrand(x,y)},main = paste('integrand', y), n = 10000)
      integrate(integrand, -Inf, Inf, y = y, rel.tol = tol, abs.tol = tol)$value *
        prod(choose(N, y))
    })

    lik_mu <- apply(ys, 1, function(y) {
      # plot(function(x){integrand_mu(x,y)},main = paste('integrand_mu', y), n = 10000)
      integrate(integrand_mu, -Inf, Inf, y = y, rel.tol = tol, abs.tol = tol)$value *
        prod(choose(N, y))
    })

    lik_sigma <- apply(ys, 1, function(y) {
      # plot(function(x){integrand_sigma(x,y)},main = paste('integrand_sigma', y), n = 10000)
      integrate(integrand_sigma, -Inf, Inf, y = y, rel.tol = tol, abs.tol = tol)$value *
        prod(choose(N, y))
    })

    # Check that integration is sufficiently accurate
    tol <- 1e-6
    if (abs(sum(lik) - 1) > tol ||
      abs(sum(lik_mu)) > tol ||
      abs(sum(lik_sigma)) > tol) {
      print(lik)
      print(sum(lik) - 1)
      print(lik_mu)
      print(sum(lik_mu))
      print(lik_sigma)
      print(sum(lik_sigma))
      print(s)
      print(N)
      print(theta)
      print(varphi)
      print(psi)
      print(correlation)
      print(form)
      stop("Error in integration of likelihoods. Likelihoods do not add to 1 or derivatives of likelihood with respect to parameters do not sum to 0")
    }

    FI <- matrix(
      c(
        sum(lik_mu^2 / lik),
        rep(sum(lik_mu * lik_sigma / lik), 2),
        sum(lik_sigma^2 / lik)
      ),
      nrow = 2
    )
    if (real.scale) { # If we want FI for the parameters on the real scale {mu,sigma} we're done!
      return(FI)
    } else { # calculate Jacobian for the change of parameters from {mu,sigma} -> {theta, rho}
      integrand_dtheta_dmu <- function(z) {
        dnorm(z, mean = mu, sd = sigma) * (z - mu) / sigma * invlink(z)
      }

      integrand_dtheta_dsigma <- function(z) {
        dnorm(z, mean = mu, sd = sigma) * ((z - mu)^2 - sigma^2) / sigma^3 * invlink(z)
      }

      # note that this function is the integrand for the dv/dmu where v non-centred second moment. To calculate drho/dmu we make corrections later in terms of theta, rho and dtheta/dmu
      integrand_drho_dmu <- function(z) {
        dnorm(z, mean = mu, sd = sigma) * (z - mu) / sigma * invlink(z)^2
      }
      # note that this function is the integrand for the dv/dsigma where v non-centred second moment. To calculate drho/dsigma we make corrections later in terms of theta, rho and dtheta/dsigma
      integrand_drho_dsigma <- function(z) {
        dnorm(z, mean = mu, sd = sigma) * ((z - mu)^2 - sigma^2) / sigma^3 * invlink(z)^2
      }

      tol <- .Machine$double.eps^0.8
      dtheta_dmu <- integrate(integrand_dtheta_dmu, lower = -Inf, upper = Inf, abs.tol = tol)$value
      dtheta_dsigma <- integrate(integrand_dtheta_dsigma, lower = -Inf, upper = Inf, abs.tol = tol)$value
      drho_dmu <- integrate(integrand_drho_dmu, lower = -Inf, upper = Inf, abs.tol = tol)$value
      drho_dmu <- (drho_dmu - dtheta_dmu * (rho * (1 - 2 * theta) + 2 * theta)) / (theta * (1 - theta))
      drho_dsigma <- integrate(integrand_drho_dsigma, lower = -Inf, upper = Inf, abs.tol = tol)$value
      drho_dsigma <- (drho_dsigma - dtheta_dsigma * (rho * (1 - 2 * theta) + 2 * theta)) / (theta * (1 - theta))

      # We actually want the the Jacobian for the inverse transformation so we invert the matrix
      J <- solve(matrix(
        c(
          dtheta_dmu, dtheta_dsigma,
          drho_dmu, drho_dsigma
        ),
        nrow = 2, byrow = TRUE
      ))
      # FI on the {theta, rho} parameterisation
      return(t(J) %*% FI %*% J)
    }
  } else {
    stop("accepted forms of the site prevalence distribution (argument form) are logitnorm, cloglognorm, beta, and discrete.")
  }
}
