# assumes you have already done rstan::expose_stan_functions("quantile_functions.stan")
GLD_solver <- function(lower_quartile, median, upper_quartile, 
                       other_quantile, alpha, check = TRUE) {
  if (!exists("GLD_icdf")) stop('call rstan::expose_stan_functions("quantile_functions.stan")')
  obj_xi <- function(xi) {   # chi will be local
    GLD_icdf(alpha, median, IQR = upper_quartile - lower_quartile, chi, xi) - other_quantile
  }
  obj_chi <- function(chi) { # xi will be local (as will skewness)
    S_75 <- S_(0.75, chi, xi)
    S_25 <- S_(0.25, chi, xi)
    (S_75 + S_25 - 2 * S_(0.5, chi, xi) ) / (S_75 - S_25) - skewness
  }
  skewness <- (upper_quartile + lower_quartile - 2 * median) / 
              (upper_quartile - lower_quartile)
  chi <- skewness # good starting value if skewness is close to zero
  names(chi) <- NULL
  old_chi <- 2
  while (abs(chi - old_chi) > 1e-8) {
    xi <- suppressWarnings(try(uniroot(obj_xi, lower = 1e-2, upper = 1 - 1e-2)$root, 
                               silent = TRUE))
    if (!is.numeric(xi)) {
      if (check) stop("no GLD is possible with these quantiles")
      return(c(asymmetry = NA_real_, steepness = NA_real_))
    }
    old_chi <- chi
    chi <- suppressWarnings(try(uniroot(obj_chi, lower = -1 + 1e-6, upper = 1 - 1e-6)$root,
                                silent = TRUE))
    if (!is.numeric(chi)) {
      if (check) stop("no GLD is possible with these quantiles")
      return(c(asymmetry = NA_real_, steepness = NA_real_))
    }
  }
  if (alpha == 0 || alpha == 1) xi <- xi - .Machine$double.eps
  a_s <- c(asymmetry = chi, steepness = xi)
  if (check) {
    low <- GLD_icdf(0, median, IQR = upper_quartile - lower_quartile,
                  asymmetry = a_s[1], steepness = a_s[2])
    if (alpha > 0 && low > -Inf) warning("solution implies a bounded lower tail at ", low)
    high <- GLD_icdf(1, median, IQR = upper_quartile - lower_quartile, 
                     asymmetry = a_s[1], steepness = a_s[2])
    if (alpha < 1 && high < Inf) warning("solution implies a bounded upper tail at ", high)
  }
  return(a_s)
}

GLD_solver_bounded <- function(bounds, median, IQR, check = TRUE, ...) {
  if (!exists("GLD_icdf")) stop('call rstan::expose_stan_functions("quantile_functions.stan")')
  obj <- function(theta) {
    chi <- theta[1]
    xi  <- theta[2]
    if (abs(chi) > 1 || xi < 0 || xi > 1 || 
        xi > 0.5 * (1 + chi) || xi > 0.5 * (1 - chi)) return(NA_real_)
    out <- (bounds[1] - GLD_icdf(0, median, IQR, chi, xi)) ^ 2 +
           (bounds[2] - GLD_icdf(1, median, IQR, chi, xi)) ^ 2
    return(out)
  }
  start <- c(0, 0.5 - 1 / sqrt(5))
  # start is a solution that implies the GLD is uniform(a, b) where 
  # median == 0.5 * (a + b) and IQR == 0.5 * (b - a)
  sln <- suppressWarnings(nlm(obj, p = start, fscale = 0, ...))
  if (check && sln$code >= 4) {
    warning("solution for asymmetry and steepness is dubious",
            " try passing optional arguments to nlm via ... ")
  }
  a_s <- c(asymmetry = sln$estimate[1], steepness = sln$estimate[2])
  if (check && sln$minimum > 1e-15) {
    low  <- GLD_icdf(0, median, IQR, asymmetry = a_s[1], steepness = a_s[2])
    high <- GLD_icdf(1, median, IQR, asymmetry = a_s[1], steepness = a_s[2])
    warning("no asymmetry and steepness values achieve the bounds exactly; ",
            paste("actual bounds are", low, "and", high))
  }
  return(a_s)
}
