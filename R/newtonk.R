#' Compute Wave Number via Newton-Raphson Iteration
#'
#' Solves the wave dispersion equation for wave number (`k`) using the
#' Newton-Raphson method, given water depth and an initial deep-water wave number estimate.
#'
#' @param h Numeric. Water depth (in meters).
#' @param ko Numeric. Deep-water wave number estimate (radians per meter).
#'
#' @return Numeric. Converged wave number `k` (radians per meter).
#'
#' @details
#' This function iteratively solves the dispersion relation:
#' \deqn{k = \frac{\omega^2}{g \tanh(kh)}}
#' using Newton-Raphson iteration. The input `ko` is used as an initial guess
#' and refined to yield a more accurate wave number accounting for finite depth.
#' The method terminates if the relative change in successive estimates is below `1e-6`
#' or after 20 iterations.
#'
newtonk <- function(h, ko) {
  eps <- 1e-6

  k <- ko * tanh(ko * h)

  if (abs(k - ko) < eps) {
    return(k)
  }

  for (i in 1:20) {
    f <- ko * h - k * h * tanh(k * h)
    fp <- -h * (h * k * 1 / cosh((h * k) ^ 2) + tanh(h * k))
    kn <- k - f / fp

    if (abs(kn - k) / kn < eps) {
      return(kn)
    }

    k <- kn
  }

  return(kn)
}
