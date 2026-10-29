#' Calibrate the VM (Vale-Maurelli) model
#'
#' Performs the calibration step of the Vale-Maurelli algorithm: for each
#' variable the Fleishman polynomial \eqn{a + bZ + cZ^2 + dZ^3} (with
#' \eqn{a = -c}) matching the target skewness and excess kurtosis is found,
#' and for each pair of variables the intermediate normal correlation that
#' yields the target correlation after transformation is solved from the
#' Vale-Maurelli cubic. Calibration is decoupled from simulation so the same
#' calibration object can be reused to draw many samples without recomputation.
#'
#' When the Fleishman equations have several solutions, a monotone increasing
#' polynomial is preferred. For many common skewness/kurtosis combinations
#' (e.g. skewness 1 and excess kurtosis 1.5) only polynomials that turn down
#' in the far tails exist; these are used, and flagged in \code{monotone}.
#'
#' @param sigma.target Target population covariance matrix.
#' @param skewness Target skewness vector.
#' @param excesskurtosis Target excess kurtosis vector.
#' @return An object of class \code{c("covsim_calib_vm", "covsim_calib")} that can
#'   be passed to \code{\link[stats]{simulate}} to draw samples. The element
#'   \code{coef} holds the Fleishman coefficients (one row per variable),
#'   \code{monotone} whether each polynomial is monotone increasing, and
#'   \code{z.corr} the intermediate correlation matrix.
#' @seealso \code{\link{rVM}}, \code{\link[stats]{simulate}}
#' @references Fleishman, A. I. (1978). A method for simulating non-normal
#'   distributions. Psychometrika, 43(4), 521-532.
#'
#'   Vale, C. D. and Maurelli, V. A. (1983). Simulating multivariate nonnormal
#'   distributions. Psychometrika, 48(3), 465-471.
#' @examples
#' sigma.target <- matrix(c(1, 0.5, 0.5, 1), 2)
#' cal <- calibrate_vm(sigma.target, skewness = c(1, 2), excesskurtosis = c(2, 7))
#' samples <- stats::simulate(cal, nsim = 5, N = 1000)
#' @export
calibrate_vm <- function(sigma.target, skewness, excesskurtosis) {
  if (is.null(sigma.target))
    stop("Please specify sigma.target")
  sigma.target <- as.matrix(sigma.target)
  nvar <- ncol(sigma.target)
  if (length(skewness) != nvar || length(excesskurtosis) != nvar)
    stop("skewness and excesskurtosis must have one entry per variable.")
  corr.target <- stats::cov2cor(sigma.target)

  coef <- t(vapply(seq_len(nvar), function(i) {
    fit <- fit_fleishman(skewness[i], excesskurtosis[i])
    if (is.null(fit))
      stop("No Fleishman polynomial exists for variable ", i, " (skewness ",
           skewness[i], ", excess kurtosis ", excesskurtosis[i],
           "). The excess kurtosis is too small relative to the skewness.")
    fit
  }, numeric(3)))
  colnames(coef) <- c("b", "c", "d")
  monotone <- apply(coef, 1, fleishman_monotone)
  coef <- cbind(a = -coef[, "c"], coef)

  z.corr <- diag(nvar)
  for (i in seq_len(nvar - 1)) {
    for (j in (i + 1):nvar) {
      rho <- solve_vm_corr(coef[i, ], coef[j, ], corr.target[i, j])
      if (is.na(rho))
        stop("No intermediate correlation exists for variables ", i, " and ", j,
             " (target correlation ", round(corr.target[i, j], 4), ").")
      z.corr[i, j] <- z.corr[j, i] <- rho
    }
  }

  min.eigen <- min(eigen(z.corr, symmetric = TRUE, only.values = TRUE)$values)
  if (min.eigen <= 0)
    stop("The intermediate correlation matrix is not positive definite ",
         "(smallest eigenvalue ", signif(min.eigen, 3), "). ",
         "The target covariance matrix cannot be reached with VM for these marginals.")

  dimnames(z.corr) <- dimnames(sigma.target)
  structure(
    list(
      coef = coef,
      monotone = monotone,
      z.corr = z.corr,
      sd = sqrt(diag(sigma.target)),
      sigma.target = sigma.target,
      skewness = skewness,
      excesskurtosis = excesskurtosis,
      nvar = as.integer(nvar),
      min.eigen = min.eigen
    ),
    class = c("covsim_calib_vm", "covsim_calib")
  )
}


#' Simulate from a calibrated VM model
#'
#' S3 method for \code{\link[stats]{simulate}} that draws samples from a
#' \code{covsim_calib_vm} object produced by \code{\link{calibrate_vm}}.
#'
#' @param object A \code{covsim_calib_vm} object.
#' @param nsim Number of independent samples (datasets) to return.
#' @param seed Optional seed passed to \code{\link{set.seed}}.
#' @param N Number of observations per simulated dataset.
#' @param ... Unused.
#' @return A list of \code{nsim} matrices, each with \code{N} rows.
#' @export
simulate.covsim_calib_vm <- function(object, nsim = 1, seed = NULL, N, ...) {
  if (missing(N))
    stop("Please specify N (number of observations per dataset).")
  if (!is.null(seed))
    set.seed(seed)

  cf <- object$coef
  lapply(seq_len(nsim), function(i) {
    z <- MASS::mvrnorm(N, rep(0, object$nvar), Sigma = object$z.corr)
    z <- matrix(z, nrow = N)
    x <- sapply(seq_len(object$nvar), function(j)
      cf[j, "a"] + z[, j] * (cf[j, "b"] + z[, j] * (cf[j, "c"] + z[, j] * cf[j, "d"])))
    x <- matrix(x, nrow = N) %*% diag(object$sd, nrow = object$nvar)
    colnames(x) <- colnames(object$sigma.target)
    x
  })
}


#' Simulation of non-normal data with the Vale-Maurelli method
#'
#' Convenience wrapper that combines \code{\link{calibrate_vm}} and
#' \code{\link[stats]{simulate}} in one call. Prefer the calibrate/simulate
#' split when drawing many samples from the same target.
#'
#' @inheritParams calibrate_vm
#' @param N Number of observations to simulate.
#' @param reps Number of simulated samples.
#' @return A list of simulated samples.
#' @references Vale, C. D. and Maurelli, V. A. (1983). Simulating multivariate
#'   nonnormal distributions. Psychometrika, 48(3), 465-471.
#' @examples
#' sigma.target <- matrix(c(1, 0.5, 0.5, 1), 2)
#' my.samples <- rVM(N = 1000, sigma.target, skewness = c(1, 2),
#'                   excesskurtosis = c(2, 7), reps = 5)
#' @export
rVM <- function(N, sigma.target, skewness, excesskurtosis, reps = 1) {
  cal <- calibrate_vm(sigma.target, skewness, excesskurtosis)
  stats::simulate(cal, nsim = reps, N = N)
}


## Fleishman helpers

# Moment equations for X = -c + bZ + cZ^2 + dZ^3: unit variance, skewness,
# excess kurtosis (Fleishman, 1978, eqs. 5, 11, 18)
fleishman_equations <- function(x, skew, kurt) {
  b <- x[1]; c <- x[2]; d <- x[3]
  c(b^2 + 6 * b * d + 2 * c^2 + 15 * d^2 - 1,
    2 * c * (b^2 + 24 * b * d + 105 * d^2 + 2) - skew,
    24 * (b * d + c^2 * (1 + b^2 + 28 * b * d) +
            d^2 * (12 + 48 * b * d + 141 * c^2 + 225 * d^2)) - kurt)
}

# b + 2cz + 3dz^2 > 0 for all z
fleishman_monotone <- function(x) {
  b <- x[1]; c <- x[2]; d <- x[3]
  b > 0 && ((d > 0 && c^2 < 3 * b * d) || (d == 0 && c == 0))
}

# Solve the Fleishman equations from several starting values. Returns c(b, c, d)
# with b > 0, preferring a monotone solution, or NULL if none is found.
fit_fleishman <- function(skew, kurt) {
  if (skew == 0 && kurt == 0)
    return(c(1, 0, 0))
  starts <- expand.grid(b = c(1, 0.9, 0.7, 0.5, 0.3), d = c(0, 0.02, 0.05, 0.1, 0.2, -0.05))
  starts$c <- skew / (2 * (starts$b^2 + 24 * starts$b * starts$d +
                             105 * starts$d^2 + 2))
  fallback <- NULL
  for (k in seq_len(nrow(starts))) {
    out <- nleqslv::nleqslv(c(starts$b[k], starts$c[k], starts$d[k]),
                            fleishman_equations, skew = skew, kurt = kurt,
                            control = list(ftol = 1e-12, xtol = 1e-14,
                                           maxit = 500))
    if (max(abs(out$fvec)) > 1e-9)
      next
    x <- out$x
    if (x[1] < 0) x[c(1, 3)] <- -x[c(1, 3)]  # (b, d) -> (-b, -d) is also a solution
    if (fleishman_monotone(x))
      return(x)
    if (is.null(fallback))
      fallback <- x
  }
  fallback
}

# Intermediate correlation from the Vale-Maurelli cubic (eq. 11):
# r = rho (b1 b2 + 3 b1 d2 + 3 d1 b2 + 9 d1 d2) + 2 c1 c2 rho^2 + 6 d1 d2 rho^3
solve_vm_corr <- function(coef1, coef2, target) {
  b1 <- coef1[["b"]]; c1 <- coef1[["c"]]; d1 <- coef1[["d"]]
  b2 <- coef2[["b"]]; c2 <- coef2[["c"]]; d2 <- coef2[["d"]]
  roots <- polyroot(c(-target,
                      b1 * b2 + 3 * b1 * d2 + 3 * d1 * b2 + 9 * d1 * d2,
                      2 * c1 * c2,
                      6 * d1 * d2))
  real <- Re(roots)[abs(Im(roots)) < 1e-8 & abs(Re(roots)) <= 1]
  if (length(real) == 0)
    return(NA_real_)
  real[which.min(abs(real - target))]
}
