fleishman_moments <- function(coef) {
  p <- function(z) coef[["a"]] + coef[["b"]] * z + coef[["c"]] * z^2 + coef[["d"]] * z^3
  m <- sapply(1:4, function(k)
    stats::integrate(function(z) p(z)^k * stats::dnorm(z), -Inf, Inf,
                     rel.tol = 1e-12)$value)
  c(mean = m[1], var = m[2] - m[1]^2, skew = m[3] / m[2]^1.5,
    exkurt = m[4] / m[2]^2 - 3)
}

test_that("Fleishman coefficients reproduce the target moments exactly", {
  targets <- list(c(1.75, 3.75), c(1, 1.5), c(2, 7), c(3, 25), c(0, -1), c(-1, 2))
  for (tg in targets) {
    x <- fit_fleishman(tg[1], tg[2])
    coef <- c(a = -x[2], b = x[1], c = x[2], d = x[3])
    expect_equal(unname(fleishman_moments(coef)), c(0, 1, tg), tolerance = 1e-8)
    expect_gt(coef[["b"]], 0)
  }
  expect_null(fit_fleishman(2, 2))
})

test_that("a monotone Fleishman polynomial is preferred when one exists", {
  x <- fit_fleishman(2, 7)
  expect_true(fleishman_monotone(x))
})

test_that("VM intermediate correlation gives the target correlation", {
  cal <- calibrate_vm(matrix(c(1, 0.5, 0.5, 1), 2), c(1, 2), c(2, 7))
  rho <- cal$z.corr[1, 2]
  # bivariate Gauss-Hermite check of corr(p1(Z1), p2(Z2))
  n <- 60
  i <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  x <- e$values; w <- e$vectors[1, ]^2
  z1 <- rep(x, n); u <- rep(x, each = n); ww <- rep(w, n) * rep(w, each = n)
  z2 <- rho * z1 + sqrt(1 - rho^2) * u
  p <- function(cf, z) cf[["a"]] + cf[["b"]] * z + cf[["c"]] * z^2 + cf[["d"]] * z^3
  expect_equal(sum(ww * p(cal$coef[1, ], z1) * p(cal$coef[2, ], z2)), 0.5,
               tolerance = 1e-8)
})

test_that("calibrate_vm + simulate matches target moments", {
  hsmodel <- "visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6"
  f <- lavaan::cfa(hsmodel, lavaan::HolzingerSwineford1939)
  sigma0 <- lavaan::lavInspect(f, "sigma.hat")
  s <- rep(2, 6)
  k <- rep(7, 6)

  cal <- calibrate_vm(sigma0, s, k)
  expect_s3_class(cal, "covsim_calib_vm")
  expect_s3_class(cal, "covsim_calib")

  set.seed(1)
  sample1 <- stats::simulate(cal, nsim = 1, N = 10^6)[[1]]
  expect_equal(colnames(sample1), colnames(sigma0))
  expect_equal(lavaan::lav_matrix_vech(cov(sample1)),
               lavaan::lav_matrix_vech(sigma0), tolerance = 0.01)
  expect_equal(unname(psych::skew(sample1)), s, tolerance = 0.1)
  expect_equal(unname(psych::kurtosi(sample1)), k, tolerance = 0.1)
})

test_that("rVM returns a list of samples and VM errors are informative", {
  out <- rVM(N = 200, diag(3), skewness = c(0, 1, 1), excesskurtosis = c(0, 2, 2),
             reps = 3)
  expect_length(out, 3)
  expect_true(all(sapply(out, dim) == c(200, 3)))

  expect_error(calibrate_vm(diag(2), c(2, 0), c(2, 0)), "No Fleishman polynomial")
})
