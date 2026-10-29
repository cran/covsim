test_that("PLSIM works", {
  hsmodel <- "visual  =~ x1 + x2
              textual =~ x4 + x5"
  f <- lavaan::cfa(hsmodel, lavaan::HolzingerSwineford1939)
  sigma0 <- lavaan::lavInspect(f, "sigma.hat")
  target0 <- lavaan::lav_matrix_vech(sigma0)

  s <- rep(2, ncol(sigma0))
  k <- rep(7, ncol(sigma0))
  set.seed(1)
  sample1 <- rPLSIM(10^6,sigma0, s, k, reps=1 )[[1]][[1]]

  cov1 <- lavaan::lav_matrix_vech(cov(sample1))
  skew1 <- psych::skew(sample1)
  kurt1 <- unname(psych::kurtosi(sample1))


  expect_equal(cov1, target0, tolerance=0.01)
  expect_equal(skew1, s, tolerance=0.1)
  expect_equal(kurt1, k, tolerance=0.1)

})




test_that("PL covariance (Price's theorem) matches the truncated-moment formula", {
  # reference: the pre-1.2.0 implementation via rectangle probabilities and
  # truncated bivariate normal moments, here evaluated with stats::integrate
  # piecewise integration so each integrand is smooth
  pint <- function(f, cuts) {
    cuts <- c(-Inf, sort(cuts), Inf)
    sum(sapply(seq_len(length(cuts) - 1), function(k)
      stats::integrate(f, cuts[k], cuts[k + 1], rel.tol = 1e-11)$value))
  }
  ref_cov <- function(a1, b1, g1, a2, b2, g2, rho) {
    s <- sqrt(1 - rho^2)
    inner <- function(z1) sapply(z1, function(z)
      pint(function(u) pl_fun(rho * z + s * u, a2, b2, g2) * stats::dnorm(u),
           (g2 - rho * z) / s))
    pint(function(z) pl_fun(z, a1, b1, g1) * inner(z) * stats::dnorm(z), g1)
  }
  g <- get_gamma(4)
  a1 <- c(0.3, 0.8, 1.2, 2); a2 <- c(1.5, 0.4, 0.7, -0.2)
  b1 <- get_bs(a1, g); b2 <- get_bs(a2, g)
  for (rho in c(-0.95, -0.4, 0.3, 0.9))
    expect_equal(get_cov(a1, b1, g, a2, b2, g, rho),
                 ref_cov(a1, b1, g, a2, b2, g, rho), tolerance = 1e-6)
  # near |rho| = 1 the nested integrate() reference breaks down; compare with
  # the value from the old tmvtnorm-based implementation instead
  expect_equal(get_cov(a1, b1, g, a2, b2, g, 0.998), 0.492789, tolerance = 1e-5)
})


test_that("PLSIM intermediate correlations hit the target exactly", {
  set.seed(1)
  sigma.target <- cov(MASS::mvrnorm(10, rep(0, 4), diag(4)))
  cal <- calibrate_plsim(sigma.target, rep(1, 4), rep(4, 4), verbose = FALSE)
  corr.target <- cov2cor(sigma.target)
  for (i in 1:3) for (j in (i + 1):4)
    expect_equal(get_cov(cal$alist[[i]], cal$blist[[i]], cal$gammalist[[i]],
                         cal$alist[[j]], cal$blist[[j]], cal$gammalist[[j]],
                         cal$z.corr[i, j]),
                 corr.target[i, j], tolerance = 1e-8)
  expect_gt(cal$min.eigen, 0)
})


test_that("PLSIM reports the attainable range for infeasible correlations", {
  sigma.target <- matrix(c(1, -0.95, -0.95, 1), 2)
  expect_error(calibrate_plsim(sigma.target, c(2, 2), c(7, 7), verbose = FALSE),
               "attainable range")
})


test_that("univariate rPLSIM keeps the covsim 1.1.0 return format", {
  set.seed(1)
  res <- rPLSIM(200, matrix(2), 1, 3, reps = 2, verbose = FALSE)
  expect_named(res$model, c("a", "b", "gamma"))
  expect_true(is.numeric(res$samples[[1]]) && is.null(dim(res$samples[[1]])))
  expect_length(res$samples[[1]], 200)
})

test_that("PLSIM samples carry the variable names of sigma.target", {
  s <- matrix(0.5, 3, 3); diag(s) <- 1
  dimnames(s) <- list(c("a", "b", "c"), c("a", "b", "c"))
  cal <- calibrate_plsim(s, rep(1, 3), rep(4, 3), verbose = FALSE)
  expect_equal(colnames(simulate(cal, N = 10)[[1]]), c("a", "b", "c"))
  res <- rPLSIM(10, s, rep(1, 3), rep(4, 3), verbose = FALSE)
  expect_equal(colnames(res$samples[[1]]), c("a", "b", "c"))
  unnamed <- calibrate_plsim(unname(s), rep(1, 3), rep(4, 3), verbose = FALSE)
  expect_null(colnames(simulate(unnamed, N = 10)[[1]]))
})
