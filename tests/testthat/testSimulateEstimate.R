# estimates for both lavaan (S4 coef) and magmaan fits
est <- function(f) if (inherits(f, "lavaan")) lavaan::coef(f) else stats::coef(f)

se_setup <- function() {
  sigma <- matrix(0.5, 4, 4); diag(sigma) <- 1
  dimnames(sigma) <- list(paste0("x", 1:4), paste0("x", 1:4))
  list(model = "F =~ x1 + x2 + x3 + x4",
       dists = list(IG = calibrate_ig(sigma, rep(1, 4), rep(3, 4)),
                    VM = calibrate_vm(sigma, rep(1, 4), rep(3, 4))))
}

test_that("simulate_estimate returns one row per replication", {
  s <- se_setup()
  res <- simulate_estimate(s$dists, s$model, n = c(50, 100), reps = 3,
                           extractor = est, seed = 1)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2 * 2 * 3)
  expect_equal(names(res)[1:6], c("dist", "n", "rep", "converged", "error", "warning"))
  expect_setequal(unique(res$dist), c("IG", "VM"))
  expect_true(all(res$converged))
  expect_true("F=~x2" %in% names(res))
  expect_equal(attr(res, "info")$seed, 1)
})

test_that("simulate_estimate is reproducible and independent of cores", {
  skip_on_os("windows")
  s <- se_setup()
  chisq <- function(fit) lavaan::fitMeasures(fit, c("chisq", "pvalue"))
  r1 <- simulate_estimate(s$dists, s$model, n = 60, reps = 4, extractor = chisq,
                          seed = 7)
  r2 <- simulate_estimate(s$dists, s$model, n = 60, reps = 4, extractor = chisq,
                          seed = 7, cores = 2)
  expect_identical(r1, r2)
  r3 <- simulate_estimate(s$dists, s$model, n = 60, reps = 4, extractor = chisq,
                          seed = 8)
  expect_false(identical(r1$chisq, r3$chisq))
})

test_that("simulate_estimate leaves the caller's RNG state untouched", {
  s <- se_setup()
  set.seed(123)
  before <- .Random.seed
  kind <- RNGkind()
  simulate_estimate(s$dists$IG, s$model, n = 30, reps = 2, extractor = est, seed = 1)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
})

test_that("errors are recorded per replication and custom fitters work", {
  s <- se_setup()
  res <- simulate_estimate(s$dists$VM, s$model, n = 30, reps = 2, seed = 1,
                           fit = function(data, model) stop("boom"), extractor = est)
  expect_equal(res$error, c("boom", "boom"))
  expect_false(any(res$converged))

  res <- simulate_estimate(s$dists$VM, s$model, n = 30, reps = 2, seed = 1,
                           fit = function(data, model) colMeans(data),
                           extractor = function(f) f)
  expect_true(all(is.na(res$converged)))
  expect_equal(names(res)[-(1:6)], paste0("x", 1:4))
})

test_that("non-scalar extractor output becomes a list column", {
  s <- se_setup()
  res <- simulate_estimate(s$dists$IG, s$model, n = 30, reps = 2, seed = 1,
                           extractor = function(f) list(fit = f))
  expect_true(inherits(res$fit[[1]], "lavaan"))
})

test_that("simulate_estimate works with magmaan", {
  skip_if_not_installed("magmaan", minimum_version = "0.2.0")
  s <- se_setup()
  res <- simulate_estimate(s$dists, s$model, n = 100, reps = 2, seed = 1,
                           fit = "magmaan", extractor = est)
  expect_true(all(res$converged))
  expect_true(all(is.na(res$error)))
  # magmaan and lavaan ML estimates agree on the same data
  lav <- simulate_estimate(s$dists, s$model, n = 100, reps = 2, extractor = est,
                           seed = 1)
  expect_false(is.null(res[["F=~x2"]]))
  expect_equal(res[["F=~x2"]], lav[["F=~x2"]], tolerance = 1e-4)
})

test_that("simulate_estimate requires an extractor function", {
  s <- se_setup()
  expect_error(simulate_estimate(s$dists$IG, s$model, n = 30, reps = 1), "extractor")
  expect_error(simulate_estimate(s$dists$IG, s$model, n = 30, reps = 1,
                                 extractor = "chisq"), "extractor")
})

test_that("lavaan arguments reach lavaan and are not taken as seed", {
  s <- se_setup()
  opts <- function(fit) {
    o <- lavaan::lavInspect(fit, "options")
    list(se = o$se, start = o$start, iter.max = o$control$iter.max)
  }
  res <- simulate_estimate(s$dists$VM, s$model, n = 50, reps = 2, extractor = opts,
                           se = "robust.sem", start = "simple",
                           control = list(iter.max = 123))
  expect_true(all(is.na(res$error)))
  expect_equal(res$se, c("robust.sem", "robust.sem"))
  expect_equal(res$start, c("simple", "simple"))
  expect_equal(res$iter.max, c(123, 123))
})

test_that("arguments in ... with a fit function give a warning", {
  s <- se_setup()
  fit_ml <- function(data, model) lavaan::sem(model, data = data)
  expect_warning(simulate_estimate(s$dists$IG, s$model, n = 30, reps = 1, extractor = est,
                                   fit = fit_ml, estimator = "MLR"), "ignored")
})
