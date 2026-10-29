## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE
)
has_magmaan <- requireNamespace("magmaan", quietly = TRUE) &&
  utils::packageVersion("magmaan") >= "0.2.0"

## -----------------------------------------------------------------------------
library(covsim)

model <- "F1 =~ x1 + x2 + x3
          F2 =~ x4 + x5 + x6"

L <- matrix(0, 6, 2)
L[1:3, 1] <- L[4:6, 2] <- 0.7
sigma <- L %*% matrix(c(1, 0.4, 0.4, 1), 2) %*% t(L)
diag(sigma) <- 1
dimnames(sigma) <- list(paste0("x", 1:6), paste0("x", 1:6))

skew <- rep(2, 6)
kurt <- rep(7, 6)   # excess kurtosis

## -----------------------------------------------------------------------------
dists <- list(
  IG    = calibrate_ig(sigma, skew, kurt),
  PLSIM = calibrate_plsim(sigma, skew, kurt, verbose = FALSE),
  VM    = calibrate_vm(sigma, skew, kurt)
)

## -----------------------------------------------------------------------------
moments <- function(x) c(skew = mean(scale(x)^3), exkurt = mean(scale(x)^4) - 3)
set.seed(1)
round(sapply(dists, function(d) moments(simulate(d, N = 1e5)[[1]][, 1])), 2)

## ----eval = FALSE-------------------------------------------------------------
# margins <- rep(list(list(distr = "gamma", shape = 1)), 6)  # skewness 2
# dists$VITA <- calibrate_vita(margins, sigma, Nmax = 1e5, verbose = FALSE)

## -----------------------------------------------------------------------------
samples <- simulate(dists$PLSIM, nsim = 50, N = 200, seed = 1)
length(samples)
head(samples[[1]], 3)

## -----------------------------------------------------------------------------
set.seed(1)
old_way <- rPLSIM(200, sigma, skew, kurt, reps = 50, verbose = FALSE)
length(old_way$samples)

## -----------------------------------------------------------------------------
fit_tests <- function(x) {
  fit <- tryCatch(lavaan::sem(model, data = as.data.frame(x), estimator = "MLR"),
                  error = function(e) NULL)
  if (is.null(fit))
    return(c(pvalue = NA, pvalue.scaled = NA))
  lavaan::fitMeasures(fit, c("pvalue", "pvalue.scaled"))
}

by_hand <- lapply(dists, function(cal) {
  samples <- simulate(cal, nsim = 50, N = 200, seed = 1)
  t(sapply(samples, fit_tests))
})
head(by_hand$VM, 3)

## -----------------------------------------------------------------------------
round(sapply(by_hand, function(p) colMeans(p < 0.05, na.rm = TRUE)), 2)

## -----------------------------------------------------------------------------
chisq_tests <- function(fit) {
  lavaan::fitMeasures(fit, c("chisq", "pvalue", "chisq.scaled", "pvalue.scaled"))
}

res <- simulate_estimate(dists, model, n = c(200, 500), reps = 50,
                         extractor = chisq_tests, estimator = "MLR", seed = 1)
head(res)

## -----------------------------------------------------------------------------
table(res$converged, useNA = "ifany")

## -----------------------------------------------------------------------------
reject <- function(p) mean(p < 0.05, na.rm = TRUE)
aggregate(cbind(ML = pvalue, MLR = pvalue.scaled) ~ dist + n, data = res,
          FUN = reject)

## -----------------------------------------------------------------------------
fit_indices <- function(fit) {
  lavaan::fitMeasures(fit, c("rmsea", "rmsea.robust", "cfi", "cfi.robust"))
}

res_fi <- simulate_estimate(dists, model, n = 200, reps = 50,
                            extractor = fit_indices, estimator = "MLR", seed = 1)
aggregate(cbind(rmsea, rmsea.robust, cfi, cfi.robust) ~ dist, data = res_fi,
          FUN = mean)

## -----------------------------------------------------------------------------
loading_x2 <- function(fit) {
  pe <- lavaan::parameterEstimates(fit)
  pe <- pe[pe$lhs == "F1" & pe$op == "=~" & pe$rhs == "x2", ]
  c(est = pe$est, se = pe$se, lower = pe$ci.lower, upper = pe$ci.upper)
}

res_ML  <- simulate_estimate(dists, model, n = 200, reps = 50,
                             extractor = loading_x2, seed = 2)
res_MLR <- simulate_estimate(dists, model, n = 200, reps = 50,
                             extractor = loading_x2, estimator = "MLR", seed = 2)

summarize <- function(r) {
  c(bias = mean(r$est) - 1,
    empirical_sd = sd(r$est),
    mean_se = mean(r$se),
    coverage = mean(r$lower < 1 & 1 < r$upper))
}
tab <- rbind(sapply(split(res_ML, res_ML$dist), summarize),
             sapply(split(res_MLR, res_MLR$dist), summarize))
rownames(tab) <- paste(rep(c("ML", "MLR"), each = 4), rownames(tab))
round(tab, 3)

## -----------------------------------------------------------------------------
attr(res, "info")

## ----eval = FALSE-------------------------------------------------------------
# res4 <- simulate_estimate(dists, model, n = c(200, 500), reps = 50,
#                           extractor = chisq_tests, estimator = "MLR", seed = 1,
#                           cores = 4)
# all.equal(res, res4)   # TRUE (on Windows, cores is set to 1)

## -----------------------------------------------------------------------------
has_magmaan

## ----eval = has_magmaan-------------------------------------------------------
# mag_tests <- function(fit) {
#   tests <- summary(fit)$tests
#   setNames(tests$pvalue, paste(tests$test, tests$reference, sep = "_"))
# }
# 
# res_mag <- simulate_estimate(dists, model, n = c(200, 500), reps = 50,
#                              extractor = mag_tests, fit = "magmaan", seed = 1)
# head(res_mag)
# aggregate(cbind(score_peba4, lr_peba4) ~ dist + n, data = res_mag,
#           FUN = reject)

## ----eval = has_magmaan-------------------------------------------------------
# both <- merge(res[, c("dist", "n", "rep", "pvalue.scaled")],
#               res_mag[, c("dist", "n", "rep", "score_peba4")])
# aggregate(cbind(MLR = pvalue.scaled, PEBA4 = score_peba4) ~ dist, data = both,
#           FUN = reject)

## ----eval = has_magmaan-------------------------------------------------------
# mag_loading_x2 <- function(fit) {
#   cf <- summary(fit)$coefficients
#   cf <- cf[cf$lhs == "F1" & cf$op == "=~" & cf$rhs == "x2", ]
#   c(est = cf$est, se = cf$se, lower = cf$ci.lower, upper = cf$ci.upper)
# }
# mag_fit_indices <- function(fit) {
#   fm <- magmaan::fit_measures(fit, lavaan_compat = "MLR")
#   setNames(fm$estimate, fm$index)[c("rmsea", "rmsea.robust", "cfi", "cfi.robust")]
# }
# 
# res_mag2 <- simulate_estimate(dists, model, n = 200, reps = 50, fit = "magmaan",
#                               extractor = function(fit) c(mag_loading_x2(fit),
#                                                           mag_fit_indices(fit)),
#                               seed = 2)
# round(sapply(split(res_mag2, res_mag2$dist), summarize), 3)

## ----eval = has_magmaan-------------------------------------------------------
# t_lavaan <- system.time(
#   simulate_estimate(dists$VM, model, n = 500, reps = 50,
#                     extractor = chisq_tests, estimator = "MLR", seed = 1))[["elapsed"]]
# t_magmaan <- system.time(
#   simulate_estimate(dists$VM, model, n = 500, reps = 50,
#                     extractor = mag_tests, fit = "magmaan", seed = 1))[["elapsed"]]
# c(lavaan = t_lavaan, magmaan = t_magmaan, ratio = t_lavaan / t_magmaan)

## ----eval = FALSE-------------------------------------------------------------
# remotes::install_github("JonasMoss/magmaan", subdir = "r-package")
# remotes::install_github("JonasMoss/magmaan", subdir = "r-magmaan")

## -----------------------------------------------------------------------------
fit_wls <- function(data, model) lavaan::sem(model, data = data, estimator = "WLS")
res_wls <- simulate_estimate(dists$IG, model, n = 500, reps = 20,
                             extractor = chisq_tests, fit = fit_wls, seed = 1)
reject(res_wls$pvalue)

