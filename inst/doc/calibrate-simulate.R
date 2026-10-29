## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
# library(covsim)
# 
# # 1. calibrate once
# cal <- calibrate_ig(sigma.target, skewness, excesskurtosis)
# #   ... or calibrate_plsim(...)
# #   ... or calibrate_vm(...)
# #   ... or calibrate_vita(margins, sigma.target)
# 
# # 2. simulate one or more datasets
# draws <- simulate(cal, nsim = 100, N = 1000)   # list of 100 data matrices

## -----------------------------------------------------------------------------
library(covsim)
set.seed(1)
d <- 4
sigma.target <- cov(MASS::mvrnorm(50, mu = rep(0, d), Sigma = diag(d)))

cal_ig <- calibrate_ig(sigma.target,
                       skewness        = rep(1, d),
                       excesskurtosis  = rep(3, d))
class(cal_ig)

# draw five independent datasets from the same calibration
draws <- simulate(cal_ig, nsim = 5, N = 10000)
length(draws)
dim(draws[[1]])
round(cov(draws[[1]]) - sigma.target, 3)

## -----------------------------------------------------------------------------
set.seed(1)
sigma.target <- cov(MASS::mvrnorm(50, mu = rep(0, 4), Sigma = diag(4)))
cal_pl <- calibrate_plsim(sigma.target,
                          skewness       = rep(1, 4),
                          excesskurtosis = rep(3, 4),
                          verbose        = FALSE)

# Re-use the calibration for many replications
replications <- simulate(cal_pl, nsim = 50, N = 500)
length(replications)

## -----------------------------------------------------------------------------
set.seed(1)
sigma.target <- cov(MASS::mvrnorm(50, mu = rep(0, 4), Sigma = diag(4)))
cal_vm <- calibrate_vm(sigma.target,
                       skewness       = rep(1, 4),
                       excesskurtosis = rep(3, 4))
cal_vm$coef      # Fleishman coefficients, one row per variable
draws <- simulate(cal_vm, nsim = 5, N = 10000)
round(cov(draws[[1]]) - sigma.target, 3)

## -----------------------------------------------------------------------------
sigma <- matrix(0.5, 4, 4); diag(sigma) <- 1
dimnames(sigma) <- list(paste0("x", 1:4), paste0("x", 1:4))
dists <- list(IG = calibrate_ig(sigma, rep(2, 4), rep(7, 4)),
              VM = calibrate_vm(sigma, rep(2, 4), rep(7, 4)))
res <- simulate_estimate(dists, model = "F =~ x1 + x2 + x3 + x4",
                         n = 200, reps = 50, seed = 1,
                         extractor = function(fit) lavaan::fitMeasures(fit, "pvalue"))
aggregate(pvalue ~ dist, res, function(p) mean(p < 0.05))

## ----eval=FALSE---------------------------------------------------------------
# set.seed(1)
# d <- 5
# sigma.target <- cov(MASS::mvrnorm(50, mu = rep(0, d), Sigma = diag(d)))
# margins <- lapply(sqrt(diag(sigma.target)),
#                   function(s) list(distr = "norm", sd = s))
# 
# cal_vita <- calibrate_vita(margins, sigma.target,
#                            Nmax    = 1e5,
#                            cores   = parallel::detectCores(),
#                            verbose = FALSE)
# 
# draws <- simulate(cal_vita, nsim = 10, N = 1000)

