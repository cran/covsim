#' Simulate data and estimate a model repeatedly
#'
#' Runs a Monte Carlo study: for every calibrated distribution in \code{dist}
#' and every sample size in \code{n}, draws \code{reps} samples with
#' \code{\link[stats]{simulate}}, fits \code{model} to each sample, and
#' collects what \code{extractor} returns into one data frame.
#'
#' Estimation is done with lavaan (\code{fit = "lavaan"}, the default), with
#' magmaan (\code{fit = "magmaan"}, if installed), or with any function
#' \code{function(data, model)}. With magmaan the model is constructed once
#' per distribution with \code{magmaan::magmaan_model()} and reused for every
#' sample, which avoids repeated parsing and setup.
#'
#' Each replication draws its data from its own L'Ecuyer-CMRG random number
#' stream, determined by \code{seed} and the position of the replication in
#' the design. Results are therefore identical for any value of \code{cores}.
#' Errors and warnings from fitting or extraction are caught and recorded per
#' replication rather than stopping the run.
#'
#' @param dist A calibration object (from \code{\link{calibrate_ig}},
#'   \code{\link{calibrate_plsim}}, \code{\link{calibrate_vm}} or
#'   \code{\link{calibrate_vita}}), or a named list of such objects to compare
#'   several distributions.
#' @param model Model syntax (lavaan syntax for \code{fit = "lavaan"} and
#'   \code{fit = "magmaan"}), or whatever a user-supplied \code{fit} function
#'   expects. Variable names must match the column names of
#'   \code{sigma.target} used in calibration.
#' @param n Vector of sample sizes.
#' @param reps Number of replications per distribution and sample size.
#' @param extractor A function of one fitted object that returns what is
#'   to be kept from that replication, as a named vector or list of scalars.
#'   Each element becomes a column. The researcher decides what to extract;
#'   see the examples for parameter estimates, standard errors, p-values, fit
#'   measures and test statistics with lavaan and magmaan. Use
#'   \code{function(fit) list(fit = fit)} to keep the fitted objects
#'   themselves (as a list column).
#' @param ... Further arguments passed to \code{lavaan::sem()} or
#'   \code{magmaan::magmaan()}, for example \code{estimator = "MLR"},
#'   \code{se = "robust.sem"}, \code{start = "simple"} or
#'   \code{control = list(iter.max = 500)}. Not used with a user-supplied
#'   \code{fit} function, which can set such options itself.
#' @param fit \code{"lavaan"}, \code{"magmaan"}, or a function
#'   \code{function(data, model)} returning a fitted object.
#' @param cores Number of cores. Values above 1 use
#'   \code{\link[parallel]{mclapply}} (not available on Windows, where the
#'   replications run sequentially).
#' @param seed Optional integer seed. If \code{NULL}, a seed is drawn from the
#'   current random number generator.
#' @return A data frame with one row per replication and columns
#'   \code{dist}, \code{n}, \code{rep}, \code{converged}, \code{error},
#'   \code{warning}, followed by the extracted values. The seed and package
#'   versions are stored in \code{attr(, "info")}.
#' @section Extractors for magmaan:
#' A magmaan fit \code{f} gives the parameter table with
#' \code{summary(f)$coefficients} (columns \code{lhs}, \code{op}, \code{rhs},
#' \code{free}, \code{est}, \code{se}, \code{pvalue}, ...), its score and
#' likelihood ratio tests with \code{summary(f)$tests}, and fit measures with
#' \code{magmaan::fit_measures(f)}. Note that by default \code{fit_measures()}
#' reports robust RMSEA and CFI and a bias-corrected SRMR under the plain
#' names; \code{magmaan::fit_measures(f, lavaan_compat = "MLR")} gives
#' lavaan's definitions and names (\code{rmsea}, \code{rmsea.robust},
#' \code{cfi.robust}, ...). magmaan's default standard errors are robust
#' (sandwich), matching lavaan with \code{estimator = "MLR"} rather than
#' plain ML. magmaan fits a mean structure, so its parameter table also
#' contains the intercepts.
#' @seealso \code{\link{calibrate_ig}}, \code{\link{calibrate_plsim}},
#'   \code{\link{calibrate_vm}}, \code{\link{calibrate_vita}}
#' @examples
#' model <- "F =~ x1 + x2 + x3 + x4"
#' sigma <- matrix(0.5, 4, 4); diag(sigma) <- 1
#' dimnames(sigma) <- list(paste0("x", 1:4), paste0("x", 1:4))
#' dists <- list(
#'   IG = calibrate_ig(sigma, skewness = rep(2, 4), excesskurtosis = rep(7, 4)),
#'   VM = calibrate_vm(sigma, skewness = rep(2, 4), excesskurtosis = rep(7, 4)))
#'
#' # rejection rate of the ML chi-square test
#' chisq <- function(fit) lavaan::fitMeasures(fit, c("chisq", "df", "pvalue"))
#' res <- simulate_estimate(dists, model, n = 200, reps = 20,
#'                          extractor = chisq, seed = 1)
#' aggregate(pvalue ~ dist, data = res, FUN = function(p) mean(p < 0.05))
#'
#' # RMSEA and CFI, ordinary and robust (MLR)
#' fitind <- function(fit)
#'   lavaan::fitMeasures(fit, c("rmsea", "cfi", "rmsea.robust", "cfi.robust"))
#' res <- simulate_estimate(dists, model, n = 200, reps = 20, estimator = "MLR",
#'                          extractor = fitind, seed = 1)
#' aggregate(cbind(rmsea, rmsea.robust) ~ dist, data = res, FUN = mean)
#'
#' # estimate, standard error and p-value of one loading
#' loading <- function(fit) {
#'   pe <- lavaan::parameterEstimates(fit)
#'   pe <- pe[pe$lhs == "F" & pe$op == "=~" & pe$rhs == "x2", ]
#'   c(est = pe$est, se = pe$se, pvalue = pe$pvalue)
#' }
#' res <- simulate_estimate(dists$VM, model, n = 200, reps = 20,
#'                          extractor = loading, seed = 1)
#' colMeans(res[, c("est", "se")])   # mean estimate vs. mean standard error
#' sd(res$est)                       # empirical standard error
#'
#' # the same with magmaan
#' if (requireNamespace("magmaan", quietly = TRUE) &&
#'     utils::packageVersion("magmaan") >= "0.2.0") {
#'   mag_loading <- function(fit) {
#'     cf <- summary(fit)$coefficients
#'     cf <- cf[cf$lhs == "F" & cf$op == "=~" & cf$rhs == "x2", ]
#'     c(est = cf$est, se = cf$se, pvalue = cf$pvalue)
#'   }
#'   mag_tests <- function(fit) {
#'     tests <- summary(fit)$tests
#'     stats::setNames(tests$pvalue, paste(tests$test, tests$reference, sep = "_"))
#'   }
#'   res <- simulate_estimate(dists$VM, model, n = 200, reps = 20, fit = "magmaan",
#'                            extractor = function(fit) c(mag_loading(fit), mag_tests(fit)),
#'                            seed = 1)
#'   colMeans(res[, c("est", "se", "score_peba4")])
#' }
#' @export
simulate_estimate <- function(dist, model, n, reps, extractor, ...,
                              fit = "lavaan", cores = 1, seed = NULL) {
  # fit, cores and seed come after `...` so that they must be named in full;
  # otherwise lavaan's `se` would partially match `seed`
  if (missing(extractor) || !is.function(extractor))
    stop("`extractor` must be a function of a fitted object, for example ",
         "function(fit) lavaan::fitMeasures(fit, c(\"chisq\", \"pvalue\")). ",
         "See ?simulate_estimate for examples.")
  if (inherits(dist, "covsim_calib"))
    dist <- list(dist)
  if (!is.list(dist) || length(dist) == 0 ||
      !all(vapply(dist, inherits, logical(1), "covsim_calib")))
    stop("`dist` must be a calibration object or a list of calibration objects.")
  if (is.null(names(dist)))
    names(dist) <- if (length(dist) == 1) class(dist[[1]])[1] else
      paste0("dist", seq_along(dist))
  if (anyNA(n) || any(n < 1) || length(reps) != 1 || reps < 1)
    stop("`n` must be positive sample sizes and `reps` a positive number.")

  fitters <- lapply(dist, function(d) make_fitter(fit, model, d, ...))

  design <- expand.grid(rep = seq_len(reps), n = n, dist = names(dist),
                        stringsAsFactors = FALSE)[, c("dist", "n", "rep")]

  # one L'Ecuyer-CMRG stream per replication, independent of `cores`
  if (is.null(seed))
    seed <- sample.int(.Machine$integer.max, 1)
  old.kind <- RNGkind()
  old.seed <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    RNGkind(old.kind[1], old.kind[2], old.kind[3])
    if (is.null(old.seed)) rm(".Random.seed", envir = globalenv())
    else assign(".Random.seed", old.seed, envir = globalenv())
  })
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  streams <- vector("list", nrow(design))
  s <- get(".Random.seed", envir = globalenv())
  for (k in seq_len(nrow(design))) {
    streams[[k]] <- s
    s <- parallel::nextRNGStream(s)
  }

  run_one <- function(k) {
    assign(".Random.seed", streams[[k]], envir = globalenv())
    d <- dist[[design$dist[k]]]
    data <- as_sim_data(stats::simulate(d, nsim = 1, N = design$n[k])[[1]], d)
    warn <- character(0)
    out <- withCallingHandlers(
      tryCatch({
        f <- fitters[[design$dist[k]]](data)
        list(converged = fit_converged(f), values = as.list(extractor(f)),
             error = NA_character_)
      }, error = function(e)
        list(converged = FALSE, values = list(), error = conditionMessage(e))),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    out$warning <- if (length(warn)) paste(unique(warn), collapse = "; ") else NA_character_
    out
  }

  cores <- usable_cores(cores)
  if (cores > 1) {
    results <- parallel::mclapply(seq_len(nrow(design)), run_one,
                                  mc.cores = cores, mc.preschedule = TRUE)
  } else {
    results <- lapply(seq_len(nrow(design)), run_one)
  }

  out <- cbind(design,
               converged = vapply(results, function(r) r$converged, logical(1)),
               error = vapply(results, function(r) r$error, character(1)),
               warning = vapply(results, function(r) r$warning, character(1)),
               stringsAsFactors = FALSE)
  values <- bind_values(lapply(results, function(r) r$values))
  if (ncol(values) > 0)
    out <- cbind(out, values)
  attr(out, "info") <- list(
    seed = seed, fit = if (is.function(fit)) "function" else fit,
    versions = package_versions(fit))
  out
}


## helpers

# Returns function(data) that fits `model` to one simulated data frame.
make_fitter <- function(fit, model, dist, ...) {
  if (is.function(fit)) {
    if (...length() > 0)
      warning("Arguments in `...` are ignored when `fit` is a function; ",
              "set such options inside the function.", call. = FALSE)
    return(function(data) fit(data, model))
  }
  fit <- match.arg(fit, c("lavaan", "magmaan"))
  if (fit == "lavaan")
    return(function(data) lavaan::sem(model, data = data, ...))

  if (!requireNamespace("magmaan", quietly = TRUE) ||
      utils::packageVersion("magmaan") < "0.2.0")
    stop("fit = \"magmaan\" requires magmaan (>= 0.2.0). ",
         "Use fit = \"lavaan\", or install a recent magmaan.")
  # construct the model once, with a zero-row prototype for the schema
  vars <- colnames(dist$sigma.target)
  if (is.null(vars))
    stop("The calibration's sigma.target needs column names matching the model.")
  prototype <- as.data.frame(matrix(numeric(0), 0, length(vars),
                                    dimnames = list(NULL, vars)))
  m <- magmaan::magmaan_model(model, prototype = prototype)
  function(data) magmaan::magmaan(m, data, ...)
}

# Simulated matrix -> data frame with the variable names of sigma.target
as_sim_data <- function(x, dist) {
  x <- as.data.frame(x)
  vars <- colnames(dist$sigma.target)
  if (!is.null(vars) && length(vars) == ncol(x))
    names(x) <- vars
  x
}

fit_converged <- function(f) {
  if (inherits(f, "lavaan"))
    return(isTRUE(lavaan::lavInspect(f, "converged")))
  if (inherits(f, "magmaan"))
    return(isTRUE(magmaan::as_lab_fit(f)$converged))
  NA
}

# List of named lists -> data frame, filling missing names with NA. Scalar
# entries become ordinary columns, others become list columns.
bind_values <- function(values) {
  keys <- unique(unlist(lapply(values, names)))
  cols <- lapply(keys, function(key) {
    col <- lapply(values, function(v) if (key %in% names(v)) v[[key]] else NA)
    scalar <- all(vapply(col, function(x) is.atomic(x) && length(x) == 1, logical(1)))
    if (scalar) unlist(col) else I(col)
  })
  names(cols) <- keys
  as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
}

package_versions <- function(fit) {
  pkgs <- c("covsim", if (identical(fit, "lavaan")) "lavaan",
            if (identical(fit, "magmaan")) c("magmaan", "magmaanlab"))
  vapply(pkgs, function(p) tryCatch(as.character(utils::packageVersion(p)),
                                    error = function(e) NA_character_),
         character(1))
}
