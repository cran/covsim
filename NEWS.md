# covsim 1.2.0

New `simulate_estimate()` runs a complete Monte Carlo study: for one or more
calibrated distributions and sample sizes it draws `reps` samples, fits a
model to each, and returns one row per replication with convergence, errors,
warnings and the values returned by a user-supplied `extractor`. Fitting uses
lavaan by default, magmaan (`fit = "magmaan"`, optional, much faster) or any
user function. Every replication has its own L'Ecuyer-CMRG stream, so results
are identical for any number of `cores`. A new vignette,
`vignette("simulation-study")`, walks through a complete study with lavaan
and magmaan.

New VM (Vale-Maurelli) simulation: `calibrate_vm()` with a `simulate()`
method, and the `rVM()` convenience wrapper. Fleishman coefficients are
solved from the exact moment equations (preferring a monotone polynomial when
several solutions exist) and intermediate correlations from the
Vale-Maurelli cubic.

PLSIM samples now carry the column names of `sigma.target`, as IG samples
already did.

PLSIM calibration is 40-100x faster. The covariance of two piecewise linear
transforms is now computed exactly from Price's theorem with a
one-dimensional Gauss-Legendre quadrature, replacing the truncated bivariate
normal moments from `tmvtnorm` (no longer a dependency). Intermediate
correlations are found by bracketed `uniroot()` to 1e-12 instead of
`nlminb()` on a squared discrepancy, and an infeasible target now reports the
attainable correlation range. The calibration object also stores the
smallest eigenvalue of the intermediate correlation matrix (`min.eigen`).

Seeded output of `rIG()` and `rPLSIM()` agrees with covsim 1.1.0 up to the
small numerical differences of the exact solvers (about 1e-7). Occasionally
such a difference flips the sign of an eigenvector of the intermediate
correlation matrix inside `MASS::mvrnorm()`, giving different draws from the
same distribution, as can also happen across platforms.

IG: the linear generator skewness and kurtosis systems are solved directly
with `solve()` instead of iteratively with `nleqslv()`, whose convergence was
not checked. The solution is the same up to numerical tolerance. Infeasible
generator moments are reported per variable before Pearson fitting, and the
generator moments are stored in the calibration object (`IGskew`,
`IGkurt.excess`).

New calibrate/simulate API: `calibrate_ig()`, `calibrate_plsim()`, and
`calibrate_vita()` return reusable calibration objects. Sampling is done via
S3 methods on `stats::simulate()`, so a single (potentially expensive)
calibration can drive many simulations without recomputation. The existing
`rIG()`, `rPLSIM()`, and `vita()` functions are unchanged in signature and
return shape; they now delegate to the new API.

`vita()` is unchanged from covsim 1.1.0 and gives identical seeded output.
A faster deterministic root search was tried during development, but a
matched-accuracy study (`validation/` on the GitHub repository) showed no
consistent gain over the 1.1.0 search, so it was not adopted.

On Windows, where forking is unavailable, `simulate_estimate()` runs on one
core, and under CRAN's core limit (`_R_CHECK_LIMIT_CORES_`) on at most two.

# covsim 1.1.0

rIG: The A matrix is now triangular by default

rPLSIM: Added tests for whether an identical pair has been calibrated. 
This speeds up rPLSIM if sigma.target has many identical values, and skewness and kurtosis values are identical. 

vita: Now checks margins are correctly named. 

# covsim 1.0.0

Third release. Accompanies journal publication in JSS.
