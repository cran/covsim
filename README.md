covsim
==========

This package contains  implementations of four simulation procedures
that simulate from a multivariate distributions with a pre-specified
covariance matrix. In addition, the user may pre-specify marginal information. The main function is vita, which implements the very general VIne-To-Anything (VITA) algorithm of Grønneberg and Foldnes (Psychometrika, 2017).


How to install
--------------


You can install:

-   the stable release on CRAN:

    ``` r
    install.packages("covsim")
    ```

-   the latest development version:

    ``` r
    devtools::install_github("njaalf/covsim")
    ```

------------------------------------------------------------------------

The calibrate / simulate workflow (since 1.2.0)
------------------------------------------------

Each of the four algorithms is now split into a **calibration** step (which
can be slow but only needs to be run once for a given target) and a
**simulation** step (which is cheap and can be repeated as many times as you
like). The calibration object is reusable: calibrate once, draw many
independent samples without rerunning the expensive part.

``` r
# Step 1 — calibrate once
cal <- calibrate_ig(sigma.target, skewness, excesskurtosis)
# or: calibrate_plsim(sigma.target, skewness, excesskurtosis)
# or: calibrate_vm(sigma.target, skewness, excesskurtosis)
# or: calibrate_vita(margins, sigma.target)

# Step 2 — simulate as many datasets as you want
samples <- simulate(cal, nsim = 100, N = 1000)   # list of 100 data matrices
```

The old `rIG()`, `rPLSIM()` and `vita()` functions are unchanged in signature
and return shape; they are now thin wrappers over the new API. `rVM()` is the
corresponding wrapper for VM.

### From calibration to a Monte Carlo study

`simulate_estimate()` draws samples from one or more calibrations, fits a
model to each sample and collects the results in one data frame:

``` r
dists <- list(IG = calibrate_ig(sigma, skew, kurt),
              VM = calibrate_vm(sigma, skew, kurt))
res <- simulate_estimate(dists, model = "F =~ x1 + x2 + x3 + x4",
                         n = c(100, 500), reps = 1000, cores = 4,
                         extractor = function(fit) lavaan::fitMeasures(fit, "pvalue"))
aggregate(pvalue ~ dist + n, res, function(p) mean(p < 0.05))
```

Estimation uses lavaan by default. With `fit = "magmaan"` the much faster
[magmaan](https://github.com/JonasMoss/magmaan) estimator is used if it is
installed, and any `function(data, model)` can be supplied instead.

See `vignette("calibrate-simulate", package = "covsim")` for a tutorial and
`vignette("simulation-study", package = "covsim")` for a complete Monte Carlo
study with lavaan and magmaan.

Package overview
----------------



### VITA

The user must completely specify the marginals from the distributions available in the stats package. Also the user must specify the population covariance matrix. Note that the variances in this matrix must equal the population variances of the marginal distributions. 
A regular vine may also be provided, that is, a list of bicopulas and a sequence for their calibrations, i.e. a sequence of trees. If no vine
is provided, a default D-vine is assumed. Note that in many cases the vine, the marginals and the covariance matrix are not compatible. That is, no vine distribution exists whose covariance equals the target covariance. Also note that calibration of the vine is computationally burdensome. To speed up calibration, one could reduce the numpoints and/or the Nmax parameters. This would reduce the precision of the adherence of the vita vine to the target population covariance matrix.

### IG

In addition to vita, a simpler and much more restricted simulation algorithm is provided in the rIG function, which implements the independent generator approach of Foldnes and Olsson (MBR, 2016)
The IG approach is comparable to the Vale-Maurelli approach, in that
it is fast and precise. Only skewness and kurtosis of the marginals is provided. The marginal distributions are therefore not fully controlled using IG. 

### PLSIM 

has a similar interface as that of IG. Skewness and kurtosis is pre-specified. PLSIM is flexible beyond the default values. Numsegments can be increased, or the breakpoints manually set by providing gammalist.  

### VM

The Vale-Maurelli method (Vale and Maurelli, 1983) transforms correlated normal variables with Fleishman polynomials. Like IG and PLSIM, skewness and excess kurtosis are pre-specified.


References
----------
Grønneberg, S., Foldnes, N. Covariance Model Simulation Using Regular Vines. Psychometrika 82, 1035–1051 (2017). https://doi.org/10.1007/s11336-017-9569-6
 
Grønneberg, S., Foldnes, N., & Marcoulides, K. M. (2022). covsim: An r package for simulating non-normal data for structural equation models using copulas. Journal of Statistical Software. DOI: 10.18637/jss.v102.i03

Njål Foldnes & Ulf Henning Olsson (2016) A Simple Simulation Technique for Nonnormal Data with Prespecified Skewness, Kurtosis, and Covariance Matrix, Multivariate Behavioral Research, 51:2-3, 207-219, DOI: 10.1080/00273171.2015.1133274

Foldnes, N., & Grønneberg, S. (2022). Non-normal data simulation using piecewise linear transforms. Structural Equation Modeling: A Multidisciplinary Journal, 29(1), 36-46.

Vale, C. D., & Maurelli, V. A. (1983). Simulating multivariate nonnormal distributions. Psychometrika, 48(3), 465-471.

