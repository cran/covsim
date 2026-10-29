# Number of worker processes for parallel::mclapply. Forking is not available
# on Windows, and CRAN checks (_R_CHECK_LIMIT_CORES_) allow at most two cores.
usable_cores <- function(cores) {
  cores <- suppressWarnings(as.integer(cores))
  if (length(cores) != 1 || is.na(cores) || cores < 1L)
    return(1L)
  if (.Platform$OS.type == "windows")
    return(1L)
  limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_"))
  if (nzchar(limit) && limit != "false")
    cores <- min(cores, 2L)
  cores
}

create.submatrix <- function(I, Matrix)
{
  # I <- unique(I)
  d <- dim(Matrix)[1]
  l <- length(I)

  for (i in (d:1))
  {
    if(sum(Matrix[, i] %in% I) == l)
      break
  }

  sub.matrix <- Matrix[(i:d), (i:d)]

  remove.ind <- which(!(diag(sub.matrix) %in% I))
  if (length(remove.ind) != 0)
  {
    sub.matrix <- sub.matrix[, -remove.ind]
  }

  new.sub.matrix <- NULL

  for (i in (1:(dim(sub.matrix)[2])))
  {
    indx <- which(sub.matrix[, i] %in% I)

    tmp.mat <- sub.matrix[indx, i]
    new.sub.matrix <- cbind(new.sub.matrix, c(rep(0, length(I) - length(tmp.mat)),
                                              tmp.mat))
  }

  sub.matrix <- new.sub.matrix
  return(list(sub.matrix = sub.matrix))
}

## functions that return the order of pcs in vines structure recursive!
get_uncond_pcs <- function(M)
{
  a <- cbind(diag(M[nrow(M):1, ]), M[1, ])
  return(a[-nrow(a), ])
}
get_pair_idx <- function(M)
{
  res <- matrix(get_uncond_pcs(M), ncol = 2)
  tree <- 1
  res <- cbind(res, tree)
  if (nrow(M) > 2)
  {
    for (tree in 2:(nrow(M) - 1))
    {
      tmpM <- M[tree:nrow(M), 1:length(tree:nrow(M))]
      tmpM <- matrix(get_uncond_pcs(tmpM), ncol = 2)
      tmpM <- cbind(tmpM, tree)
      res <- rbind(res, tmpM)
    }
  }
  res
}
# un-flatten
unflatten <- function(longlist)
{
  position <- length(longlist)
  ll <- 1
  mylist <- list()
  while (position > 0)
  {
    mylist <- c(mylist, list(longlist[(position - ll + 1):position]))
    position <- position - ll
    ll <- ll + 1
  }
  rev(mylist)
}

## adf test 2. We first estimate W Matrix from large sample, and use it
## repeatedly.
# test_covariance2 <- function(sigma.target, vc, n = 10^5, bign = 10^6, reps = 100)
# {
#   model <- NULL
#   for (i in 1:ncol(sigma.target)) model <- paste(model, paste0("x", i,
#                                                                "~~", diag(sigma.target)[i], "*x", i, "\n"))
#   for (i in 2:nrow(sigma.target))
#   {
#     for (j in 1:(i - 1)) model <- paste(model, paste0("x", i, "~~",
#                                                       sigma.target[i, j], "*x", j, "\n"))
#   }
#   ## large WLS weight matrix
#   bigsample <- data.frame(rvine(bign, vc, cores=parallel::detectCores()))
#   wls.V <- solve(lavaan:::lavGamma(bigsample)[[1]])
#   wls.0 <- lavaan::lav_matrix_vech(sigma.target)
#   pvalues <- numeric(reps)
#   set.seed(1)
#   pb <- txtProgressBar(min = 0, max = reps, style = 3)
#   for (i in 1:reps)
#   {
#     wls.obs <- lavaan::lav_matrix_vech(cov(rvinecopulib::rvine(n, vc,cores=parallel::detectCores())))
#
#     wls.res <- as.matrix(wls.obs - wls.0)
#     TestStat <- (n - 1) * as.numeric(t(wls.res) %*% wls.V %*% wls.res)
#     pvalues[i] <- 1 - pchisq(TestStat, df = length(wls.obs))
#     setTxtProgressBar(pb, i)
#   }
#   # pvalues should be uniformly distributed
#   res <- ks.test(pvalues, y = "punif")
#   res$p.value
# }



# one-parameter copulas
get_lowerupper <- function(name)
{
  names <- c("gaussian", "clayton", "gumbel", "frank", "joe", "indep")
  bound_list <- list(c(-1, 1), c(0.0001, 28), c(1, 50), c(-35, 34.8), c(1,
                                                                      30), c(0,0))
  bound_list[[grep(name, names, fixed = T)]]
}

#####
## PLSIM functions
####
get_truncated_moments <- function(g1,g2){# orjebins recursive formula
  # second term in recursive formula should be zero at infinite endpoint
  infg1 <- is.infinite(g1); infg2 <-  is.infinite(g2)
  Diff <- stats::pnorm(g2)-stats::pnorm(g1)
  m <- c(0,1, (stats::dnorm(g1)-stats::dnorm(g2))/Diff)
  m <- c(m, (length(m)-2)*m[length(m)-1]-(ifelse(infg2,0, g2)^(length(m)-2)*stats::dnorm(g2)-ifelse(infg1,0, g1)^(length(m)-2)*stats::dnorm(g1))/Diff)
  m <- c(m, (length(m)-2)*m[length(m)-1]-(ifelse(infg2,0, g2)^(length(m)-2)*stats::dnorm(g2)-ifelse(infg1,0, g1)^(length(m)-2)*stats::dnorm(g1))/Diff)
  m <- c(m, (length(m)-2)*m[length(m)-1]-(ifelse(infg2,0, g2)^(length(m)-2)*stats::dnorm(g2)-ifelse(infg1,0, g1)^(length(m)-2)*stats::dnorm(g1))/Diff)
  m[3:length(m)]
}

# (a_iZ+b_i)*I( g1 < Z < g2)
get_segment_moments <- function(a, b, g1, g2){
  #moments of Z*I( g1 < Z < g2)
  mom <- c(1, get_truncated_moments(g1, g2))*(stats::pnorm(g2)-stats::pnorm(g1))

  koefs <- matrix(c(0, 0, 0, a, b,
                    0, 0, a^2, 2*a*b, b^2,
                    0, a^3, 3*a^2*b, 3*a*b^2, b^3,
                    a^4, 4*a^3*b, 6*a^2*b^2,4*a*b^3, b^4), byrow = T, nrow=4)
  as.vector(koefs %*% matrix(rev(mom)))
}



#skew and kurtosis
skewkurt_discrepancy <- function(a, gamma, skew, kurt){#
  ##  continuity.zero mean.
  b <- get_bs(a, gamma)
  mom <- get_pl_moments(a,b, gamma)
  var <- mom[2]
  s <- mom[3]/var^1.5
  k <- mom[4]/var^2-3
  (s-skew)^2+(k-kurt)^2
}
#calibrate b for zero mean and continuity
get_bs <- function(a, gamma){
  b <- rep(0, length(a))
  for(i in 2:length(b))
    b[i] <- b[i-1]+(a[i-1]-a[i])*gamma[i-1]
  b-get_pl_mean(a,b,gamma)
}

get_pl_mean <- function(a,b,gamma){
  ddiff <- diff(stats::dnorm(c(-Inf, gamma, Inf)))
  pdiff <- diff(stats::pnorm(c(-Inf, gamma, Inf)))
  -sum(ddiff*a)+sum(pdiff*b)#in paper
}

get_pl_moments <- function(a, b, gamma){
  g <- c(-Inf, gamma, Inf)
  res <- sapply(1:length(a), function(i){
    get_segment_moments(a[i],b[i], g[i], g[i+1])
  })
  rowSums(res)
}

## fit piecewise linear functions
fit_univariate <- function(gammalist, skew, kurt, scale=FALSE, monot){

  alist <- lapply( 1:length(skew), function(i){
    out <- stats::nlminb(start=rep(1, length(gammalist[[i]])+1),
                  objective =skewkurt_discrepancy,
                  gamma=gammalist[[i]], skew=skew[i], kurt=kurt[i], lower=ifelse(monot, 0.05,-Inf))
    if(out$objective< 1e-4)
      return(out$par)
    else
      return(NA)
  })

  if(sum(is.na(alist)) > 0)
     return(NA)

  if(scale){
    alist <- lapply(1:length(alist), function(i){
      a <- alist[[i]]; gamma <- gammalist[[i]]
      b <- get_bs(a, gamma)
      mom <- get_pl_moments(a, b, gamma)
      a/sqrt(mom[2])
    })
  }
  return(alist)
}


# regular gammas
get_gamma <- function(k) stats::qnorm((1:k)/k)[1:(k-1)]

#
pl_fun <- function(x, a, b, gamma){
  k <- findInterval(x, gamma)+1
  a[k]*x+b[k]
}


## COVARIANCE
# Covariance of two continuous piecewise linear transforms f, g of a standard
# bivariate normal pair with correlation rho. By Price's theorem,
# d/drho E[f(Z1)g(Z2)] = E[f'(Z1)g'(Z2)], and since f'' and g'' are point masses
# at the breakpoints (slope jumps diff(a)),
#   Cov = rho*E[f']E[g'] + sum_kl da1_k da2_l int_0^rho (rho-t) phi2(gk, gl; t) dt.
# The integral is smooth after substituting t = sin(theta), and Gauss-Legendre
# quadrature in theta is accurate to machine precision even for |rho| near 1.
# b1, b2 are unused (the covariance depends on the slopes and breakpoints only)
# but kept so the signature matches the old tmvtnorm-based implementation.
get_cov <- function(a1, b1, gamma1, a2, b2, gamma2, rho){
  p1 <- diff(stats::pnorm(c(-Inf, gamma1, Inf)))
  p2 <- diff(stats::pnorm(c(-Inf, gamma2, Inf)))
  cov <- rho * sum(a1 * p1) * sum(a2 * p2)
  if (rho == 0 || length(gamma1) == 0 || length(gamma2) == 0)
    return(cov)

  gl <- gauss_legendre()
  theta <- asin(rho) * (gl$nodes + 1) / 2
  w <- gl$weights * asin(rho) / 2 * cos(theta) * (rho - sin(theta))
  t <- sin(theta)

  x <- rep(gamma1, times = length(gamma2))
  y <- rep(gamma2, each = length(gamma1))
  jumps <- as.vector(outer(diff(a1), diff(a2)))
  s2 <- 1 - t^2
  q <- outer(x^2 + y^2, rep(1, length(t))) - 2 * outer(x * y, t)
  phi2 <- exp(-q / rep(2 * s2, each = length(x))) /
    rep(2 * pi * sqrt(s2), each = length(x))
  cov + sum(jumps * as.vector(phi2 %*% w))
}

# Intermediate normal correlation giving correlation `target` after the
# (unit-variance) piecewise linear transforms. Cov(0) = 0, so for a positive
# target we bracket on [0, bound] and for a negative one on [-bound, 0]. If that
# fails (possible for non-monotone transforms) a grid over [-bound, bound] is
# scanned for a sign change. Returns NA if no root exists.
solve_pl_corr <- function(a1, b1, gamma1, a2, b2, gamma2, target, bound = 0.999){
  if (target == 0)
    return(0)
  f <- function(rho) get_cov(a1, b1, gamma1, a2, b2, gamma2, rho) - target
  if (target > 0) {
    lims <- c(0, bound); f.lims <- c(-target, f(bound))
  } else {
    lims <- c(-bound, 0); f.lims <- c(f(-bound), -target)
  }
  if (prod(f.lims) <= 0)
    return(stats::uniroot(f, lims, f.lower = f.lims[1], f.upper = f.lims[2],
                          tol = 1e-12)$root)
  grid <- seq(-bound, bound, length.out = 41)
  f.grid <- vapply(grid, f, numeric(1))
  k <- which(f.grid[-1] * f.grid[-length(grid)] <= 0)
  if (length(k) == 0)
    return(NA_real_)
  k <- k[which.min(abs(grid[k] - target))]
  stats::uniroot(f, grid[c(k, k + 1)], f.lower = f.grid[k],
                 f.upper = f.grid[k + 1], tol = 1e-12)$root
}

# Attainable correlation range for a pair of unit-variance PL transforms
pl_corr_range <- function(a1, b1, gamma1, a2, b2, gamma2, bound = 0.999){
  grid <- seq(-bound, bound, length.out = 41)
  r <- vapply(grid, function(rho) get_cov(a1, b1, gamma1, a2, b2, gamma2, rho),
              numeric(1))
  range(r)
}

# 32-point Gauss-Legendre rule on [-1, 1] (Golub-Welsch), computed once
gauss_legendre <- local({
  rule <- NULL
  function() {
    if (is.null(rule)) {
      n <- 32
      i <- seq_len(n - 1)
      J <- matrix(0, n, n)
      J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- i / sqrt(4 * i^2 - 1)
      e <- eigen(J, symmetric = TRUE)
      rule <<- list(nodes = e$values, weights = 2 * e$vectors[1, ]^2)
    }
    rule
  }
})






