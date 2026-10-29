test_that("IG works", {
  hsmodel <- "visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6"
  f <- lavaan::cfa(hsmodel, lavaan::HolzingerSwineford1939)
  sigma0 <- lavaan::lavInspect(f, "sigma.hat")
  target0 <- lavaan::lav_matrix_vech(sigma0)

  s <- rep(2, 6)
  k <- rep(7, 6)
  set.seed(1)
  sample1 <- rIG(10^6,sigma0, s, k, reps=1, typeA="triang")[[1]]
  sample2 <- rIG(10^6,sigma0, s, k, reps=1, typeA="symm")[[1]]

  cov1 <- lavaan::lav_matrix_vech(cov(sample1))
  cov2 <- lavaan::lav_matrix_vech(cov(sample2))
  skew1 <- psych::skew(sample1)
  skew2 <- psych::skew(sample2)
  kurt1 <- unname(psych::kurtosi(sample1))
  kurt2 <- unname(psych::kurtosi(sample2))

  expect_equal(cov1, target0, tolerance=0.01)
  expect_equal(cov2, target0, tolerance=0.01)
  expect_equal(skew1, s, tolerance=0.1)
  expect_equal(skew2, s, tolerance=0.1)
  expect_equal(kurt1, k, tolerance=0.1)
  expect_equal(kurt2, k, tolerance=0.1)

})




test_that("IG generator moments solve the linear moment systems exactly", {
  set.seed(1)
  sigma0 <- cov(MASS::mvrnorm(20, rep(0, 4), diag(4)))
  s <- c(0.5, 1, 1.5, 2); k <- c(2, 4, 6, 8)
  for (typeA in c("triang", "symm")) {
    cal <- calibrate_ig(sigma0, s, k, typeA = typeA)
    A <- cal$A; sc <- sqrt(rowSums(A^2))
    expect_equal(as.vector((A^3 / sc^3) %*% cal$IGskew), s, tolerance = 1e-10)
    expect_equal(as.vector((A^4 / sc^4) %*% cal$IGkurt.excess), k, tolerance = 1e-10)
  }
})


test_that("IG reports infeasible generator moments", {
  sigma0 <- matrix(c(1, 0.9, 0.9, 1), 2)
  expect_error(calibrate_ig(sigma0, c(0, 3), c(0, 5)), "No valid IG generator")
})
