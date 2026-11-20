test_that("a factor's dummy columns form a single covariate group", {
  set.seed(1)
  site <- data.frame(
    elev = rnorm(40),
    habitat = factor(sample(c("a", "b", "c"), 40, replace = TRUE),
                     levels = c("a", "b", "c"))
  )
  d <- bvsDesign(~ elev + habitat, site)

  expect_s3_class(d, "bvsDesign")
  expect_equal(d$numVars, 4L)          # intercept, elev, two dummies
  expect_equal(d$ncov, 2L)             # elev and habitat
  expect_equal(d$indexes_covariates, c(1L, 2L, 3L, 3L))
  expect_equal(d$labels, c("(Intercept)", "elev", "habitat"))
})

test_that("an interaction is treated as its own covariate group", {
  set.seed(2)
  site <- data.frame(a = rnorm(30), b = factor(sample(1:2, 30, TRUE)))
  d <- bvsDesign(~ a * b, site)
  # intercept, a, b2, a:b2 -> groups 1, 2, 3, 4
  expect_equal(d$indexes_covariates, c(1L, 2L, 3L, 4L))
  expect_equal(d$ncov, 3L)
})

test_that("a formula without an intercept is rejected", {
  site <- data.frame(a = rnorm(10))
  expect_error(bvsDesign(~ a - 1, site), "intercept")
})

test_that("missing covariate values are rejected", {
  site <- data.frame(a = c(rnorm(9), NA))
  expect_error(bvsDesign(~ a, site), "missing values")
})

test_that("array designs put the covariate index in the second dimension", {
  set.seed(3)
  M <- 12; J <- 4
  visits <- data.frame(
    wind = rnorm(M * J),
    obs = factor(sample(c("x", "y"), M * J, replace = TRUE), levels = c("x", "y"))
  )
  d <- bvsDesignArray(~ wind + obs, visits, dims = c(M, J))

  expect_equal(dim(d$X), c(M, 3L, J))
  expect_equal(d$indexes_covariates, c(1L, 2L, 3L))

  # The slice for visit j must equal the model matrix rows for that visit.
  mm <- model.matrix(~ wind + obs, visits)
  for (j in seq_len(J)) {
    rows <- (j - 1L) * M + seq_len(M)
    expect_equal(d$X[, , j], mm[rows, ], ignore_attr = TRUE)
  }
})

test_that("array designs check the row count", {
  visits <- data.frame(wind = rnorm(10))
  expect_error(bvsDesignArray(~ wind, visits, dims = c(4, 4)), "rows")
})

test_that("four-dimensional designs are laid out correctly", {
  set.seed(4)
  M <- 6; J <- 3; Tn <- 2
  dat <- data.frame(w = rnorm(M * J * Tn))
  d <- bvsDesignArray(~ w, dat, dims = c(M, J, Tn))
  expect_equal(dim(d$X), c(M, 2L, J, Tn))
  # element (i, 2, j, t) is the covariate for row i + (j-1)M + (t-1)MJ
  i <- 2; j <- 3; t <- 2
  expect_equal(d$X[i, 2, j, t], dat$w[i + (j - 1) * M + (t - 1) * M * J])
})

test_that("the prior gives a factor's columns an exchangeable block", {
  ic <- c(1L, 2L, 3L, 3L)
  p <- bvsPrior(as.bvsDesign(matrix(0, 5, 4), ic), phi.mu = 4, phi.beta = 0.25)

  expect_equal(p$b, c(0.5, 0, 0, 0))
  expect_equal(p$B[1, 1], 4)
  expect_equal(p$B[2, 2], 0.25)                  # continuous covariate
  expect_equal(p$B[3, 3], 0.25)                  # factor variance
  expect_equal(p$B[3, 4], 0.125)                 # correlation of one half
  expect_equal(p$B[2, 3], 0)                     # no correlation across groups
  expect_true(all(eigen(p$B, only.values = TRUE)$values > 0))
})

test_that("malformed covariate groupings are rejected", {
  X <- matrix(0, 5, 3)
  expect_error(as.bvsDesign(X, c(2L, 3L, 4L)), "start at 1")
  expect_error(as.bvsDesign(X, c(1L, 3L, 2L)), "non-decreasing")
  expect_error(as.bvsDesign(X, c(1L, 3L, 3L)), "no gaps")
  expect_error(as.bvsDesign(X, c(1L, 1L, 2L)), "only the intercept")
  expect_error(as.bvsDesign(X, c(1L, 2L)), "length 2")
})

test_that("indexes_covariates is stored as double, not integer", {
  ## The nimbleFunctions that consume this vector (dLgamma, compute_predictor,
  ## build_block_cov) declare it as double(1). An integer vector compiles to
  ## NimArr<1,int>, the generated C++ then fails to match the double signature,
  ## and the model builds fine in R but dies at compileNimble(). Keep it double.
  d <- bvsDesign(~ x1 + f, data.frame(x1 = rnorm(20),
                                      f = factor(rep(c("a", "b", "c"), length.out = 20))))
  expect_type(d$indexes_covariates, "double")
  expect_false(is.integer(d$indexes_covariates))
  expect_equal(d$indexes_covariates, c(1, 2, 3, 3))

  a <- bvsDesignArray(~ x1, data.frame(x1 = rnorm(20 * 3)), dims = c(20, 3))
  expect_type(a$indexes_covariates, "double")

  ## also when the user hands in integers directly
  X <- cbind(1, matrix(rnorm(40), 20, 2))
  b <- as.bvsDesign(X, indexes_covariates = c(1L, 2L, 3L))
  expect_type(b$indexes_covariates, "double")

  ## group counts stay integer
  expect_type(d$ncov, "integer")
  expect_type(b$ncov, "integer")
})
