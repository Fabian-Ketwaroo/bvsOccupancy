test_that("node names are parsed into an index matrix", {
  ix <- bvsOccupancy:::.parseNodeIndices(c("y[1, 2, 3]", "y[10, 1, 4]"))
  expect_equal(ix, matrix(c(1L, 10L, 2L, 1L, 3L, 4L), ncol = 3))

  ix1 <- bvsOccupancy:::.parseNodeIndices(c("z[1]", "z[7]"))
  expect_equal(ix1, matrix(c(1L, 7L), ncol = 1))

  expect_error(bvsOccupancy:::.parseNodeIndices(c("y[1, 2]", "y[3]")), "inconsistent")
})

test_that("a design array is flattened to match a given node order", {
  set.seed(5)
  M <- 4; V <- 3; J <- 2
  X <- array(rnorm(M * V * J), dim = c(M, V, J))

  # rows in an arbitrary order, as the model might report them
  idx <- cbind(c(3L, 1L, 4L), c(2L, 1L, 2L))   # (site, visit)
  flat <- bvsOccupancy:::.flattenDesign(X, idx, map = c(1L, 2L))

  expect_equal(dim(flat), c(3L, V))
  for (r in seq_len(nrow(idx))) {
    expect_equal(flat[r, ], X[idx[r, 1], , idx[r, 2]])
  }
})

test_that("a design slice keeps its matrix shape", {
  X <- array(seq_len(4 * 3 * 2), dim = c(4, 3, 2))
  s <- bvsOccupancy:::.slice(X, 2)
  expect_equal(dim(s), c(4L, 3L))
  expect_equal(s, X[, , 2])
})

test_that("surveyed visits are moved to the front and covariates follow", {
  y <- matrix(c(1, NA, 0,
                NA, 1, 1), nrow = 2, byrow = TRUE)
  X <- array(0, dim = c(2, 1, 3))
  X[1, 1, ] <- c(10, 20, 30)
  X[2, 1, ] <- c(40, 50, 60)

  out <- bvsOccupancy:::.leftPack(y, X)

  expect_true(out$changed)
  expect_equal(as.vector(out$y[1, ]), c(1, 0, NA))
  expect_equal(as.vector(out$y[2, ]), c(1, 1, NA))
  # covariates travel with their visit
  expect_equal(as.vector(out$X[1, 1, ]), c(10, 30, 20))
  expect_equal(as.vector(out$X[2, 1, ]), c(50, 60, 40))
})

test_that("left packing leaves complete data alone", {
  y <- matrix(c(1, 0, 1, 1), 2, 2)
  X <- array(rnorm(2 * 2 * 2), dim = c(2, 2, 2))
  out <- bvsOccupancy:::.leftPack(y, X)
  expect_false(out$changed)
  expect_equal(out$y, y)
  expect_equal(out$X, X)
})

test_that("Polya-Gamma starting values have the right mean", {
  set.seed(6)
  draws <- bvsOccupancy:::.rpg0(4000)
  expect_true(all(draws > 0))
  expect_equal(mean(draws), 0.25, tolerance = 0.02)   # PG(1, 0) has mean 1/4
})

test_that("detection histories must be binary", {
  expect_error(bvsOccupancy:::.checkY(matrix(c(0, 1, 2, 1), 2, 2), 2L), "only 0, 1 and NA")
  expect_error(bvsOccupancy:::.checkY(array(0, c(2, 2, 2)), 2L), "2-dimensional")
})

test_that("MCMC settings are validated", {
  expect_error(bvsOccupancy:::.mcmcControl(100, 100, 1, 2, TRUE, TRUE, NULL, FALSE), "burnin")
  expect_error(bvsOccupancy:::.mcmcControl(100, 0, 0, 2, TRUE, TRUE, NULL, FALSE), "thin")
  expect_error(bvsOccupancy:::.mcmcControl(100, 90, 50, 2, TRUE, TRUE, NULL, FALSE), "fewer than 2")
  ok <- bvsOccupancy:::.mcmcControl(1000, 500, 5, 2, TRUE, TRUE, 1, FALSE)
  expect_equal(ok$n.kept, 100)
})

test_that("chains after the first are jittered", {
  set.seed(7)
  base <- list(beta_psi = rep(0, 3), gamma_psi = rep(1, 3), Omega_psi = rep(0.25, 5))
  ini <- bvsOccupancy:::.makeInits(base, n.chains = 3)

  expect_length(ini, 3)
  expect_equal(ini[[1]]$beta_psi, rep(0, 3))       # first chain is the base
  expect_false(isTRUE(all.equal(ini[[2]]$beta_psi, rep(0, 3))))
  expect_equal(ini[[2]]$gamma_psi, rep(1, 3))      # indicators are not jittered
})

test_that("the potential scale reduction factor is sane", {
  set.seed(8)
  same <- matrix(rnorm(2000), ncol = 2)
  expect_lt(bvsOccupancy:::.rhat(same), 1.05)

  split <- cbind(rnorm(1000), rnorm(1000, mean = 5))
  expect_gt(bvsOccupancy:::.rhat(split), 1.5)

  expect_true(is.na(bvsOccupancy:::.rhat(matrix(1, 100, 2))))  # a stuck indicator
})
