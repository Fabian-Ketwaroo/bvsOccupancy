## End-to-end tests. Everything that needs a compiled model is skipped unless
## nimble is actually installed, so the suite still runs on a bare machine.

test_that("simulators return objects of the advertised shape", {
  sim <- simSSOM(M = 60, J = 3, n.cont = 2, n.cat = 1, n.levels = 3, seed = 1)
  expect_equal(dim(sim$y), c(60L, 3L))
  expect_true(all(sim$y %in% c(0L, 1L)))
  expect_s3_class(sim$occ.design, "bvsDesign")
  expect_s3_class(sim$det.design, "bvsDesign")
  expect_equal(length(sim$beta.psi), ncol(sim$occ.design$X))
  expect_equal(length(sim$beta.p), dim(sim$det.design$X)[2L])

  sdm <- simSDOM(M = 40, J = 3, n.seasons = 4, n.cont = 2, n.cat = 1, seed = 1)
  expect_equal(dim(sdm$y), c(40L, 3L, 4L))
  expect_equal(dim(sdm$z), c(40L, 4L))
  expect_equal(dim(sdm$per.design$X)[c(1L, 3L)], c(40L, 3L))

  msm <- simMSOM(N = 4, S = 50, nvisits = 3, n.cont = 2, n.cat = 1, seed = 1)
  expect_equal(dim(msm$y), c(4L, 50L, 3L))
  expect_equal(dim(msm$beta.z), c(4L, ncol(msm$occ.design$X)))
})

test_that("simulators honour unequal replication", {
  sim <- simSSOM(M = 20, J = 4, seed = 2, Jindex = rep(c(2L, 4L), 10))
  expect_equal(as.integer(rowSums(!is.na(sim$y))), rep(c(2L, 4L), 10))

  sdm <- simSDOM(M = 10, J = 4, n.seasons = 3, seed = 2,
                 Jindex = matrix(2L, 10, 3))
  expect_true(all(is.na(sdm$y[, 3:4, ])))
})

test_that("supplied coefficients are used verbatim", {
  d <- bvsDesign(~ x1 + x2, data.frame(x1 = rnorm(30), x2 = rnorm(30)))
  b <- c(0.5, -1, 0)
  sim <- simSSOM(M = 30, J = 2, n.cont = 2, n.cat = 0, beta.psi = b, seed = 3)
  expect_equal(sim$beta.psi, b)
  expect_equal(length(sim$psi), 30L)
})

test_that("models refuse malformed input before touching nimble", {
  sim <- simSSOM(M = 30, J = 3, seed = 4)

  expect_error(bvsSSOM(sim$y[1:10, ], sim$occ.design, sim$det.design, run = FALSE),
               "sites")
  expect_error(bvsSSOM(sim$y, sim$det.design, sim$det.design, run = FALSE),
               "2 dimensions")
  expect_error(bvsSSOM(sim$y, sim$occ.design, sim$occ.design, run = FALSE),
               "3 dimensions")

  y2 <- sim$y
  y2[1, 1] <- 3
  expect_error(bvsSSOM(y2, sim$occ.design, sim$det.design, run = FALSE),
               "0|1|NA")

  sdm <- simSDOM(M = 20, J = 3, n.seasons = 2, seed = 4)
  expect_error(
    bvsSDOM(sdm$y[, , 1, drop = FALSE], sdm$occ.design, sdm$per.design,
            sdm$col.design, sdm$det.design, run = FALSE),
    "two seasons"
  )

  msm <- simMSOM(N = 3, S = 30, nvisits = 3, seed = 4)
  expect_error(
    bvsMSOM(msm$y[1, , , drop = FALSE], msm$occ.design, msm$det.design, run = FALSE),
    "two species"
  )
  expect_error(
    bvsMSOM(msm$y, msm$occ.design, msm$det.design,
            tau.occ = list(intercept = c(10, 10)), run = FALSE),
    "covariate"
  )
})

test_that("priors are checked against the design they belong to", {
  d <- bvsDesign(~ x1 + x2, data.frame(x1 = rnorm(20), x2 = rnorm(20)))
  sim <- simSSOM(M = 20, J = 2, seed = 5)
  expect_error(
    bvsSSOM(sim$y, sim$occ.design, sim$det.design, prior.occ = bvsPrior(d),
            run = FALSE),
    "prior"
  )
})

## ---------------------------------------------------------------------------
## Everything below needs nimble.
## ---------------------------------------------------------------------------

test_that("bvsSSOM builds an uncompiled model", {
  skip_without_nimble()
  sim <- simSSOM(M = 40, J = 3, n.cont = 2, n.cat = 1, seed = 10)
  fit <- bvsSSOM(sim$y, sim$occ.design, sim$det.design,
                 run = FALSE, verbose = FALSE)
  expect_s3_class(fit, "bvsOccupancy")
  expect_equal(fit$type, "SSOM")
  expect_true(is.null(fit$samples))
})

test_that("bvsSSOM runs a short chain and returns usable output", {
  skip_on_cran()
  skip_without_nimble()
  sim <- simSSOM(M = 60, J = 3, n.cont = 2, n.cat = 1, seed = 11)
  fit <- bvsSSOM(sim$y, sim$occ.design, sim$det.design,
                 n.iter = 400, n.burnin = 100, n.thin = 1, n.chains = 2,
                 seed = 11, verbose = FALSE)

  expect_s3_class(fit, "bvsOccupancy")
  expect_s3_class(as.mcmc.list(fit), "mcmc.list")
  expect_equal(coda::nchain(as.mcmc.list(fit)), 2L)

  ip <- inclusionProbs(fit)
  expect_named(ip, c("psi", "p"))
  expect_true(all(ip$psi >= 0 & ip$psi <= 1))
  expect_equal(length(ip$psi), length(sim$occ.design$labels))
  expect_equal(unname(ip$psi[1L]), 1)
  expect_type(inclusionProbs(fit, "p"), "double")

  mm <- medianModel(fit)
  expect_type(mm, "list")
  expect_true(is.character(mm$psi))

  expect_output(print(fit), "single-season")
  s <- summary(fit)
  expect_true(is.list(s))
})

test_that("bvsSDOM runs a short chain", {
  skip_on_cran()
  skip_without_nimble()
  sim <- simSDOM(M = 40, J = 3, n.seasons = 3, n.cont = 2, n.cat = 1, seed = 12)
  fit <- bvsSDOM(sim$y, sim$occ.design, sim$per.design, sim$col.design,
                 sim$det.design,
                 n.iter = 400, n.burnin = 100, n.thin = 1, n.chains = 1,
                 seed = 12, verbose = FALSE)
  expect_equal(fit$type, "SDOM")
  ip <- inclusionProbs(fit)
  expect_named(ip, c("psi", "phi", "eta", "p"))
  expect_true(all(vapply(ip, function(v) all(v >= 0 & v <= 1), logical(1))))
})

test_that("bvsMSOM runs a short chain", {
  skip_on_cran()
  skip_without_nimble()
  sim <- simMSOM(N = 4, S = 50, nvisits = 3, n.cont = 2, n.cat = 1, seed = 13)
  fit <- bvsMSOM(sim$y, sim$occ.design, sim$det.design,
                 n.iter = 400, n.burnin = 100, n.thin = 1, n.chains = 1,
                 seed = 13, verbose = FALSE)
  expect_equal(fit$type, "MSOM")
  ip <- inclusionProbs(fit)
  expect_named(ip, c("z", "p"))
  expect_true(is.matrix(ip$z))
  expect_equal(nrow(ip$z), 4L)
  expect_true(all(ip$z >= 0 & ip$z <= 1))
})

test_that("a long SSOM run recovers the generating model", {
  skip_on_cran()
  skip_without_nimble()
  skip_if(Sys.getenv("BVSOCCUPANCY_SLOW_TESTS") == "", "slow test")

  sim <- simSSOM(M = 400, J = 5, n.cont = 3, n.cat = 1, n.levels = 3, seed = 20)
  fit <- bvsSSOM(sim$y, sim$occ.design, sim$det.design,
                 n.iter = 8000, n.burnin = 3000, n.chains = 2,
                 seed = 20, verbose = FALSE)
  ip <- inclusionProbs(fit, "psi")
  active <- tapply(sim$beta.psi, sim$occ.design$indexes_covariates,
                   function(b) any(abs(b) > 0))
  expect_equal(unname(ip > 0.5), unname(as.logical(active)))
})

## ---------------------------------------------------------------------------
## Multi-species dynamic model.
## ---------------------------------------------------------------------------

test_that("simMSDOM returns the shape bvsMSDOM expects", {
  sim <- simMSDOM(N = 4, M = 40, J = 3, n.seasons = 3,
                  n.cont = 2, n.cat = 1, seed = 30)
  expect_equal(dim(sim$y), c(4L, 40L, 3L, 3L))
  expect_equal(dim(sim$z), c(4L, 40L, 3L))
  expect_true(all(sim$y %in% c(0L, 1L)))

  ## transitions, not seasons
  expect_equal(dim(sim$per.design$X)[3L], 2L)
  expect_equal(dim(sim$col.design$X)[3L], 2L)
  expect_equal(dim(sim$det.design$X)[c(1L, 3L, 4L)], c(40L, 3L, 3L))

  expect_equal(dim(sim$beta.z), c(4L, ncol(sim$occ.design$X)))
  expect_equal(dim(sim$beta.phi), c(4L, dim(sim$per.design$X)[2L]))
  expect_equal(length(sim$mu.eta), max(sim$col.design$indexes_covariates))
})

test_that("simMSDOM honours unequal replication across site-seasons", {
  Ji <- matrix(2L, 20, 3)
  Ji[1:5, 2] <- 1L
  sim <- simMSDOM(N = 3, M = 20, J = 3, n.seasons = 3, seed = 31, Jindex = Ji)
  expect_true(all(is.na(sim$y[, , 3L, ])))
  expect_true(all(is.na(sim$y[, 1:5, 2L, 2L])))
  ## missingness is identical across species
  for (k in 2:3) expect_identical(is.na(sim$y[k, , , ]), is.na(sim$y[1L, , , ]))
})

test_that("bvsMSDOM rejects malformed input before touching nimble", {
  sim <- simMSDOM(N = 3, M = 20, J = 2, n.seasons = 3, seed = 32)

  expect_error(
    bvsMSDOM(sim$y[1, , , , drop = FALSE], sim$occ.design, sim$per.design,
             sim$col.design, sim$det.design, run = FALSE),
    "two species"
  )
  expect_error(
    bvsMSDOM(sim$y[, , , 1, drop = FALSE], sim$occ.design, sim$per.design,
             sim$col.design, sim$det.design, run = FALSE),
    "two seasons"
  )
  expect_error(
    bvsMSDOM(sim$y[, 1:10, , , drop = FALSE], sim$occ.design, sim$per.design,
             sim$col.design, sim$det.design, run = FALSE),
    "sites"
  )
  expect_error(
    bvsMSDOM(sim$y, sim$occ.design, sim$occ.design, sim$col.design,
             sim$det.design, run = FALSE),
    "3 dimensions"
  )
  expect_error(
    bvsMSDOM(sim$y, sim$occ.design, sim$per.design, sim$col.design,
             sim$det.design, tau.occ = list(intercept = c(10, 10)), run = FALSE),
    "covariate"
  )

  ## a species recorded on different visits from the others
  y2 <- sim$y
  y2[2L, 1L, 2L, 1L] <- NA
  expect_error(
    bvsMSDOM(y2, sim$occ.design, sim$per.design, sim$col.design,
             sim$det.design, run = FALSE),
    "same for every species"
  )
})

test_that("bvsMSDOM runs a short chain", {
  skip_on_cran()
  skip_without_nimble()
  sim <- simMSDOM(N = 3, M = 40, J = 2, n.seasons = 3,
                  n.cont = 2, n.cat = 1, seed = 33)
  fit <- bvsMSDOM(sim$y, sim$occ.design, sim$per.design, sim$col.design,
                  sim$det.design,
                  n.iter = 400, n.burnin = 100, n.thin = 1, n.chains = 1,
                  seed = 33, verbose = FALSE)

  expect_equal(fit$type, "MSDOM")
  ip <- inclusionProbs(fit)
  expect_named(ip, c("z", "phi", "eta", "p"))
  for (v in ip) {
    expect_true(is.matrix(v))
    expect_equal(nrow(v), 3L)
    expect_true(all(v >= 0 & v <= 1, na.rm = TRUE))
    expect_true(all(v[, 1L] == 1))
  }
  expect_output(print(fit), "multi-species dynamic")
  expect_type(medianModel(fit), "list")
})

test_that("bvsMSDOM carries a separate V, G and prior per component", {
  ## The four components may have entirely different candidate covariates, so
  ## nothing about their dimensions may be shared. Build four designs that
  ## differ in every way that matters and check they survive intact.
  M <- 30L; J <- 2L; nT <- 3L; N <- 3L

  occ <- bvsDesign(~ a1, data.frame(a1 = rnorm(M)))
  per <- bvsDesignArray(~ b1 + b2, data.frame(b1 = rnorm(M * (nT - 1L)),
                                              b2 = rnorm(M * (nT - 1L))),
                        dims = c(M, nT - 1L))
  col <- bvsDesignArray(~ c1 + g1,
                        data.frame(c1 = rnorm(M * (nT - 1L)),
                                   g1 = factor(rep(c("x", "y", "z"),
                                                   length.out = M * (nT - 1L)))),
                        dims = c(M, nT - 1L))
  det <- bvsDesignArray(~ d1 + d2 + d3,
                        data.frame(d1 = rnorm(M * J * nT), d2 = rnorm(M * J * nT),
                                   d3 = rnorm(M * J * nT)),
                        dims = c(M, J, nT))

  ## V (design columns), G (covariate groups) and ncov all differ.
  expect_equal(c(occ$numVars, per$numVars, col$numVars, det$numVars),
               c(2L, 3L, 4L, 4L))
  expect_equal(c(occ$ncov, per$ncov, col$ncov, det$ncov), c(1L, 2L, 2L, 3L))
  ## The factor in the colonisation design occupies two columns but one group.
  expect_equal(col$indexes_covariates, c(1, 2, 3, 3))

  ## Each prior is sized from its own design.
  expect_equal(length(bvsPrior(occ)$b), 2L)
  expect_equal(dim(bvsPrior(det)$B), c(4L, 4L))

  ## A prior built for one component is rejected for another.
  expect_error(
    bvsMSDOM(array(0, dim = c(N, M, J, nT)), occ, per, col, det,
             prior.occ = bvsPrior(det), run = FALSE),
    "prior.occ"
  )

  skip_on_cran()
  skip_without_nimble()
  sim <- simMSDOM(N = N, M = M, J = J, n.seasons = nT, seed = 40)
  fit <- bvsMSDOM(sim$y, occ, per, col, det,
                  n.iter = 300, n.burnin = 100, n.thin = 1, n.chains = 1,
                  seed = 40, verbose = FALSE)

  ip <- inclusionProbs(fit)
  expect_equal(vapply(ip, ncol, integer(1)),
               c(z = 2L, phi = 3L, eta = 3L, p = 4L))
  expect_equal(colnames(ip$eta), col$labels)
  expect_equal(fit$constants$V_z, 2L)
  expect_equal(fit$constants$V_p, 4L)
  expect_false(identical(fit$constants$ic_phi, fit$constants$ic_eta))
})

test_that("the Polya-Gamma update covers every design column, not just the groups", {
  ## Regression test for a bug in the four PG samplers. `index_present` was
  ## allocated and looped over with length(gamma), the number of covariate
  ## GROUPS, rather than length(indexes_covariates), the number of design
  ## COLUMNS. Whenever a categorical covariate is present the design matrix has
  ## more columns than groups, so the trailing columns were silently left out of
  ## the linear predictor used to draw Omega. Everything downstream conditions
  ## on Omega, so the whole fit was affected, and only when a factor was in the
  ## model. With all-continuous covariates columns and groups coincide and the
  ## bug was invisible.
  ##
  ## The test puts a strongly active factor in the LAST covariate group, which
  ## is exactly the group the old code dropped.
  skip_on_cran()
  skip_without_nimble()

  ## Truth: x1 null, f1 strongly active. f1 is the last covariate group and
  ## occupies the last two design columns, which is exactly what the old code
  ## dropped.
  sim <- simSSOM(M = 300, J = 5, n.cont = 1, n.cat = 1, n.levels = 3,
                 beta.psi = c(0.4, 0, 2.2, -2.2), seed = 77)

  ## Columns outnumber groups: intercept, x1, and two dummies for f1.
  expect_equal(sim$occ.design$numVars, 4L)
  expect_equal(sim$occ.design$ncov + 1L, 3L)
  expect_equal(sim$occ.design$indexes_covariates, c(1, 2, 3, 3))

  fit <- bvsSSOM(sim$y, sim$occ.design, sim$det.design,
                 n.iter = 6000, n.burnin = 2000, n.chains = 2,
                 seed = 77, verbose = FALSE)

  ip <- inclusionProbs(fit, "psi")
  expect_gt(unname(ip[3L]), 0.9)   # the factor, which the old code could not see
  expect_lt(unname(ip[2L]), 0.5)   # the null continuous covariate
})
