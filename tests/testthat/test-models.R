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
