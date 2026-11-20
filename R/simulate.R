## Simulators. These exist mainly so that the samplers can be exercised against
## a known truth, which is what a reproducibility check on a submitted paper
## needs. Each returns data together with ready-made design objects, so a
## simulated data set can be passed straight to the corresponding fitting
## function.

#' Simulate from a single-season occupancy model
#'
#' Generates site-level occupancy covariates and visit-level detection
#' covariates, a mixture of continuous and categorical, and simulates the
#' latent occupancy state and the detection history.
#'
#' By default some covariates are given a coefficient of exactly zero, so that
#' variable selection has something to find. A categorical covariate is either
#' wholly present or wholly absent, matching the way the sampler treats it.
#'
#' @param M Number of sites.
#' @param J Maximum number of visits per site.
#' @param n.cont Number of continuous covariates in each component.
#' @param n.cat Number of categorical covariates in each component.
#' @param n.levels Number of levels of each categorical covariate.
#' @param beta.psi,beta.p Coefficient vectors, including the intercept. Default
#'   to a sparse pattern built from the design.
#' @param Jindex Optional integer vector of length \code{M} giving the number of
#'   visits actually made to each site, for simulating unequal replication.
#' @param seed Optional random seed.
#'
#' @return A list with the detection history \code{y}, the true occupancy state
#'   \code{z}, the design objects \code{occ.design} and \code{det.design}, the
#'   true coefficients, and the covariate data frames.
#'
#' @examples
#' sim <- simSSOM(M = 100, J = 3, seed = 1)
#' str(sim$y)
#' sim$beta.psi
#'
#' @export
simSSOM <- function(M = 500, J = 5, n.cont = 3, n.cat = 1, n.levels = 3,
                    beta.psi = NULL, beta.p = NULL, Jindex = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  .checkSimDims(M, J, n.cont, n.cat, n.levels)

  site <- .simCovariates(M, n.cont, n.cat, n.levels)
  occ.design <- bvsDesign(stats::reformulate(names(site)), site)
  beta.psi <- .resolveBetas(beta.psi, occ.design, "beta.psi")
  psi <- as.vector(stats::plogis(occ.design$X %*% beta.psi))
  z <- stats::rbinom(M, 1, psi)

  visit <- .simCovariates(M * J, n.cont, n.cat, n.levels)
  det.design <- bvsDesignArray(stats::reformulate(names(visit)), visit, dims = c(M, J))
  beta.p <- .resolveBetas(beta.p, det.design, "beta.p")

  p <- matrix(0, M, J)
  for (j in seq_len(J)) {
    p[, j] <- as.vector(stats::plogis(.slice(det.design$X, j) %*% beta.p))
  }
  y <- matrix(stats::rbinom(M * J, 1, as.vector(z * p)), M, J)

  if (!is.null(Jindex)) {
    Jindex <- as.integer(Jindex)
    if (length(Jindex) != M || any(Jindex < 1L) || any(Jindex > J)) {
      stop("`Jindex` must have one entry per site, between 1 and J.", call. = FALSE)
    }
    for (i in seq_len(M)) if (Jindex[i] < J) y[i, (Jindex[i] + 1L):J] <- NA
  }

  list(y = y, z = z, psi = psi, p = p,
       occ.design = occ.design, det.design = det.design,
       beta.psi = beta.psi, beta.p = beta.p,
       site.data = site, visit.data = visit)
}

#' Simulate from a dynamic occupancy model
#'
#' Generates covariates and simulates initial occupancy, season-to-season
#' persistence and colonisation, and the detection history.
#'
#' @inheritParams simSSOM
#' @param M Number of sites.
#' @param J Maximum number of visits per site and season.
#' @param n.seasons Number of seasons.
#' @param beta.psi,beta.phi,beta.eta,beta.p Coefficient vectors, including the
#'   intercept. Default to a sparse pattern built from each design.
#' @param Jindex Optional sites by seasons integer matrix giving the number of
#'   visits actually made, for unequal replication.
#'
#' @return A list with the detection history \code{y}, the true occupancy
#'   states \code{z}, the four design objects, and the true coefficients.
#'
#' @examples
#' sim <- simSDOM(M = 60, J = 3, n.seasons = 4, seed = 1)
#' dim(sim$y)
#'
#' @export
simSDOM <- function(M = 200, J = 5, n.seasons = 10, n.cont = 3, n.cat = 1, n.levels = 3,
                    beta.psi = NULL, beta.phi = NULL, beta.eta = NULL, beta.p = NULL,
                    Jindex = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  .checkSimDims(M, J, n.cont, n.cat, n.levels)
  nT <- as.integer(n.seasons)
  if (is.na(nT) || nT < 2L) stop("`n.seasons` must be at least 2.", call. = FALSE)

  site <- .simCovariates(M, n.cont, n.cat, n.levels)
  occ.design <- bvsDesign(stats::reformulate(names(site)), site)
  beta.psi <- .resolveBetas(beta.psi, occ.design, "beta.psi")

  per.data <- .simCovariates(M * (nT - 1L), n.cont, n.cat, n.levels)
  per.design <- bvsDesignArray(stats::reformulate(names(per.data)), per.data,
                               dims = c(M, nT - 1L))
  beta.phi <- .resolveBetas(beta.phi, per.design, "beta.phi", intercept = 1)

  col.data <- .simCovariates(M * (nT - 1L), n.cont, n.cat, n.levels)
  col.design <- bvsDesignArray(stats::reformulate(names(col.data)), col.data,
                               dims = c(M, nT - 1L))
  beta.eta <- .resolveBetas(beta.eta, col.design, "beta.eta", intercept = -0.5)

  det.data <- .simCovariates(M * J * nT, n.cont, n.cat, n.levels)
  det.design <- bvsDesignArray(stats::reformulate(names(det.data)), det.data,
                               dims = c(M, J, nT))
  beta.p <- .resolveBetas(beta.p, det.design, "beta.p")

  phi <- eta <- matrix(0, M, nT - 1L)
  for (t in seq_len(nT - 1L)) {
    phi[, t] <- as.vector(stats::plogis(.slice(per.design$X, t) %*% beta.phi))
    eta[, t] <- as.vector(stats::plogis(.slice(col.design$X, t) %*% beta.eta))
  }

  z <- matrix(0L, M, nT)
  psi <- as.vector(stats::plogis(occ.design$X %*% beta.psi))
  z[, 1L] <- stats::rbinom(M, 1, psi)
  for (t in 2:nT) {
    pr <- z[, t - 1L] * phi[, t - 1L] + (1 - z[, t - 1L]) * eta[, t - 1L]
    z[, t] <- stats::rbinom(M, 1, pr)
  }

  y <- array(0, dim = c(M, J, nT))
  p <- array(0, dim = c(M, J, nT))
  for (t in seq_len(nT)) {
    for (j in seq_len(J)) {
      p[, j, t] <- as.vector(stats::plogis(.slice(det.design$X, j, t) %*% beta.p))
      y[, j, t] <- stats::rbinom(M, 1, z[, t] * p[, j, t])
    }
  }

  if (!is.null(Jindex)) {
    if (!is.matrix(Jindex) || !identical(dim(Jindex), c(as.integer(M), nT)) ||
        any(Jindex < 1L) || any(Jindex > J)) {
      stop("`Jindex` must be an M by n.seasons matrix with entries between 1 and J.",
           call. = FALSE)
    }
    for (i in seq_len(M)) {
      for (t in seq_len(nT)) {
        if (Jindex[i, t] < J) y[i, (Jindex[i, t] + 1L):J, t] <- NA
      }
    }
  }

  list(y = y, z = z, psi = psi, phi = phi, eta = eta, p = p,
       occ.design = occ.design, per.design = per.design,
       col.design = col.design, det.design = det.design,
       beta.psi = beta.psi, beta.phi = beta.phi,
       beta.eta = beta.eta, beta.p = beta.p)
}

#' Simulate from a multi-species occupancy model
#'
#' Draws species-level coefficients from community-level distributions, using
#' the same block covariance as the fitted model, then simulates occupancy and
#' detection for each species.
#'
#' Covariates whose community mean is zero are also given a very small
#' community variance, so that they are genuinely absent for every species.
#' This gives variable selection a clean target.
#'
#' @inheritParams simSSOM
#' @param N Number of species.
#' @param S Number of sites.
#' @param nvisits Maximum number of visits per site.
#' @param mu.z,mu.p Community means, one per covariate group including the
#'   intercept. Default to a sparse pattern.
#' @param tau.z,tau.p Community variances, one per covariate group. Default to
#'   a moderate variance for active covariates and a negligible one for the
#'   rest.
#' @param Trep Optional integer vector of length \code{S} giving the number of
#'   visits actually made to each site.
#'
#' @return A list with the detection history \code{y}, the true occupancy state
#'   \code{z}, the design objects, and the true species-level and
#'   community-level parameters.
#'
#' @examples
#' sim <- simMSOM(N = 5, S = 80, nvisits = 3, seed = 1)
#' dim(sim$y)
#' sim$mu.z
#'
#' @export
simMSOM <- function(N = 20, S = 250, nvisits = 5, n.cont = 3, n.cat = 1, n.levels = 3,
                    mu.z = NULL, tau.z = NULL, mu.p = NULL, tau.p = NULL,
                    Trep = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  .checkSimDims(S, nvisits, n.cont, n.cat, n.levels)
  N <- as.integer(N)
  if (is.na(N) || N < 2L) stop("`N` must be at least 2.", call. = FALSE)

  site <- .simCovariates(S, n.cont, n.cat, n.levels)
  occ.design <- bvsDesign(stats::reformulate(names(site)), site)
  visit <- .simCovariates(S * nvisits, n.cont, n.cat, n.levels)
  det.design <- bvsDesignArray(stats::reformulate(names(visit)), visit,
                               dims = c(S, nvisits))

  comm.z <- .communityTruth(occ.design, mu.z, tau.z, "z")
  comm.p <- .communityTruth(det.design, mu.p, tau.p, "p")

  beta.z <- .drawSpecies(N, comm.z$mu.col, comm.z$Sigma)
  beta.p <- .drawSpecies(N, comm.p$mu.col, comm.p$Sigma)

  psi <- matrix(0, N, S)
  z <- matrix(0L, N, S)
  y <- array(0, dim = c(N, S, nvisits))
  for (k in seq_len(N)) {
    psi[k, ] <- as.vector(stats::plogis(occ.design$X %*% beta.z[k, ]))
    z[k, ] <- stats::rbinom(S, 1, psi[k, ])
    for (t in seq_len(nvisits)) {
      pkt <- as.vector(stats::plogis(.slice(det.design$X, t) %*% beta.p[k, ]))
      y[k, , t] <- stats::rbinom(S, 1, z[k, ] * pkt)
    }
  }

  if (!is.null(Trep)) {
    Trep <- as.integer(Trep)
    if (length(Trep) != S || any(Trep < 1L) || any(Trep > nvisits)) {
      stop("`Trep` must have one entry per site, between 1 and nvisits.", call. = FALSE)
    }
    for (s in seq_len(S)) {
      if (Trep[s] < nvisits) y[, s, (Trep[s] + 1L):nvisits] <- NA
    }
  }

  list(y = y, z = z, psi = psi,
       occ.design = occ.design, det.design = det.design,
       beta.z = beta.z, beta.p = beta.p,
       mu.z = comm.z$mu, tau.z = comm.z$tau,
       mu.p = comm.p$mu, tau.p = comm.p$tau,
       site.data = site, visit.data = visit)
}

## ---------------------------------------------------------------- internal --

.checkSimDims <- function(n, J, n.cont, n.cat, n.levels) {
  if (n < 2L) stop("need at least two sites.", call. = FALSE)
  if (J < 1L) stop("need at least one visit.", call. = FALSE)
  if (n.cont < 0L || n.cat < 0L || n.cont + n.cat < 1L) {
    stop("simulate at least one covariate.", call. = FALSE)
  }
  if (n.cat > 0L && n.levels < 2L) {
    stop("`n.levels` must be at least 2.", call. = FALSE)
  }
  invisible(TRUE)
}

.simCovariates <- function(n, n.cont, n.cat, n.levels) {
  d <- list()
  for (i in seq_len(n.cont)) d[[paste0("x", i)]] <- stats::rnorm(n)
  for (i in seq_len(n.cat)) {
    d[[paste0("f", i)]] <- factor(sample.int(n.levels, n, replace = TRUE),
                                  levels = seq_len(n.levels))
  }
  as.data.frame(d, stringsAsFactors = FALSE)
}

## A sparse pattern of group-level effects: every column of a categorical
## covariate takes the same value, so the group is in or out as a whole.
.groupEffects <- function(ng, intercept = 0.5) {
  vals <- c(0.8, 0, -0.6, 0, 0.45, 0, -0.35, 0)
  c(intercept, rep(vals, length.out = max(ng - 1L, 0L)))
}

.resolveBetas <- function(beta, design, argname, intercept = 0.5) {
  ic <- design$indexes_covariates
  if (is.null(beta)) {
    return(.groupEffects(max(ic), intercept)[ic])
  }
  beta <- as.numeric(beta)
  if (length(beta) != design$numVars) {
    stop(sprintf("`%s` must have %d entries, one per column of the design matrix.",
                 argname, design$numVars), call. = FALSE)
  }
  beta
}

.blockCov <- function(tau, ic) {
  V <- length(ic)
  S <- matrix(0, V, V)
  for (a in seq_len(V)) {
    for (b in seq_len(V)) {
      if (ic[a] == ic[b]) S[a, b] <- if (a == b) tau[ic[a]] else tau[ic[a]] * 0.5
    }
  }
  S
}

.communityTruth <- function(design, mu, tau, what) {
  ic <- design$indexes_covariates
  ng <- max(ic)
  if (is.null(mu)) mu <- .groupEffects(ng, intercept = 0.5)
  mu <- as.numeric(mu)
  if (length(mu) != ng) {
    stop(sprintf("`mu.%s` must have %d entries, one per covariate group including the intercept.",
                 what, ng), call. = FALSE)
  }
  if (is.null(tau)) {
    tau <- ifelse(mu == 0, 0.01, 0.5)
    tau[1L] <- 1
  }
  tau <- as.numeric(tau)
  if (length(tau) != ng || any(tau <= 0)) {
    stop(sprintf("`tau.%s` must be %d positive variances, one per covariate group.",
                 what, ng), call. = FALSE)
  }
  list(mu = mu, tau = tau, mu.col = mu[ic], Sigma = .blockCov(tau, ic))
}

.drawSpecies <- function(N, mu, Sigma) {
  V <- length(mu)
  L <- t(chol(Sigma))
  out <- matrix(0, N, V)
  for (k in seq_len(N)) out[k, ] <- mu + as.vector(L %*% stats::rnorm(V))
  out
}
