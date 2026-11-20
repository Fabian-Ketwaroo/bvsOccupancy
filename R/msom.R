#' Multi-species occupancy model with Bayesian variable selection
#'
#' Fits a single-season multi-species (community) occupancy model in which
#' species-level coefficients are drawn from community-level distributions, and
#' variable selection is carried out separately for each species on both the
#' occupancy and the detection component. A covariate can therefore be retained
#' for some species and dropped for others, which is usually the point of
#' fitting a community model in the first place.
#'
#' @details
#' For species \eqn{k}, occupancy and detection follow the single-season model
#' with their own inclusion vectors \eqn{\gamma^{z}_k} and \eqn{\gamma^{p}_k}.
#' The coefficients are exchangeable across species,
#' \deqn{\beta_k \sim N(\mu, \Sigma(\tau)),}
#' where \eqn{\mu} is the community mean and \eqn{\Sigma(\tau)} is block
#' diagonal: one variance per covariate, shared by the dummy columns of a
#' categorical covariate with an exchangeable correlation of one half. That
#' block structure keeps the community-level effects invariant to the choice of
#' reference level, which is not true of the usual independent-variance
#' formulation.
#'
#' The community mean and variance are updated conjugately, and only those
#' species that currently include a covariate contribute to its community
#' parameters. If no species includes a covariate at a given iteration, its
#' community parameters are drawn from the prior.
#'
#' All species are assumed to be surveyed on the same visits, so the pattern of
#' missing values must be identical across species. Sites may still be visited
#' unequal numbers of times.
#'
#' @param y Detection history: a species by sites by visits array of 0, 1 and
#'   \code{NA}.
#' @param occ.design Occupancy design, one row per site; see
#'   \code{\link{bvsDesign}}. Shared by all species.
#' @param det.design Detection design, of dimension sites by columns by visits;
#'   see \code{\link{bvsDesignArray}}. Shared by all species.
#' @param prior.occ,prior.det Priors on the community means, as returned by
#'   \code{\link{bvsPrior}}.
#' @param tau.occ,tau.det Inverse-gamma priors on the community variances, each
#'   a list with elements \code{intercept} and \code{covariate}, giving the
#'   shape and scale for the intercept variance and for the remaining
#'   covariates.
#' @inheritParams bvsSSOM
#' @param monitors Additional parameters to monitor, on top of the
#'   species-level and community-level coefficients, the community variances
#'   and the inclusion indicators.
#'
#' @return An object of class \code{"bvsOccupancy"}. For this model
#'   \code{\link{inclusionProbs}} returns a species by covariate matrix for
#'   each component.
#'
#' @references
#' Dorazio, R. M. and Royle, J. A. (2005) Estimating size and composition of
#' biological communities by modeling the occurrence of species.
#' \emph{Journal of the American Statistical Association} \strong{100}, 389--398.
#'
#' @seealso \code{\link{bvsSSOM}}, \code{\link{bvsSDOM}}, \code{\link{simMSOM}}
#'
#' @examples
#' \donttest{
#' sim <- simMSOM(N = 8, S = 120, nvisits = 4, seed = 1)
#' fit <- bvsMSOM(sim$y, sim$occ.design, sim$det.design,
#'                n.iter = 2000, n.burnin = 1000, n.chains = 2, seed = 1)
#' inclusionProbs(fit)$z
#' }
#'
#' @export
bvsMSOM <- function(y, occ.design, det.design,
                    prior.occ = NULL, prior.det = NULL,
                    tau.occ = list(intercept = c(10, 10), covariate = c(4, 0.75)),
                    tau.det = list(intercept = c(10, 10), covariate = c(4, 0.75)),
                    n.iter = 10000, n.burnin = 5000, n.thin = 5, n.chains = 2,
                    inits = NULL, monitors = NULL, monitor.z = FALSE,
                    run = TRUE, summary = TRUE, seed = NULL, verbose = TRUE) {

  call <- match.call()
  ctrl <- .mcmcControl(n.iter, n.burnin, n.thin, n.chains, run, summary, seed, verbose)

  ## ---- data and designs ---------------------------------------------------
  y <- .checkY(y, ndim = 3L)
  N <- dim(y)[1L]; S <- dim(y)[2L]; nTrep <- dim(y)[3L]
  if (N < 2L) {
    stop("a community model needs at least two species; use bvsSSOM() for one species.",
         call. = FALSE)
  }

  occ.design <- .asDesign(occ.design, "occ.design", ndim = 2L)
  det.design <- .asDesign(det.design, "det.design", ndim = 3L)
  if (nrow(occ.design$X) != S) {
    stop(sprintf("`occ.design` has %d rows but `y` has %d sites.", nrow(occ.design$X), S),
         call. = FALSE)
  }
  if (dim(det.design$X)[1L] != S || dim(det.design$X)[3L] != nTrep) {
    stop(sprintf("`det.design` is %s but `y` is %d species by %d sites by %d visits.",
                 paste(dim(det.design$X), collapse = " x "), N, S, nTrep), call. = FALSE)
  }

  ## The survey design is a property of the site, not of the species.
  miss <- is.na(y)
  ref <- miss[1L, , , drop = FALSE]
  for (k in seq_len(N)) {
    if (!identical(miss[k, , ], ref[1L, , ])) {
      stop("the pattern of missing values must be the same for every species: all species ",
           "are assumed to be recorded on the same visits. Enter a site-visit that was ",
           "surveyed but where a species went undetected as 0, not NA.", call. = FALSE)
    }
  }

  ## Left-pack using species 1, then apply the same permutation to all species.
  packed <- .leftPack(matrix(y[1L, , ], S, nTrep), det.design$X)
  if (packed$changed) {
    .tick(ctrl$verbose,
          "Some sites had gaps between surveyed visits; visits were reordered within those ",
          "sites, together with their covariates.")
    for (s in seq_len(S)) {
      o <- order(is.na(y[1L, s, ]))
      if (!identical(o, seq_len(nTrep))) y[, s, ] <- y[, s, o, drop = FALSE]
    }
    det.design$X <- packed$X
  }

  Trep <- as.integer(rowSums(!is.na(matrix(y[1L, , ], S, nTrep))))
  if (any(Trep == 0L)) {
    stop("every site must be surveyed at least once.", call. = FALSE)
  }

  prior.occ <- if (is.null(prior.occ)) bvsPrior(occ.design) else
    .checkPrior(prior.occ, occ.design, "prior.occ")
  prior.det <- if (is.null(prior.det)) bvsPrior(det.design) else
    .checkPrior(prior.det, det.design, "prior.det")
  tau.occ <- .checkTau(tau.occ, "tau.occ")
  tau.det <- .checkTau(tau.det, "tau.det")

  V_z <- occ.design$numVars; G_z <- occ.design$ncov + 1L
  V_p <- det.design$numVars; G_p <- det.design$ncov + 1L

  constants <- list(
    N = N, S = S, nTrep = nTrep, Trep = Trep,
    V_z = V_z, G_z = G_z, ic_z = occ.design$indexes_covariates,
    b_z = prior.occ$b, B_z = prior.occ$B,
    V_p = V_p, G_p = G_p, ic_p = det.design$indexes_covariates,
    b_p = prior.det$b, B_p = prior.det$B,
    az1 = tau.occ$intercept[1L], bz1 = tau.occ$intercept[2L],
    az2 = tau.occ$covariate[1L], bz2 = tau.occ$covariate[2L],
    ap1 = tau.det$intercept[1L], bp1 = tau.det$intercept[2L],
    ap2 = tau.det$covariate[1L], bp2 = tau.det$covariate[2L]
  )

  ## ---- model code ---------------------------------------------------------
  code <- nimbleCode({

    ## Community-level occupancy parameters.
    beta.comm.z[1:V_z] ~ dmnorm(b_z[1:V_z], cov = B_z[1:V_z, 1:V_z])
    TBz[1] ~ dinvgamma(az1, bz1)
    for (g in 2:G_z) {
      TBz[g] ~ dinvgamma(az2, bz2)
    }
    var.comm.z[1:V_z, 1:V_z] <- build_block_cov(TBz[1:G_z], ic_z[1:V_z], V_z)

    ## Community-level detection parameters.
    beta.comm.p[1:V_p] ~ dmnorm(b_p[1:V_p], cov = B_p[1:V_p, 1:V_p])
    TBp[1] ~ dinvgamma(ap1, bp1)
    for (g in 2:G_p) {
      TBp[g] ~ dinvgamma(ap2, bp2)
    }
    var.comm.p[1:V_p, 1:V_p] <- build_block_cov(TBp[1:G_p], ic_p[1:V_p], V_p)

    for (i in 1:N) {

      ## Species-level occupancy.
      Omega_z[i, 1:S] ~ dPG()
      gamma_z[i, 1:G_z] ~ dmnorm(b_z[1:G_z], B_z[1:G_z, 1:G_z])
      beta_z[i, 1:V_z] ~ dmnorm(beta.comm.z[1:V_z], cov = var.comm.z[1:V_z, 1:V_z])
      zlinpred[i, 1:S] <- compute_predictor(X_z[1:S, 1:V_z], beta_z[i, 1:V_z],
                                            ic_z[1:V_z], gamma_z[i, 1:G_z])[1:S, 1]

      ## Species-level detection.
      gamma_p[i, 1:G_p] ~ dmnorm(b_p[1:G_p], B_p[1:G_p, 1:G_p])
      beta_p[i, 1:V_p] ~ dmnorm(beta.comm.p[1:V_p], cov = var.comm.p[1:V_p, 1:V_p])
      for (t in 1:nTrep) {
        Omega_p[i, 1:S, t] ~ dPG()
        linpred_p[i, 1:S, t] <- compute_predictor(X_p[1:S, 1:V_p, t], beta_p[i, 1:V_p],
                                                  ic_p[1:V_p], gamma_p[i, 1:G_p])[1:S, 1]
      }

      for (s in 1:S) {
        z[i, s] ~ dbern(psi[i, s])
        logit(psi[i, s]) <- zlinpred[i, s]
        for (t in 1:Trep[s]) {
          y[i, s, t] ~ dbern(p[i, s, t] * z[i, s])
          logit(p[i, s, t]) <- linpred_p[i, s, t]
        }
      }
    }
  })

  ## ---- initial values -----------------------------------------------------
  zinit <- apply(y, c(1L, 2L), max, na.rm = TRUE)
  zinit[!is.finite(zinit)] <- 0
  base.inits <- list(
    z = matrix(as.numeric(zinit), N, S),
    beta_z = matrix(0, N, V_z), gamma_z = matrix(1, N, G_z),
    beta.comm.z = rep(0, V_z), TBz = rep(1, G_z),
    Omega_z = matrix(.rpg0(N * S), N, S),
    beta_p = matrix(0, N, V_p), gamma_p = matrix(1, N, G_p),
    beta.comm.p = rep(0, V_p), TBp = rep(1, G_p),
    Omega_p = array(.rpg0(N * S * nTrep), dim = c(N, S, nTrep))
  )
  chain.inits <- .makeInits(base.inits, ctrl$n.chains, inits)

  ## ---- model --------------------------------------------------------------
  .tick(ctrl$verbose, "Building the model.")
  model <- nimbleModel(
    code, constants = constants,
    data = list(y = y, X_z = occ.design$X, X_p = det.design$X),
    inits = chain.inits[[1L]], calculate = FALSE
  )
  if (!is.finite(model$calculate())) {
    stop("the initial log posterior is not finite. Check the design matrices, or supply ",
         "your own initial values.", call. = FALSE)
  }

  ## The community samplers address the species-level nodes by arithmetic on a
  ## flattened vector, so confirm the flattening is the one they assume.
  for (nm in c("beta_z", "gamma_z", "beta_p", "gamma_p")) {
    .assertColumnMajor(model, nm, N)
  }

  ## ---- samplers -----------------------------------------------------------
  conf <- configureMCMC(model, print = FALSE)
  .removeSamplers(conf, c("Omega_z", "gamma_z", "beta_z", "TBz", "beta.comm.z",
                          "Omega_p", "gamma_p", "beta_p", "TBp", "beta.comm.p"))

  for (i in seq_len(N)) {
    .addSiteLevelBVS(conf, "z", occ.design,
                     omega = paste0("Omega_z[", i, ", ]"),
                     gamma = paste0("gamma_z[", i, ", ]"),
                     beta  = paste0("beta_z[", i, ", ]"))

    ynodes <- model$expandNodeNames(paste0("y[", i, ", , ]"), returnScalarComponents = TRUE)
    yidx <- .parseNodeIndices(ynodes)
    Xp.flat <- .flattenDesign(det.design$X, yidx, map = c(2L, 3L))
    zindex <- .linIndex(model, paste0("z[", i, ", ]"),
                        yidx[, c(1L, 2L), drop = FALSE], "z")
    omega_pos <- .linIndex(model, paste0("Omega_p[", i, ", , ]"), yidx, "Omega_p")

    for (t in seq_len(nTrep)) {
      vi <- which(!is.na(y[i, , t]))
      if (!length(vi)) next
      conf$addSampler(
        target = paste0("Omega_p[", i, ", , ", t, "]"), type = PG_sampler_p,
        control = list(designMatrix = .slice(det.design$X, t),
                       znodes = paste0("z[", i, ", ]"), valid_indices = vi,
                       indexes_covariates = det.design$indexes_covariates,
                       fixedEffects = paste0("beta_p[", i, ", ]"),
                       Gamma = paste0("gamma_p[", i, ", ]"))
      )
    }
    conf$addSampler(
      target = paste0("gamma_p[", i, ", ]"), type = gamma_sampler_p,
      control = list(ncov = det.design$ncov,
                     indexes_covariates = det.design$indexes_covariates,
                     designMatrix = Xp.flat,
                     PG = paste0("Omega_p[", i, ", , ]"),
                     fixedEffects = paste0("beta_p[", i, ", ]"),
                     znodes = paste0("z[", i, ", ]"),
                     ynodes = paste0("y[", i, ", , ]"),
                     zindex = zindex, omega_pos = omega_pos)
    )
    conf$addSampler(
      target = paste0("beta_p[", i, ", ]"), type = beta_sampler_p,
      control = list(indexes_covariates = det.design$indexes_covariates,
                     designMatrix = Xp.flat,
                     Gamma = paste0("gamma_p[", i, ", ]"),
                     PG = paste0("Omega_p[", i, ", , ]"),
                     znodes = paste0("z[", i, ", ]"),
                     ynodes = paste0("y[", i, ", , ]"),
                     zindex = zindex, omega_pos = omega_pos)
    )
  }

  ## Community means, jointly across all columns.
  conf$addSampler(
    target = "beta.comm.z", type = Comm_multivariate_mean_sampler,
    control = list(beta.species = "beta_z", gamma.species = "gamma_z",
                   sig2.comm = "TBz",
                   indexes_covariates = occ.design$indexes_covariates)
  )
  conf$addSampler(
    target = "beta.comm.p", type = Comm_multivariate_mean_sampler,
    control = list(beta.species = "beta_p", gamma.species = "gamma_p",
                   sig2.comm = "TBp",
                   indexes_covariates = det.design$indexes_covariates)
  )

  ## Community variances, one per covariate group.
  grp.z <- split(seq_len(V_z), occ.design$indexes_covariates)
  for (g in seq_len(G_z)) {
    cols <- grp.z[[g]]
    conf$addSampler(
      target = paste0("TBz[", g, "]"), type = Comm_var_sampler,
      control = list(beta.species = paste0("beta_z[ ,", cols, "]"),
                     gamma.species = paste0("gamma_z[ ,", g, "]"),
                     beta.comm = paste0("beta.comm.z[", cols, "]"))
    )
  }
  grp.p <- split(seq_len(V_p), det.design$indexes_covariates)
  for (g in seq_len(G_p)) {
    cols <- grp.p[[g]]
    conf$addSampler(
      target = paste0("TBp[", g, "]"), type = Comm_var_sampler,
      control = list(beta.species = paste0("beta_p[ ,", cols, "]"),
                     gamma.species = paste0("gamma_p[ ,", g, "]"),
                     beta.comm = paste0("beta.comm.p[", cols, "]"))
    )
  }

  conf$resetMonitors()
  conf$addMonitors(unique(c("beta_z", "gamma_z", "beta.comm.z", "TBz",
                            "beta_p", "gamma_p", "beta.comm.p", "TBp", monitors)))
  if (isTRUE(monitor.z)) conf$addMonitors2("z")

  extra <- list(
    constants = constants,
    designs = list(z = occ.design, p = det.design),
    priors = list(z = prior.occ, p = prior.det, tau.z = tau.occ, tau.p = tau.det),
    data.dims = c(species = N, sites = S, visits = nTrep),
    gamma.map = list(
      z = list(node = "gamma_z", labels = occ.design$labels, species = TRUE, n.species = N),
      p = list(node = "gamma_p", labels = det.design$labels, species = TRUE, n.species = N)
    )
  )

  .buildAndRun(model, conf, chain.inits, ctrl, type = "MSOM", extra = extra, call = call)
}

.checkTau <- function(tau, argname) {
  if (!is.list(tau) || !all(c("intercept", "covariate") %in% names(tau))) {
    stop(sprintf("`%s` must be a list with elements `intercept` and `covariate`.", argname),
         call. = FALSE)
  }
  for (nm in c("intercept", "covariate")) {
    v <- tau[[nm]]
    if (!is.numeric(v) || length(v) != 2L || any(v <= 0)) {
      stop(sprintf("`%s$%s` must be two positive numbers: the inverse-gamma shape and scale.",
                   argname, nm), call. = FALSE)
    }
  }
  tau
}
