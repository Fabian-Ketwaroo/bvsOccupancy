#' Multi-species dynamic occupancy model with Bayesian variable selection
#'
#' Fits a community dynamic occupancy model: several species observed at the
#' same sites over several seasons, with variable selection applied separately
#' to initial occupancy, persistence, colonisation and detection, and separately
#' for every species. Species are tied together by community-level means and
#' variances, so a species with little data is shrunk towards the community.
#'
#' @details
#' For species \eqn{i} and site \eqn{s}, initial occupancy is
#' \deqn{z_{is1} \sim \mathrm{Bernoulli}(\psi_{is}), \qquad
#'       \mathrm{logit}(\psi_{is}) = X^{z}_{s,\gamma^{z}_i} \beta^{z}_{i,\gamma^{z}_i}}
#' and the state then evolves through persistence \eqn{\phi} and colonisation
#' \eqn{\eta},
#' \deqn{z_{ist} \sim \mathrm{Bernoulli}\bigl(z_{is,t-1}\phi_{is,t-1} +
#'       (1 - z_{is,t-1})\eta_{is,t-1}\bigr),}
#' with detection
#' \deqn{y_{isjt} \mid z_{ist} \sim \mathrm{Bernoulli}(z_{ist} p_{isjt}).}
#' Each of the four linear predictors carries its own vector of inclusion
#' indicators per species, and each species-level coefficient vector is drawn
#' around a community mean with a block covariance built by
#' \code{build_block_cov}, so that the fit does not depend on which level of a
#' factor happens to be the reference.
#'
#' This is the combination of \code{\link{bvsSDOM}} and \code{\link{bvsMSOM}}:
#' the dynamic structure of the former, the community structure of the latter.
#' The same samplers are used throughout; nothing in the method is new.
#'
#' @section Survey design:
#' The pattern of missing visits must be identical for every species, because a
#' visit either happened or it did not. A site-visit that was surveyed but where
#' a species went undetected is a zero, not \code{NA}. Sites may be visited
#' unequal numbers of times, and the number of visits may differ between
#' seasons.
#'
#' @param y Detection history: a species by sites by visits by seasons array of
#'   0, 1 and \code{NA}.
#' @param occ.design Initial-occupancy design, a \code{"bvsDesign"} object with
#'   one row per site; see \code{\link{bvsDesign}}.
#' @param per.design,col.design Persistence and colonisation designs, of
#'   dimension sites by columns by seasons - 1; see \code{\link{bvsDesignArray}}.
#'   A design covering all seasons is accepted and its final season dropped,
#'   since the transition out of the last season is not observed. Colonisation
#'   defaults to the persistence design.
#' @param det.design Detection design, of dimension sites by columns by visits
#'   by seasons.
#' @param prior.occ,prior.per,prior.col,prior.det Priors on the community means,
#'   as returned by \code{\link{bvsPrior}}. Defaults are built from the
#'   corresponding design.
#' @param tau.occ,tau.per,tau.col,tau.det Inverse-gamma priors on the community
#'   variances, each a list with elements \code{intercept} and \code{covariate},
#'   both of length two giving the shape and scale.
#' @param n.iter,n.burnin,n.thin,n.chains MCMC settings.
#' @param inits Optional initial values; see \code{\link{bvsSSOM}}.
#' @param monitors Additional parameters to monitor.
#' @param monitor.z Monitor the latent occupancy state, stored separately in
#'   \code{samples.z}. It is large: species by sites by seasons.
#' @param run Build and compile but do not sample.
#' @param summary Compute posterior summaries.
#' @param seed Integer seed. Chains are seeded consecutively from this value.
#' @param verbose Report progress.
#'
#' @return An object of class \code{"bvsOccupancy"} with \code{type} equal to
#'   \code{"MSDOM"}. \code{\link{inclusionProbs}} returns one species by
#'   covariate matrix per component.
#'
#' @section Scale:
#' This is the largest model in the package. The number of nodes grows as
#' species times sites times seasons, and compilation time grows with it. Start
#' from a small subset of species to confirm the designs are right before
#' running the full community.
#'
#' @seealso \code{\link{bvsSDOM}}, \code{\link{bvsMSOM}}, \code{\link{simMSDOM}}
#'
#' @examples
#' \donttest{
#' sim <- simMSDOM(N = 5, M = 80, J = 3, n.seasons = 3, seed = 1)
#' fit <- bvsMSDOM(sim$y, sim$occ.design, sim$per.design, sim$col.design,
#'                 sim$det.design,
#'                 n.iter = 2000, n.burnin = 1000, n.chains = 2, seed = 1)
#' inclusionProbs(fit, "phi")
#' }
#'
#' @export
bvsMSDOM <- function(y, occ.design, per.design, col.design = per.design, det.design,
                     prior.occ = NULL, prior.per = NULL, prior.col = NULL,
                     prior.det = NULL,
                     tau.occ = list(intercept = c(10, 10), covariate = c(4, 0.75)),
                     tau.per = list(intercept = c(10, 10), covariate = c(4, 0.75)),
                     tau.col = list(intercept = c(10, 10), covariate = c(4, 0.75)),
                     tau.det = list(intercept = c(10, 10), covariate = c(4, 0.75)),
                     n.iter = 25000, n.burnin = 10000, n.thin = 5, n.chains = 2,
                     inits = NULL, monitors = NULL, monitor.z = FALSE,
                     run = TRUE, summary = TRUE, seed = NULL, verbose = TRUE) {

  call <- match.call()
  ctrl <- .mcmcControl(n.iter, n.burnin, n.thin, n.chains, run, summary, seed, verbose)

  ## ---- data and designs ---------------------------------------------------
  y <- .checkY(y, ndim = 4L)
  N <- dim(y)[1L]; S <- dim(y)[2L]; J <- dim(y)[3L]; nT <- dim(y)[4L]
  if (N < 2L) {
    stop("a community model needs at least two species; use bvsSDOM() for one species.",
         call. = FALSE)
  }
  if (nT < 2L) {
    stop("a dynamic model needs at least two seasons; use bvsMSOM() for a single season.",
         call. = FALSE)
  }

  occ.design <- .asDesign(occ.design, "occ.design", ndim = 2L)
  per.design <- .asDesign(per.design, "per.design", ndim = 3L)
  col.design <- .asDesign(col.design, "col.design", ndim = 3L)
  det.design <- .asDesign(det.design, "det.design", ndim = 4L)

  if (nrow(occ.design$X) != S) {
    stop(sprintf("`occ.design` has %d rows but `y` has %d sites.", nrow(occ.design$X), S),
         call. = FALSE)
  }
  per.design$X <- .trimSeasons(per.design$X, S, nT, "per.design", ctrl$verbose)
  col.design$X <- .trimSeasons(col.design$X, S, nT, "col.design", ctrl$verbose)
  if (!identical(dim(det.design$X)[c(1L, 3L, 4L)], c(S, J, nT))) {
    stop(sprintf("`det.design` is %s but `y` is %d species by %d sites by %d visits by %d seasons.",
                 paste(dim(det.design$X), collapse = " x "), N, S, J, nT), call. = FALSE)
  }

  ## A visit either happened or it did not, so the missingness pattern belongs
  ## to the site, not to the species.
  miss <- is.na(y)
  ref <- miss[1L, , , , drop = FALSE]
  for (k in seq_len(N)) {
    if (!identical(miss[k, , , ], ref[1L, , , ])) {
      stop("the pattern of missing values must be the same for every species: all species ",
           "are assumed to be recorded on the same visits. Enter a site-visit that was ",
           "surveyed but where a species went undetected as 0, not NA.", call. = FALSE)
    }
  }

  ## Left-pack the surveyed visits of each site-season, using species 1 to work
  ## out the permutation and then applying it to every species.
  y1 <- array(y[1L, , , ], dim = c(S, J, nT))
  packed <- .leftPack(y1, det.design$X)
  if (packed$changed) {
    .tick(ctrl$verbose,
          "Some site-seasons had gaps between surveyed visits; visits were reordered within ",
          "those site-seasons, together with their covariates.")
    for (s in seq_len(S)) {
      for (t in seq_len(nT)) {
        o <- order(is.na(y[1L, s, , t]))
        if (!identical(o, seq_len(J))) y[, s, , t] <- y[, s, o, t, drop = FALSE]
      }
    }
    det.design$X <- packed$X
    y1 <- packed$y
  }

  Jindex <- matrix(as.integer(apply(!is.na(y1), c(1L, 3L), sum)), S, nT)
  if (any(Jindex == 0L)) {
    stop("every site must be surveyed at least once in every season. Occupancy models ",
         "tolerate unsurveyed site-seasons, but this implementation does not; ",
         "consider dropping the affected sites.", call. = FALSE)
  }

  prior.occ <- if (is.null(prior.occ)) bvsPrior(occ.design) else
    .checkPrior(prior.occ, occ.design, "prior.occ")
  prior.per <- if (is.null(prior.per)) bvsPrior(per.design) else
    .checkPrior(prior.per, per.design, "prior.per")
  prior.col <- if (is.null(prior.col)) bvsPrior(col.design) else
    .checkPrior(prior.col, col.design, "prior.col")
  prior.det <- if (is.null(prior.det)) bvsPrior(det.design) else
    .checkPrior(prior.det, det.design, "prior.det")

  tau.occ <- .checkTau(tau.occ, "tau.occ")
  tau.per <- .checkTau(tau.per, "tau.per")
  tau.col <- .checkTau(tau.col, "tau.col")
  tau.det <- .checkTau(tau.det, "tau.det")

  V_z   <- occ.design$numVars; G_z   <- occ.design$ncov + 1L
  V_phi <- per.design$numVars; G_phi <- per.design$ncov + 1L
  V_eta <- col.design$numVars; G_eta <- col.design$ncov + 1L
  V_p   <- det.design$numVars; G_p   <- det.design$ncov + 1L

  constants <- list(
    N = N, S = S, J = J, nyears = nT, Jindex = Jindex,
    V_z = V_z, G_z = G_z, ic_z = occ.design$indexes_covariates,
    b_z = prior.occ$b, B_z = prior.occ$B,
    V_phi = V_phi, G_phi = G_phi, ic_phi = per.design$indexes_covariates,
    b_phi = prior.per$b, B_phi = prior.per$B,
    V_eta = V_eta, G_eta = G_eta, ic_eta = col.design$indexes_covariates,
    b_eta = prior.col$b, B_eta = prior.col$B,
    V_p = V_p, G_p = G_p, ic_p = det.design$indexes_covariates,
    b_p = prior.det$b, B_p = prior.det$B,
    az1 = tau.occ$intercept[1L], bz1 = tau.occ$intercept[2L],
    az2 = tau.occ$covariate[1L], bz2 = tau.occ$covariate[2L],
    aphi1 = tau.per$intercept[1L], bphi1 = tau.per$intercept[2L],
    aphi2 = tau.per$covariate[1L], bphi2 = tau.per$covariate[2L],
    aeta1 = tau.col$intercept[1L], ceta1 = tau.col$intercept[2L],
    aeta2 = tau.col$covariate[1L], ceta2 = tau.col$covariate[2L],
    ap1 = tau.det$intercept[1L], bp1 = tau.det$intercept[2L],
    ap2 = tau.det$covariate[1L], bp2 = tau.det$covariate[2L]
  )

  ## ---- model code ---------------------------------------------------------
  code <- nimbleCode({

    ## Community level: one mean vector and one variance per covariate group,
    ## for each of the four components.
    beta.comm.z[1:V_z] ~ dmnorm(b_z[1:V_z], cov = B_z[1:V_z, 1:V_z])
    TBz[1] ~ dinvgamma(az1, bz1)
    for (g in 2:G_z) {
      TBz[g] ~ dinvgamma(az2, bz2)
    }
    var.comm.z[1:V_z, 1:V_z] <- build_block_cov(TBz[1:G_z], ic_z[1:V_z], V_z)

    beta.comm.phi[1:V_phi] ~ dmnorm(b_phi[1:V_phi], cov = B_phi[1:V_phi, 1:V_phi])
    TBphi[1] ~ dinvgamma(aphi1, bphi1)
    for (g in 2:G_phi) {
      TBphi[g] ~ dinvgamma(aphi2, bphi2)
    }
    var.comm.phi[1:V_phi, 1:V_phi] <- build_block_cov(TBphi[1:G_phi], ic_phi[1:V_phi], V_phi)

    beta.comm.eta[1:V_eta] ~ dmnorm(b_eta[1:V_eta], cov = B_eta[1:V_eta, 1:V_eta])
    TBeta[1] ~ dinvgamma(aeta1, ceta1)
    for (g in 2:G_eta) {
      TBeta[g] ~ dinvgamma(aeta2, ceta2)
    }
    var.comm.eta[1:V_eta, 1:V_eta] <- build_block_cov(TBeta[1:G_eta], ic_eta[1:V_eta], V_eta)

    beta.comm.p[1:V_p] ~ dmnorm(b_p[1:V_p], cov = B_p[1:V_p, 1:V_p])
    TBp[1] ~ dinvgamma(ap1, bp1)
    for (g in 2:G_p) {
      TBp[g] ~ dinvgamma(ap2, bp2)
    }
    var.comm.p[1:V_p, 1:V_p] <- build_block_cov(TBp[1:G_p], ic_p[1:V_p], V_p)

    for (i in 1:N) {

      ## Species-level coefficients and indicators. The gamma declarations are
      ## placeholders; the real model prior comes from the add / delete / swap
      ## proposal in the gamma samplers.
      Omega_z[i, 1:S] ~ dPG()
      gamma_z[i, 1:G_z] ~ dmnorm(b_z[1:G_z], B_z[1:G_z, 1:G_z])
      beta_z[i, 1:V_z] ~ dmnorm(beta.comm.z[1:V_z], cov = var.comm.z[1:V_z, 1:V_z])
      zlinpred[i, 1:S] <- compute_predictor(X_z[1:S, 1:V_z], beta_z[i, 1:V_z],
                                            ic_z[1:V_z], gamma_z[i, 1:G_z])[1:S, 1]

      gamma_phi[i, 1:G_phi] ~ dmnorm(b_phi[1:G_phi], B_phi[1:G_phi, 1:G_phi])
      beta_phi[i, 1:V_phi] ~ dmnorm(beta.comm.phi[1:V_phi],
                                    cov = var.comm.phi[1:V_phi, 1:V_phi])
      gamma_eta[i, 1:G_eta] ~ dmnorm(b_eta[1:G_eta], B_eta[1:G_eta, 1:G_eta])
      beta_eta[i, 1:V_eta] ~ dmnorm(beta.comm.eta[1:V_eta],
                                    cov = var.comm.eta[1:V_eta, 1:V_eta])
      gamma_p[i, 1:G_p] ~ dmnorm(b_p[1:G_p], B_p[1:G_p, 1:G_p])
      beta_p[i, 1:V_p] ~ dmnorm(beta.comm.p[1:V_p], cov = var.comm.p[1:V_p, 1:V_p])

      ## Transitions are defined on the nyears - 1 gaps between seasons.
      for (t in 1:(nyears - 1)) {
        Omega_phi[i, 1:S, t] ~ dPG()
        Omega_eta[i, 1:S, t] ~ dPG()
        linpred_phi[i, 1:S, t] <- compute_predictor(X_phi[1:S, 1:V_phi, t],
                                                    beta_phi[i, 1:V_phi],
                                                    ic_phi[1:V_phi],
                                                    gamma_phi[i, 1:G_phi])[1:S, 1]
        linpred_eta[i, 1:S, t] <- compute_predictor(X_eta[1:S, 1:V_eta, t],
                                                    beta_eta[i, 1:V_eta],
                                                    ic_eta[1:V_eta],
                                                    gamma_eta[i, 1:G_eta])[1:S, 1]
      }

      for (t in 1:nyears) {
        for (j in 1:J) {
          Omega_p[i, 1:S, j, t] ~ dPG()
          linpred_p[i, 1:S, j, t] <- compute_predictor(X_p[1:S, 1:V_p, j, t],
                                                       beta_p[i, 1:V_p],
                                                       ic_p[1:V_p],
                                                       gamma_p[i, 1:G_p])[1:S, 1]
        }
      }

      for (s in 1:S) {
        for (t in 1:(nyears - 1)) {
          logit(phi[i, s, t]) <- linpred_phi[i, s, t]
          logit(eta[i, s, t]) <- linpred_eta[i, s, t]
        }

        z[i, s, 1] ~ dbern(psi[i, s])
        logit(psi[i, s]) <- zlinpred[i, s]

        for (t in 2:nyears) {
          z[i, s, t] ~ dbern(z[i, s, t - 1] * phi[i, s, t - 1] +
                             (1 - z[i, s, t - 1]) * eta[i, s, t - 1])
        }

        for (t in 1:nyears) {
          for (j in 1:Jindex[s, t]) {
            y[i, s, j, t] ~ dbern(z[i, s, t] * p[i, s, j, t])
            logit(p[i, s, j, t]) <- linpred_p[i, s, j, t]
          }
        }
      }
    }
  })

  ## ---- initial values -----------------------------------------------------
  zobs <- apply(y, c(1L, 2L, 4L), max, na.rm = TRUE)
  zobs[!is.finite(zobs)] <- 0
  zinit <- array(as.numeric(zobs > 0), dim = c(N, S, nT))
  base.inits <- list(
    z = zinit,
    beta_z = matrix(0, N, V_z), gamma_z = matrix(1, N, G_z),
    beta.comm.z = rep(0, V_z), TBz = rep(1, G_z),
    Omega_z = matrix(.rpg0(N * S), N, S),
    beta_phi = matrix(0, N, V_phi), gamma_phi = matrix(1, N, G_phi),
    beta.comm.phi = rep(0, V_phi), TBphi = rep(1, G_phi),
    Omega_phi = array(.rpg0(N * S * (nT - 1L)), dim = c(N, S, nT - 1L)),
    beta_eta = matrix(0, N, V_eta), gamma_eta = matrix(1, N, G_eta),
    beta.comm.eta = rep(0, V_eta), TBeta = rep(1, G_eta),
    Omega_eta = array(.rpg0(N * S * (nT - 1L)), dim = c(N, S, nT - 1L)),
    beta_p = matrix(0, N, V_p), gamma_p = matrix(1, N, G_p),
    beta.comm.p = rep(0, V_p), TBp = rep(1, G_p),
    Omega_p = array(.rpg0(N * S * J * nT), dim = c(N, S, J, nT))
  )
  chain.inits <- .makeInits(base.inits, ctrl$n.chains, inits)

  ## ---- model --------------------------------------------------------------
  .tick(ctrl$verbose, "Building the model. This model is large; expect it to take a while.")
  model <- nimbleModel(
    code, constants = constants,
    data = list(y = y, X_z = occ.design$X, X_phi = per.design$X,
                X_eta = col.design$X, X_p = det.design$X),
    inits = chain.inits[[1L]], calculate = FALSE
  )
  if (!is.finite(model$calculate())) {
    stop("the initial log posterior is not finite. Check the design matrices, or supply ",
         "your own initial values.", call. = FALSE)
  }

  ## The community samplers address species-level nodes by arithmetic on a
  ## flattened vector, which is only valid if nimble lists them column by column.
  for (nm in c("beta_z", "gamma_z", "beta_phi", "gamma_phi",
               "beta_eta", "gamma_eta", "beta_p", "gamma_p")) {
    .assertColumnMajor(model, nm, N)
  }

  ## ---- samplers -----------------------------------------------------------
  conf <- configureMCMC(model, print = FALSE)
  .removeSamplers(conf, c("Omega_z", "gamma_z", "beta_z", "TBz", "beta.comm.z",
                          "Omega_phi", "gamma_phi", "beta_phi", "TBphi", "beta.comm.phi",
                          "Omega_eta", "gamma_eta", "beta_eta", "TBeta", "beta.comm.eta",
                          "Omega_p", "gamma_p", "beta_p", "TBp", "beta.comm.p"))

  ## The flattened designs are shared by all species, because the designs are
  ## site-level and every species is surveyed on the same visits. They are built
  ## once, from the node ordering the model reports for species 1.
  zn1 <- sprintf("z[1, 1:%d, 1:%d]", S, nT - 1L)
  yn1 <- sprintf("z[1, 1:%d, 2:%d]", S, nT)
  zix <- .parseNodeIndices(model$expandNodeNames(zn1, returnScalarComponents = TRUE))
  yix <- .parseNodeIndices(model$expandNodeNames(yn1, returnScalarComponents = TRUE))
  if (nrow(zix) != nrow(yix) || !all(yix[, 2L] == zix[, 2L]) ||
      !all(yix[, 3L] == zix[, 3L] + 1L)) {
    stop("the occupancy states of consecutive seasons could not be aligned. ",
         "This indicates an unexpected node ordering in nimble.", call. = FALSE)
  }
  ## Columns 2 and 3 of zix are site and season; column 1 is the species.
  Xphi.flat <- .flattenDesign(per.design$X, zix, map = c(2L, 3L))
  Xeta.flat <- .flattenDesign(col.design$X, zix, map = c(2L, 3L))

  yidx1 <- .parseNodeIndices(
    model$expandNodeNames("y[1, , , ]", returnScalarComponents = TRUE))
  Xp.flat <- .flattenDesign(det.design$X, yidx1, map = c(2L, 3L, 4L))

  for (i in seq_len(N)) {
    sp <- function(var, idx) paste0(var, "[", i, ", ", idx, "]")

    ## Initial occupancy: a plain site-level matrix, so the shared helper does it.
    .addSiteLevelBVS(conf, "z", occ.design,
                     omega = sp("Omega_z", ""),
                     gamma = sp("gamma_z", ""),
                     beta  = sp("beta_z", ""))

    ## Persistence and colonisation.
    zn <- sprintf("z[%d, 1:%d, 1:%d]", i, S, nT - 1L)
    yn <- sprintf("z[%d, 1:%d, 2:%d]", i, S, nT)
    for (t in seq_len(nT - 1L)) {
      conf$addSampler(
        target = paste0("Omega_phi[", i, ", , ", t, "]"), type = PG_sampler_phi,
        control = list(indexes_covariates = per.design$indexes_covariates,
                       designMatrix = .slice(per.design$X, t),
                       fixedEffects = sp("beta_phi", ""), Gamma = sp("gamma_phi", ""),
                       znodes = paste0("z[", i, ", , ", t, "]"))
      )
      conf$addSampler(
        target = paste0("Omega_eta[", i, ", , ", t, "]"), type = PG_sampler_eta,
        control = list(indexes_covariates = col.design$indexes_covariates,
                       designMatrix = .slice(col.design$X, t),
                       fixedEffects = sp("beta_eta", ""), Gamma = sp("gamma_eta", ""),
                       znodes = paste0("z[", i, ", , ", t, "]"))
      )
    }
    conf$addSampler(
      target = sp("gamma_phi", ""), type = gamma_sampler_phi,
      control = list(ncov = per.design$ncov,
                     indexes_covariates = per.design$indexes_covariates,
                     designMatrix = Xphi.flat,
                     PG = paste0("Omega_phi[", i, ", , ]"),
                     fixedEffects = sp("beta_phi", ""),
                     znodes = zn, ynodes = yn)
    )
    conf$addSampler(
      target = sp("beta_phi", ""), type = beta_sampler_phi,
      control = list(indexes_covariates = per.design$indexes_covariates,
                     designMatrix = Xphi.flat,
                     Gamma = sp("gamma_phi", ""),
                     PG = paste0("Omega_phi[", i, ", , ]"),
                     znodes = zn, ynodes = yn)
    )
    conf$addSampler(
      target = sp("gamma_eta", ""), type = gamma_sampler_eta,
      control = list(ncov = col.design$ncov,
                     indexes_covariates = col.design$indexes_covariates,
                     designMatrix = Xeta.flat,
                     PG = paste0("Omega_eta[", i, ", , ]"),
                     fixedEffects = sp("beta_eta", ""),
                     znodes = zn, ynodes = yn)
    )
    conf$addSampler(
      target = sp("beta_eta", ""), type = beta_sampler_eta,
      control = list(indexes_covariates = col.design$indexes_covariates,
                     designMatrix = Xeta.flat,
                     Gamma = sp("gamma_eta", ""),
                     PG = paste0("Omega_eta[", i, ", , ]"),
                     znodes = zn, ynodes = yn)
    )

    ## Detection. zindex and omega_pos point into this species' own slices, so
    ## they are rebuilt for every species even though the design is shared.
    yidx <- .parseNodeIndices(
      model$expandNodeNames(paste0("y[", i, ", , , ]"), returnScalarComponents = TRUE))
    zindex <- .linIndex(model, paste0("z[", i, ", , ]"),
                        yidx[, c(1L, 2L, 4L), drop = FALSE], "z")
    omega_pos <- .linIndex(model, paste0("Omega_p[", i, ", , , ]"), yidx, "Omega_p")

    for (t in seq_len(nT)) {
      for (j in seq_len(J)) {
        vi <- which(!is.na(y[i, , j, t]))
        if (!length(vi)) next
        conf$addSampler(
          target = paste0("Omega_p[", i, ", , ", j, ", ", t, "]"), type = PG_sampler_p,
          control = list(designMatrix = .slice(det.design$X, j, t),
                         znodes = paste0("z[", i, ", , ", t, "]"), valid_indices = vi,
                         indexes_covariates = det.design$indexes_covariates,
                         fixedEffects = sp("beta_p", ""), Gamma = sp("gamma_p", ""))
        )
      }
    }
    conf$addSampler(
      target = sp("gamma_p", ""), type = gamma_sampler_p,
      control = list(ncov = det.design$ncov,
                     indexes_covariates = det.design$indexes_covariates,
                     designMatrix = Xp.flat,
                     PG = paste0("Omega_p[", i, ", , , ]"),
                     fixedEffects = sp("beta_p", ""),
                     znodes = paste0("z[", i, ", , ]"),
                     ynodes = paste0("y[", i, ", , , ]"),
                     zindex = zindex, omega_pos = omega_pos)
    )
    conf$addSampler(
      target = sp("beta_p", ""), type = beta_sampler_p,
      control = list(indexes_covariates = det.design$indexes_covariates,
                     designMatrix = Xp.flat,
                     Gamma = sp("gamma_p", ""),
                     PG = paste0("Omega_p[", i, ", , , ]"),
                     znodes = paste0("z[", i, ", , ]"),
                     ynodes = paste0("y[", i, ", , , ]"),
                     zindex = zindex, omega_pos = omega_pos)
    )
  }

  ## ---- community samplers -------------------------------------------------
  comps <- list(
    list(mean = "beta.comm.z",   beta = "beta_z",   gamma = "gamma_z",
         tau = "TBz",   design = occ.design),
    list(mean = "beta.comm.phi", beta = "beta_phi", gamma = "gamma_phi",
         tau = "TBphi", design = per.design),
    list(mean = "beta.comm.eta", beta = "beta_eta", gamma = "gamma_eta",
         tau = "TBeta", design = col.design),
    list(mean = "beta.comm.p",   beta = "beta_p",   gamma = "gamma_p",
         tau = "TBp",   design = det.design)
  )
  for (cp in comps) {
    conf$addSampler(
      target = cp$mean, type = Comm_multivariate_mean_sampler,
      control = list(beta.species = cp$beta, gamma.species = cp$gamma,
                     sig2.comm = cp$tau,
                     indexes_covariates = cp$design$indexes_covariates)
    )
    grp <- split(seq_len(cp$design$numVars), cp$design$indexes_covariates)
    for (g in seq_along(grp)) {
      cols <- grp[[g]]
      conf$addSampler(
        target = paste0(cp$tau, "[", g, "]"), type = Comm_var_sampler,
        control = list(beta.species = paste0(cp$beta, "[ ,", cols, "]"),
                       gamma.species = paste0(cp$gamma, "[ ,", g, "]"),
                       beta.comm = paste0(cp$mean, "[", cols, "]"))
      )
    }
  }

  conf$resetMonitors()
  conf$addMonitors(unique(c("beta_z", "gamma_z", "beta.comm.z", "TBz",
                            "beta_phi", "gamma_phi", "beta.comm.phi", "TBphi",
                            "beta_eta", "gamma_eta", "beta.comm.eta", "TBeta",
                            "beta_p", "gamma_p", "beta.comm.p", "TBp", monitors)))
  if (isTRUE(monitor.z)) conf$addMonitors2("z")

  extra <- list(
    constants = constants,
    designs = list(z = occ.design, phi = per.design, eta = col.design, p = det.design),
    priors = list(z = prior.occ, phi = prior.per, eta = prior.col, p = prior.det,
                  tau.z = tau.occ, tau.phi = tau.per, tau.eta = tau.col, tau.p = tau.det),
    data.dims = c(species = N, sites = S, visits = J, seasons = nT),
    gamma.map = list(
      z   = list(node = "gamma_z",   labels = occ.design$labels,
                 species = TRUE, n.species = N),
      phi = list(node = "gamma_phi", labels = per.design$labels,
                 species = TRUE, n.species = N),
      eta = list(node = "gamma_eta", labels = col.design$labels,
                 species = TRUE, n.species = N),
      p   = list(node = "gamma_p",   labels = det.design$labels,
                 species = TRUE, n.species = N)
    )
  )

  .buildAndRun(model, conf, chain.inits, ctrl, type = "MSDOM", extra = extra, call = call)
}
