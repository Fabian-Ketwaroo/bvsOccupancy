#' Dynamic occupancy model with Bayesian variable selection
#'
#' Fits the multi-season dynamic occupancy model of MacKenzie et al. (2003)
#' with variable selection applied to all four components: initial occupancy,
#' persistence, colonisation and detection. Each component has its own design
#' matrix and its own vector of inclusion indicators, so a covariate can drive
#' colonisation without being retained for persistence.
#'
#' @details
#' Occupancy in the first season follows the single-season model. Thereafter
#' \deqn{z_{i,t} \mid z_{i,t-1} \sim \mathrm{Bernoulli}\!\left(z_{i,t-1}\phi_{i,t-1}
#'       + (1 - z_{i,t-1})\eta_{i,t-1}\right)}
#' with \eqn{\phi} the probability that an occupied site stays occupied and
#' \eqn{\eta} the probability that an empty site is colonised. Both are modelled
#' on the logit scale with their own selected covariates, as is detection.
#'
#' Because a site contributes to the persistence likelihood only when it was
#' occupied in the preceding season, and to the colonisation likelihood only
#' when it was not, the relevant samplers subset their design matrix by the
#' current value of the latent state at every iteration.
#'
#' Sites may be surveyed unequal numbers of times in different seasons. As in
#' \code{\link{bvsSSOM}}, visits are reordered within a site and season if
#' necessary so that unsurveyed visits come last.
#'
#' @param y Detection history: a sites by visits by seasons array of 0, 1 and
#'   \code{NA}.
#' @param occ.design Design for initial occupancy, one row per site; see
#'   \code{\link{bvsDesign}}.
#' @param per.design Design for persistence, of dimension sites by columns by
#'   seasons. Only the first \code{nseasons - 1} slices are used, since the
#'   transition out of the final season is not observed. A design with either
#'   \code{nseasons} or \code{nseasons - 1} slices is accepted.
#' @param col.design Design for colonisation, in the same layout as
#'   \code{per.design}. Defaults to \code{per.design}.
#' @param det.design Design for detection, of dimension sites by columns by
#'   visits by seasons; see \code{\link{bvsDesignArray}}.
#' @param prior.occ,prior.per,prior.col,prior.det Priors on the regression
#'   coefficients, as returned by \code{\link{bvsPrior}}. Defaults are
#'   constructed from the corresponding design.
#' @inheritParams bvsSSOM
#' @param monitors Additional parameters to monitor, on top of the coefficients
#'   and indicators of all four components.
#'
#' @return An object of class \code{"bvsOccupancy"}.
#'
#' @references
#' MacKenzie, D. I., Nichols, J. D., Hines, J. E., Knutson, M. G. and Franklin,
#' A. B. (2003) Estimating site occupancy, colonization, and local extinction
#' when a species is detected imperfectly. \emph{Ecology} \strong{84}, 2200--2207.
#'
#' @seealso \code{\link{bvsSSOM}}, \code{\link{bvsMSOM}}, \code{\link{simSDOM}}
#'
#' @examples
#' \donttest{
#' sim <- simSDOM(M = 100, J = 4, n.seasons = 5, seed = 1)
#' fit <- bvsSDOM(sim$y, sim$occ.design, sim$per.design, sim$col.design,
#'                sim$det.design,
#'                n.iter = 2000, n.burnin = 1000, n.chains = 2, seed = 1)
#' inclusionProbs(fit)
#' }
#'
#' @export
bvsSDOM <- function(y, occ.design, per.design, col.design = per.design, det.design,
                    prior.occ = NULL, prior.per = NULL, prior.col = NULL, prior.det = NULL,
                    n.iter = 25000, n.burnin = 10000, n.thin = 5, n.chains = 2,
                    inits = NULL, monitors = NULL, monitor.z = FALSE,
                    run = TRUE, summary = TRUE, seed = NULL, verbose = TRUE) {

  call <- match.call()
  ctrl <- .mcmcControl(n.iter, n.burnin, n.thin, n.chains, run, summary, seed, verbose)

  ## ---- data and designs ---------------------------------------------------
  y <- .checkY(y, ndim = 3L)
  M <- dim(y)[1L]; J <- dim(y)[2L]; nT <- dim(y)[3L]
  if (nT < 2L) {
    stop("a dynamic model needs at least two seasons; use bvsSSOM() for a single season.",
         call. = FALSE)
  }

  occ.design <- .asDesign(occ.design, "occ.design", ndim = 2L)
  per.design <- .asDesign(per.design, "per.design", ndim = 3L)
  col.design <- .asDesign(col.design, "col.design", ndim = 3L)
  det.design <- .asDesign(det.design, "det.design", ndim = 4L)

  if (nrow(occ.design$X) != M) {
    stop(sprintf("`occ.design` has %d rows but `y` has %d sites.", nrow(occ.design$X), M),
         call. = FALSE)
  }
  per.design$X <- .trimSeasons(per.design$X, M, nT, "per.design", ctrl$verbose)
  col.design$X <- .trimSeasons(col.design$X, M, nT, "col.design", ctrl$verbose)
  if (!identical(dim(det.design$X)[c(1L, 3L, 4L)], c(M, J, nT))) {
    stop(sprintf("`det.design` is %s but `y` is %d sites by %d visits by %d seasons.",
                 paste(dim(det.design$X), collapse = " x "), M, J, nT), call. = FALSE)
  }

  packed <- .leftPack(y, det.design$X)
  if (packed$changed) {
    .tick(ctrl$verbose,
          "Some site-seasons had gaps between surveyed visits; visits were reordered within ",
          "those site-seasons, together with their covariates.")
  }
  y <- packed$y
  det.design$X <- packed$X

  Jindex <- matrix(as.integer(apply(!is.na(y), c(1L, 3L), sum)), M, nT)
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

  V_psi <- occ.design$numVars; G_psi <- occ.design$ncov + 1L
  V_phi <- per.design$numVars; G_phi <- per.design$ncov + 1L
  V_eta <- col.design$numVars; G_eta <- col.design$ncov + 1L
  V_p   <- det.design$numVars; G_p   <- det.design$ncov + 1L

  constants <- list(
    M = M, J = J, nyears = nT, Jindex = Jindex,
    V_psi = V_psi, G_psi = G_psi, ic_psi = occ.design$indexes_covariates,
    b_psi = prior.occ$b, B_psi = prior.occ$B,
    V_phi = V_phi, G_phi = G_phi, ic_phi = per.design$indexes_covariates,
    b_phi = prior.per$b, B_phi = prior.per$B,
    V_eta = V_eta, G_eta = G_eta, ic_eta = col.design$indexes_covariates,
    b_eta = prior.col$b, B_eta = prior.col$B,
    V_p = V_p, G_p = G_p, ic_p = det.design$indexes_covariates,
    b_p = prior.det$b, B_p = prior.det$B
  )

  ## ---- model code ---------------------------------------------------------
  code <- nimbleCode({

    ## Initial occupancy.
    Omega_psi[1:M] ~ dPG()
    gamma_psi[1:G_psi] ~ dmnorm(b_psi[1:G_psi], B_psi[1:G_psi, 1:G_psi])
    beta_psi[1:V_psi] ~ dmnorm(b_psi[1:V_psi], cov = B_psi[1:V_psi, 1:V_psi])
    linpred_psi[1:M] <- compute_predictor(X_psi[1:M, 1:V_psi], beta_psi[1:V_psi],
                                          ic_psi[1:V_psi], gamma_psi[1:G_psi])[1:M, 1]

    ## Persistence and colonisation.
    gamma_phi[1:G_phi] ~ dmnorm(b_phi[1:G_phi], B_phi[1:G_phi, 1:G_phi])
    beta_phi[1:V_phi] ~ dmnorm(b_phi[1:V_phi], cov = B_phi[1:V_phi, 1:V_phi])
    gamma_eta[1:G_eta] ~ dmnorm(b_eta[1:G_eta], B_eta[1:G_eta, 1:G_eta])
    beta_eta[1:V_eta] ~ dmnorm(b_eta[1:V_eta], cov = B_eta[1:V_eta, 1:V_eta])

    for (t in 1:(nyears - 1)) {
      Omega_phi[1:M, t] ~ dPG()
      Omega_eta[1:M, t] ~ dPG()
      linpred_phi[1:M, t] <- compute_predictor(X_phi[1:M, 1:V_phi, t], beta_phi[1:V_phi],
                                               ic_phi[1:V_phi], gamma_phi[1:G_phi])[1:M, 1]
      linpred_eta[1:M, t] <- compute_predictor(X_eta[1:M, 1:V_eta, t], beta_eta[1:V_eta],
                                               ic_eta[1:V_eta], gamma_eta[1:G_eta])[1:M, 1]
    }

    ## Detection.
    gamma_p[1:G_p] ~ dmnorm(b_p[1:G_p], B_p[1:G_p, 1:G_p])
    beta_p[1:V_p] ~ dmnorm(b_p[1:V_p], cov = B_p[1:V_p, 1:V_p])
    for (t in 1:nyears) {
      for (j in 1:J) {
        Omega_p[1:M, j, t] ~ dPG()
        linpred_p[1:M, j, t] <- compute_predictor(X_p[1:M, 1:V_p, j, t], beta_p[1:V_p],
                                                  ic_p[1:V_p], gamma_p[1:G_p])[1:M, 1]
      }
    }

    for (i in 1:M) {
      for (t in 1:(nyears - 1)) {
        logit(phi[i, t]) <- linpred_phi[i, t]
        logit(eta[i, t]) <- linpred_eta[i, t]
      }
    }

    for (i in 1:M) {
      z[i, 1] ~ dbern(psi[i])
      logit(psi[i]) <- linpred_psi[i]
    }

    for (i in 1:M) {
      for (t in 2:nyears) {
        z[i, t] ~ dbern(z[i, t - 1] * phi[i, t - 1] + (1 - z[i, t - 1]) * eta[i, t - 1])
      }
    }

    for (i in 1:M) {
      for (t in 1:nyears) {
        for (j in 1:Jindex[i, t]) {
          y[i, j, t] ~ dbern(z[i, t] * p[i, j, t])
          logit(p[i, j, t]) <- linpred_p[i, j, t]
        }
      }
    }
  })

  ## ---- initial values -----------------------------------------------------
  zobs <- apply(y, c(1L, 3L), max, na.rm = TRUE)
  zobs[!is.finite(zobs)] <- 0
  ## A site seen at any point must be treated as occupied throughout, otherwise
  ## the transition likelihood can start at an impossible state.
  zinit <- matrix(as.numeric(zobs > 0), M, nT)
  base.inits <- list(
    z = zinit,
    beta_psi = rep(0, V_psi), gamma_psi = rep(1, G_psi), Omega_psi = .rpg0(M),
    beta_phi = rep(0, V_phi), gamma_phi = rep(1, G_phi),
    Omega_phi = matrix(.rpg0(M * (nT - 1L)), M, nT - 1L),
    beta_eta = rep(0, V_eta), gamma_eta = rep(1, G_eta),
    Omega_eta = matrix(.rpg0(M * (nT - 1L)), M, nT - 1L),
    beta_p = rep(0, V_p), gamma_p = rep(1, G_p),
    Omega_p = array(.rpg0(M * J * nT), dim = c(M, J, nT))
  )
  chain.inits <- .makeInits(base.inits, ctrl$n.chains, inits)

  ## ---- model --------------------------------------------------------------
  .tick(ctrl$verbose, "Building the model.")
  model <- nimbleModel(
    code, constants = constants,
    data = list(y = y, X_psi = occ.design$X, X_phi = per.design$X,
                X_eta = col.design$X, X_p = det.design$X),
    inits = chain.inits[[1L]], calculate = FALSE
  )
  if (!is.finite(model$calculate())) {
    stop("the initial log posterior is not finite. Check the design matrices, or supply ",
         "your own initial values.", call. = FALSE)
  }

  ## ---- samplers -----------------------------------------------------------
  conf <- configureMCMC(model, print = FALSE)
  .removeSamplers(conf, c("Omega_psi", "gamma_psi", "beta_psi",
                          "Omega_phi", "gamma_phi", "beta_phi",
                          "Omega_eta", "gamma_eta", "beta_eta",
                          "Omega_p", "gamma_p", "beta_p"))

  .addSiteLevelBVS(conf, "psi", occ.design,
                   omega = "Omega_psi", gamma = "gamma_psi", beta = "beta_psi")

  ## Persistence and colonisation compare the state in season t with the state
  ## in season t + 1, so the two node vectors must be in the same order and the
  ## flattened design must follow them row for row.
  zn.str <- sprintf("z[1:%d, 1:%d]", M, nT - 1L)
  yn.str <- sprintf("z[1:%d, 2:%d]", M, nT)
  zn <- model$expandNodeNames(zn.str, returnScalarComponents = TRUE)
  yn <- model$expandNodeNames(yn.str, returnScalarComponents = TRUE)
  zix <- .parseNodeIndices(zn)
  yix <- .parseNodeIndices(yn)
  if (nrow(zix) != nrow(yix) || !all(yix[, 1L] == zix[, 1L]) ||
      !all(yix[, 2L] == zix[, 2L] + 1L)) {
    stop("the occupancy states of consecutive seasons could not be aligned. ",
         "This indicates an unexpected node ordering in nimble.", call. = FALSE)
  }

  for (nm in c("Omega_phi", "Omega_eta")) {
    pix <- .parseNodeIndices(model$expandNodeNames(nm, returnScalarComponents = TRUE))
    if (nrow(pix) != nrow(zix) || !all(pix == zix)) {
      stop(sprintf("`%s` is not ordered like the occupancy states it belongs to.", nm),
           call. = FALSE)
    }
  }

  Xphi.flat <- .flattenDesign(per.design$X, zix, map = c(1L, 2L))
  Xeta.flat <- .flattenDesign(col.design$X, zix, map = c(1L, 2L))

  for (t in seq_len(nT - 1L)) {
    conf$addSampler(
      target = paste0("Omega_phi[, ", t, "]"), type = PG_sampler_phi,
      control = list(indexes_covariates = per.design$indexes_covariates,
                     designMatrix = .slice(per.design$X, t),
                     fixedEffects = "beta_phi", Gamma = "gamma_phi",
                     znodes = paste0("z[, ", t, "]"))
    )
    conf$addSampler(
      target = paste0("Omega_eta[, ", t, "]"), type = PG_sampler_eta,
      control = list(indexes_covariates = col.design$indexes_covariates,
                     designMatrix = .slice(col.design$X, t),
                     fixedEffects = "beta_eta", Gamma = "gamma_eta",
                     znodes = paste0("z[, ", t, "]"))
    )
  }

  conf$addSampler(
    target = "gamma_phi", type = gamma_sampler_phi,
    control = list(ncov = per.design$ncov,
                   indexes_covariates = per.design$indexes_covariates,
                   designMatrix = Xphi.flat, PG = "Omega_phi", fixedEffects = "beta_phi",
                   znodes = zn.str, ynodes = yn.str)
  )
  conf$addSampler(
    target = "beta_phi", type = beta_sampler_phi,
    control = list(indexes_covariates = per.design$indexes_covariates,
                   designMatrix = Xphi.flat, Gamma = "gamma_phi", PG = "Omega_phi",
                   znodes = zn.str, ynodes = yn.str)
  )
  conf$addSampler(
    target = "gamma_eta", type = gamma_sampler_eta,
    control = list(ncov = col.design$ncov,
                   indexes_covariates = col.design$indexes_covariates,
                   designMatrix = Xeta.flat, PG = "Omega_eta", fixedEffects = "beta_eta",
                   znodes = zn.str, ynodes = yn.str)
  )
  conf$addSampler(
    target = "beta_eta", type = beta_sampler_eta,
    control = list(indexes_covariates = col.design$indexes_covariates,
                   designMatrix = Xeta.flat, Gamma = "gamma_eta", PG = "Omega_eta",
                   znodes = zn.str, ynodes = yn.str)
  )

  ## Detection.
  ynodes <- model$expandNodeNames("y", returnScalarComponents = TRUE)
  yidx <- .parseNodeIndices(ynodes)
  Xp.flat <- .flattenDesign(det.design$X, yidx, map = c(1L, 2L, 3L))
  zindex <- .linIndex(model, "z", yidx[, c(1L, 3L), drop = FALSE], "z")
  omega_pos <- .linIndex(model, "Omega_p", yidx, "Omega_p")

  for (t in seq_len(nT)) {
    for (j in seq_len(J)) {
      vi <- which(!is.na(y[, j, t]))
      if (!length(vi)) next
      conf$addSampler(
        target = paste0("Omega_p[, ", j, ", ", t, "]"), type = PG_sampler_p,
        control = list(designMatrix = .slice(det.design$X, j, t),
                       znodes = paste0("z[, ", t, "]"), valid_indices = vi,
                       indexes_covariates = det.design$indexes_covariates,
                       fixedEffects = "beta_p", Gamma = "gamma_p")
      )
    }
  }
  conf$addSampler(
    target = "gamma_p", type = gamma_sampler_p,
    control = list(ncov = det.design$ncov,
                   indexes_covariates = det.design$indexes_covariates,
                   designMatrix = Xp.flat, PG = "Omega_p", fixedEffects = "beta_p",
                   znodes = "z", ynodes = "y", zindex = zindex, omega_pos = omega_pos)
  )
  conf$addSampler(
    target = "beta_p", type = beta_sampler_p,
    control = list(indexes_covariates = det.design$indexes_covariates,
                   designMatrix = Xp.flat, Gamma = "gamma_p", PG = "Omega_p",
                   znodes = "z", ynodes = "y", zindex = zindex, omega_pos = omega_pos)
  )

  conf$resetMonitors()
  conf$addMonitors(unique(c("beta_psi", "gamma_psi", "beta_phi", "gamma_phi",
                            "beta_eta", "gamma_eta", "beta_p", "gamma_p", monitors)))
  if (isTRUE(monitor.z)) conf$addMonitors2("z")

  extra <- list(
    constants = constants,
    designs = list(psi = occ.design, phi = per.design, eta = col.design, p = det.design),
    priors = list(psi = prior.occ, phi = prior.per, eta = prior.col, p = prior.det),
    data.dims = c(sites = M, visits = J, seasons = nT),
    gamma.map = list(
      psi = list(node = "gamma_psi", labels = occ.design$labels, species = FALSE),
      phi = list(node = "gamma_phi", labels = per.design$labels, species = FALSE),
      eta = list(node = "gamma_eta", labels = col.design$labels, species = FALSE),
      p   = list(node = "gamma_p",   labels = det.design$labels, species = FALSE)
    )
  )

  .buildAndRun(model, conf, chain.inits, ctrl, type = "SDOM", extra = extra, call = call)
}

## Persistence and colonisation only need the first nseasons - 1 slices.
.trimSeasons <- function(X, M, nT, argname, verbose) {
  d <- dim(X)
  if (d[1L] != M) {
    stop(sprintf("`%s` has %d sites but `y` has %d.", argname, d[1L], M), call. = FALSE)
  }
  if (d[3L] == nT) {
    .tick(verbose, sprintf("`%s` has %d seasons; using the first %d, since the transition ",
                           argname, nT, nT - 1L),
          "out of the final season is not observed.")
    X <- X[, , seq_len(nT - 1L), drop = FALSE]
  } else if (d[3L] != nT - 1L) {
    stop(sprintf("`%s` must have %d or %d seasons but has %d.",
                 argname, nT - 1L, nT, d[3L]), call. = FALSE)
  }
  X
}
