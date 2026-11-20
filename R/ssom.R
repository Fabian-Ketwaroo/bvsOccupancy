#' Single-season occupancy model with Bayesian variable selection
#'
#' Fits the standard single-season, single-species occupancy model of
#' MacKenzie et al. (2002) with variable selection applied to both the
#' occupancy and the detection component. Each component carries its own vector
#' of binary inclusion indicators, and the covariates entering occupancy need
#' not be the same as those entering detection.
#'
#' @details
#' The model is
#' \deqn{z_i \sim \mathrm{Bernoulli}(\psi_i), \qquad
#'       \mathrm{logit}(\psi_i) = X^{\psi}_{i,\gamma^{\psi}} \beta^{\psi}_{\gamma^{\psi}}}
#' \deqn{y_{ij} \mid z_i \sim \mathrm{Bernoulli}(z_i p_{ij}), \qquad
#'       \mathrm{logit}(p_{ij}) = X^{p}_{ij,\gamma^{p}} \beta^{p}_{\gamma^{p}}}
#' where \eqn{\gamma} is a vector of inclusion indicators, one per covariate,
#' and the subscript \eqn{\gamma} denotes restriction to the active columns.
#'
#' Both Bernoulli likelihoods are written in Polya-Gamma augmented form, which
#' makes the coefficients conditionally Gaussian. They are therefore drawn in a
#' single conjugate block, and can be integrated out analytically when the
#' inclusion vector is updated, so the add / delete / swap move compares models
#' on their marginal likelihood directly.
#'
#' Sites may be surveyed unequal numbers of times. Missing visits are simply
#' left out of the likelihood; if the missing visits at a site are not the last
#' ones, the visits of that site are reordered, and the detection covariates are
#' reordered with them, before the model is built. Visits within a site are
#' exchangeable given their covariates, so this has no effect on inference.
#'
#' @param y Detection history: a sites by visits matrix of 0, 1 and \code{NA}.
#' @param occ.design Occupancy design, a \code{"bvsDesign"} object with one row
#'   per site; see \code{\link{bvsDesign}}.
#' @param det.design Detection design, a \code{"bvsDesign"} object of dimension
#'   sites by columns by visits; see \code{\link{bvsDesignArray}}.
#' @param prior.occ,prior.det Priors on the regression coefficients, as
#'   returned by \code{\link{bvsPrior}}. Defaults are constructed from the
#'   corresponding design.
#' @param n.iter,n.burnin,n.thin,n.chains MCMC settings.
#' @param inits Optional initial values. Either a named list applied to every
#'   chain, or a list of such lists with one element per chain. Anything not
#'   supplied is filled in automatically: the latent occupancy state starts at
#'   the observed maximum, coefficients at zero, and all covariates start in
#'   the model.
#' @param monitors Additional parameters to monitor, on top of
#'   \code{beta_psi}, \code{gamma_psi}, \code{beta_p} and \code{gamma_p}.
#' @param monitor.z Monitor the latent occupancy state. Stored separately in
#'   the \code{samples.z} component, since it is usually large.
#' @param run Build and compile the model but do not sample. Useful for
#'   inspecting or modifying the sampler configuration before running it.
#' @param summary Compute posterior summaries.
#' @param seed Integer seed. Chains are seeded consecutively from this value.
#' @param verbose Report progress.
#'
#' @return An object of class \code{"bvsOccupancy"}. See
#'   \code{\link{inclusionProbs}} for posterior inclusion probabilities and
#'   \code{\link{summary.bvsOccupancy}} for parameter summaries. The
#'   uncompiled \code{model}, \code{conf} and \code{mcmc} objects are returned
#'   as well, so the configuration can be inspected or extended.
#'
#' @references
#' MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle, J. A.
#' and Langtimm, C. A. (2002) Estimating site occupancy rates when detection
#' probabilities are less than one. \emph{Ecology} \strong{83}, 2248--2255.
#'
#' @seealso \code{\link{bvsSDOM}}, \code{\link{bvsMSOM}}, \code{\link{simSSOM}}
#'
#' @examples
#' \donttest{
#' sim <- simSSOM(M = 150, J = 4, seed = 1)
#' fit <- bvsSSOM(sim$y, sim$occ.design, sim$det.design,
#'                n.iter = 2000, n.burnin = 1000, n.chains = 2, seed = 1)
#' inclusionProbs(fit)
#' }
#'
#' @export
bvsSSOM <- function(y, occ.design, det.design,
                    prior.occ = NULL, prior.det = NULL,
                    n.iter = 15000, n.burnin = 5000, n.thin = 5, n.chains = 2,
                    inits = NULL, monitors = NULL, monitor.z = FALSE,
                    run = TRUE, summary = TRUE, seed = NULL, verbose = TRUE) {

  call <- match.call()
  ctrl <- .mcmcControl(n.iter, n.burnin, n.thin, n.chains, run, summary, seed, verbose)

  ## ---- data and designs ---------------------------------------------------
  y <- .checkY(y, ndim = 2L)
  occ.design <- .asDesign(occ.design, "occ.design", ndim = 2L)
  det.design <- .asDesign(det.design, "det.design", ndim = 3L)

  M <- nrow(y)
  J <- ncol(y)
  if (nrow(occ.design$X) != M) {
    stop(sprintf("`occ.design` has %d rows but `y` has %d sites.", nrow(occ.design$X), M),
         call. = FALSE)
  }
  if (dim(det.design$X)[1L] != M || dim(det.design$X)[3L] != J) {
    stop(sprintf("`det.design` is %s but `y` is %d sites by %d visits.",
                 paste(dim(det.design$X), collapse = " x "), M, J), call. = FALSE)
  }

  packed <- .leftPack(y, det.design$X)
  if (packed$changed) {
    .tick(ctrl$verbose,
          "Some sites had gaps between surveyed visits; visits were reordered within those ",
          "sites, together with their covariates, so that unsurveyed visits come last.")
  }
  y <- packed$y
  det.design$X <- packed$X

  Jindex <- as.integer(rowSums(!is.na(y)))
  if (any(Jindex == 0L)) {
    stop("every site must be surveyed at least once; drop the sites with no visits.", call. = FALSE)
  }

  prior.occ <- if (is.null(prior.occ)) bvsPrior(occ.design) else
    .checkPrior(prior.occ, occ.design, "prior.occ")
  prior.det <- if (is.null(prior.det)) bvsPrior(det.design) else
    .checkPrior(prior.det, det.design, "prior.det")

  V_psi <- occ.design$numVars; G_psi <- occ.design$ncov + 1L
  V_p <- det.design$numVars;   G_p <- det.design$ncov + 1L

  constants <- list(
    nsites = M, J = J, Jindex = Jindex,
    V_psi = V_psi, G_psi = G_psi, ic_psi = occ.design$indexes_covariates,
    b_psi = prior.occ$b, B_psi = prior.occ$B,
    V_p = V_p, G_p = G_p, ic_p = det.design$indexes_covariates,
    b_p = prior.det$b, B_p = prior.det$B
  )

  ## ---- model code ---------------------------------------------------------
  code <- nimbleCode({

    ## Occupancy: Polya-Gamma latent variables, inclusion indicators and
    ## coefficients. The indicator prior is a placeholder; the add / delete /
    ## swap sampler supplies the model prior through its proposal.
    Omega_psi[1:nsites] ~ dPG()
    gamma_psi[1:G_psi] ~ dmnorm(b_psi[1:G_psi], B_psi[1:G_psi, 1:G_psi])
    beta_psi[1:V_psi] ~ dmnorm(b_psi[1:V_psi], cov = B_psi[1:V_psi, 1:V_psi])
    linpred_psi[1:nsites] <- compute_predictor(X_psi[1:nsites, 1:V_psi],
                                               beta_psi[1:V_psi],
                                               ic_psi[1:V_psi],
                                               gamma_psi[1:G_psi])[1:nsites, 1]

    ## Detection.
    gamma_p[1:G_p] ~ dmnorm(b_p[1:G_p], B_p[1:G_p, 1:G_p])
    beta_p[1:V_p] ~ dmnorm(b_p[1:V_p], cov = B_p[1:V_p, 1:V_p])
    for (j in 1:J) {
      Omega_p[1:nsites, j] ~ dPG()
      linpred_p[1:nsites, j] <- compute_predictor(X_p[1:nsites, 1:V_p, j],
                                                  beta_p[1:V_p],
                                                  ic_p[1:V_p],
                                                  gamma_p[1:G_p])[1:nsites, 1]
    }

    ## Likelihood.
    for (s in 1:nsites) {
      z[s] ~ dbern(psi[s])
      logit(psi[s]) <- linpred_psi[s]
      for (j in 1:Jindex[s]) {
        y[s, j] ~ dbern(z[s] * p[s, j])
        logit(p[s, j]) <- linpred_p[s, j]
      }
    }
  })

  ## ---- initial values -----------------------------------------------------
  zinit <- as.numeric(apply(y, 1L, max, na.rm = TRUE))
  base.inits <- list(
    z = zinit,
    beta_psi = rep(0, V_psi), gamma_psi = rep(1, G_psi), Omega_psi = .rpg0(M),
    beta_p = rep(0, V_p), gamma_p = rep(1, G_p),
    Omega_p = matrix(.rpg0(M * J), M, J)
  )
  chain.inits <- .makeInits(base.inits, ctrl$n.chains, inits)

  ## ---- model --------------------------------------------------------------
  .tick(ctrl$verbose, "Building the model.")
  model <- nimbleModel(
    code, constants = constants,
    data = list(y = y, X_psi = occ.design$X, X_p = det.design$X),
    inits = chain.inits[[1L]], calculate = FALSE
  )
  if (!is.finite(model$calculate())) {
    stop("the initial log posterior is not finite. Check the design matrices for extreme ",
         "values, or supply your own initial values.", call. = FALSE)
  }

  ## ---- samplers -----------------------------------------------------------
  conf <- configureMCMC(model, print = FALSE)
  .removeSamplers(conf, c("Omega_psi", "gamma_psi", "beta_psi",
                          "Omega_p", "gamma_p", "beta_p"))

  .addSiteLevelBVS(conf, "psi", occ.design,
                   omega = "Omega_psi", gamma = "gamma_psi", beta = "beta_psi")

  ## Detection works on the flattened vector of observed site-visits. Both the
  ## row order of the flattened design and the two index vectors are taken from
  ## the node names the model reports, so they cannot drift out of step with
  ## `values(model, ynodes)`.
  ynodes <- model$expandNodeNames("y", returnScalarComponents = TRUE)
  yidx <- .parseNodeIndices(ynodes)
  Xp.flat <- .flattenDesign(det.design$X, yidx, map = c(1L, 2L))
  zindex <- .linIndex(model, "z", yidx[, 1L, drop = FALSE], "z")
  omega_pos <- .linIndex(model, "Omega_p", yidx[, c(1L, 2L), drop = FALSE], "Omega_p")

  for (j in seq_len(J)) {
    vi <- which(!is.na(y[, j]))
    if (!length(vi)) next
    conf$addSampler(
      target = paste0("Omega_p[, ", j, "]"), type = PG_sampler_p,
      control = list(designMatrix = .slice(det.design$X, j),
                     znodes = "z", valid_indices = vi,
                     indexes_covariates = det.design$indexes_covariates,
                     fixedEffects = "beta_p", Gamma = "gamma_p")
    )
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
  conf$addMonitors(unique(c("beta_psi", "gamma_psi", "beta_p", "gamma_p", monitors)))
  if (isTRUE(monitor.z)) conf$addMonitors2("z")

  extra <- list(
    constants = constants,
    designs = list(psi = occ.design, p = det.design),
    priors = list(psi = prior.occ, p = prior.det),
    data.dims = c(sites = M, visits = J),
    gamma.map = list(
      psi = list(node = "gamma_psi", labels = occ.design$labels, species = FALSE),
      p   = list(node = "gamma_p",   labels = det.design$labels, species = FALSE)
    )
  )

  .buildAndRun(model, conf, chain.inits, ctrl, type = "SSOM", extra = extra, call = call)
}
