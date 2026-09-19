## Shared machinery behind bvsSSOM(), bvsSDOM() and bvsMSOM(). Keeping the
## build / compile / run / package steps in one place means the three model
## functions differ only in their model code and sampler assignment.

.removeSamplers <- function(conf, nodes) {
  ok <- tryCatch({
    conf$removeSamplers(nodes)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) conf$removeSampler(nodes)
  invisible(NULL)
}

.mcmcControl <- function(n.iter, n.burnin, n.thin, n.chains, run, summary, seed, verbose) {
  n.iter <- as.integer(n.iter)
  n.burnin <- as.integer(n.burnin)
  n.thin <- as.integer(n.thin)
  n.chains <- as.integer(n.chains)
  if (is.na(n.iter) || n.iter < 2L) stop("`n.iter` must be at least 2.", call. = FALSE)
  if (is.na(n.burnin) || n.burnin < 0L || n.burnin >= n.iter) {
    stop("`n.burnin` must be non-negative and smaller than `n.iter`.", call. = FALSE)
  }
  if (is.na(n.thin) || n.thin < 1L) stop("`n.thin` must be at least 1.", call. = FALSE)
  if (is.na(n.chains) || n.chains < 1L) stop("`n.chains` must be at least 1.", call. = FALSE)
  kept <- (n.iter - n.burnin) %/% n.thin
  if (kept < 2L) {
    stop("the settings keep fewer than 2 draws per chain; reduce `n.thin` or `n.burnin`.",
         call. = FALSE)
  }
  list(n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.chains = n.chains,
       run = isTRUE(run), summary = isTRUE(summary), seed = seed,
       verbose = isTRUE(verbose), n.kept = kept)
}

.buildAndRun <- function(model, conf, inits, ctrl, type, extra, call) {
  t0 <- Sys.time()
  .tick(ctrl$verbose, "Building the MCMC.")
  mcmc <- buildMCMC(conf)

  cmodel <- cmcmc <- samples <- samples2 <- smry <- NULL

  if (ctrl$run) {
    .tick(ctrl$verbose, "Compiling the model and the MCMC. This is the slow step.")
    cmodel <- compileNimble(model)
    cmcmc <- compileNimble(mcmc, project = model, resetFunctions = TRUE)

    setSeed <- if (is.null(ctrl$seed)) FALSE else
      as.numeric(ctrl$seed)[1L] + seq_len(ctrl$n.chains) - 1

    .tick(ctrl$verbose, sprintf("Running %d chain%s of %d iterations.",
                                ctrl$n.chains, if (ctrl$n.chains > 1L) "s" else "", ctrl$n.iter))
    out <- runMCMC(
      cmcmc,
      niter = ctrl$n.iter, nburnin = ctrl$n.burnin, thin = ctrl$n.thin,
      nchains = ctrl$n.chains,
      inits = if (ctrl$n.chains > 1L) inits else inits[[1L]],
      samplesAsCodaMCMC = TRUE, summary = FALSE,
      setSeed = setSeed, progressBar = ctrl$verbose
    )

    if (is.list(out) && !is.null(out$samples)) {
      samples <- out$samples
      samples2 <- out$samples2
    } else {
      samples <- out
    }
    if (!inherits(samples, "mcmc.list")) samples <- coda::mcmc.list(samples)
    if (!is.null(samples2) && !inherits(samples2, "mcmc.list")) {
      samples2 <- coda::mcmc.list(samples2)
    }
    if (ctrl$summary) smry <- .mcmcSummary(samples)
  }

  structure(
    list(
      samples = samples,
      samples.z = samples2,
      summary = smry,
      type = type,
      call = call,
      model = model,
      conf = conf,
      mcmc = mcmc,
      compiled.model = cmodel,
      compiled.mcmc = cmcmc,
      constants = extra$constants,
      designs = extra$designs,
      priors = extra$priors,
      gamma.map = extra$gamma.map,
      data.dims = extra$data.dims,
      mcmc.settings = ctrl,
      run.time = difftime(Sys.time(), t0, units = "mins")
    ),
    class = "bvsOccupancy"
  )
}

## Assign the trio of samplers that drives variable selection on one component
## whose design is a plain site-level matrix (occupancy in every model).
.addSiteLevelBVS <- function(conf, suffix, design, omega, gamma, beta) {
  conf$addSampler(
    target = omega, type = PG_sampler_psi,
    control = list(indexes_covariates = design$indexes_covariates,
                   designMatrix = design$X,
                   fixedEffects = beta, Gamma = gamma)
  )
  conf$addSampler(
    target = gamma, type = gamma_sampler_psi,
    control = list(ncov = design$ncov,
                   indexes_covariates = design$indexes_covariates,
                   designMatrix = design$X,
                   PG = omega, fixedEffects = beta)
  )
  conf$addSampler(
    target = beta, type = beta_sampler_psi,
    control = list(indexes_covariates = design$indexes_covariates,
                   designMatrix = design$X,
                   Gamma = gamma, PG = omega)
  )
  invisible(NULL)
}
