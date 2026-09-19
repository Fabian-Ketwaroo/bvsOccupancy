## Methods for fitted models.

#' Posterior inclusion probabilities
#'
#' Extracts the posterior probability that each covariate is in the model, that
#' is, the posterior mean of its inclusion indicator. The intercept is always
#' retained, so its probability is one by construction.
#'
#' These are the quantities to report as the output of variable selection. A
#' common convention treats covariates with probability above 0.5 as belonging
#' to the median probability model, though the whole distribution is more
#' informative than any single threshold.
#'
#' @param object A fitted model from \code{\link{bvsSSOM}},
#'   \code{\link{bvsSDOM}} or \code{\link{bvsMSOM}}.
#' @param component Which component to return. The default returns all of them.
#'   Single-season models have \code{"psi"} and \code{"p"}; dynamic models add
#'   \code{"phi"} and \code{"eta"}; multi-species models have \code{"z"} and
#'   \code{"p"}.
#' @param ... Ignored.
#'
#' @return For single-species models, a named list of numeric vectors, one per
#'   component, with the covariates as names. For multi-species models, a list
#'   of species by covariate matrices. If \code{component} names a single
#'   component, that element is returned directly.
#'
#' @examples
#' \donttest{
#' sim <- simSSOM(M = 150, J = 4, seed = 1)
#' fit <- bvsSSOM(sim$y, sim$occ.design, sim$det.design,
#'                n.iter = 2000, n.burnin = 1000, seed = 1)
#' inclusionProbs(fit, "psi")
#' }
#'
#' @export
inclusionProbs <- function(object, component = NULL, ...) {
  if (!inherits(object, "bvsOccupancy")) {
    stop("`object` must be a fitted bvsOccupancy model.", call. = FALSE)
  }
  if (is.null(object$samples)) {
    stop("this model has not been run, so there are no posterior samples. ",
         "Refit with run = TRUE.", call. = FALSE)
  }
  map <- object$gamma.map
  if (!is.null(component)) {
    unknown <- setdiff(component, names(map))
    if (length(unknown)) {
      stop(sprintf("unknown component%s %s. Available: %s.",
                   if (length(unknown) > 1L) "s" else "",
                   paste(sQuote(unknown), collapse = ", "),
                   paste(sQuote(names(map)), collapse = ", ")), call. = FALSE)
    }
    map <- map[component]
  }

  draws <- as.matrix(object$samples)
  out <- lapply(map, function(m) {
    if (isTRUE(m$species)) {
      nsp <- m$n.species
      ng <- length(m$labels)
      res <- matrix(NA_real_, nsp, ng,
                    dimnames = list(paste0("species", seq_len(nsp)), m$labels))
      for (k in seq_len(nsp)) {
        for (g in seq_len(ng)) {
          nm <- sprintf("%s[%d, %d]", m$node, k, g)
          if (nm %in% colnames(draws)) res[k, g] <- mean(draws[, nm])
        }
      }
      res
    } else {
      ng <- length(m$labels)
      res <- setNames(rep(NA_real_, ng), m$labels)
      for (g in seq_len(ng)) {
        nm <- sprintf("%s[%d]", m$node, g)
        if (nm %in% colnames(draws)) res[g] <- mean(draws[, nm])
      }
      res
    }
  })
  if (length(out) == 1L) out[[1L]] else out
}

#' @export
print.bvsOccupancy <- function(x, ...) {
  lbl <- switch(x$type,
                SSOM = "single-season occupancy model",
                SDOM = "dynamic occupancy model",
                MSOM = "multi-species occupancy model",
                x$type)
  cat("<bvsOccupancy>", lbl, "with Bayesian variable selection\n\n")
  cat("Call:\n  ", paste(deparse(x$call), collapse = "\n   "), "\n\n", sep = "")
  cat("Data: ", paste(sprintf("%d %s", x$data.dims, names(x$data.dims)), collapse = ", "),
      "\n", sep = "")
  cat("Candidate covariates:\n")
  for (nm in names(x$designs)) {
    d <- x$designs[[nm]]
    cat(sprintf("  %-4s %2d covariate%s, %d columns\n", nm, d$ncov,
                if (d$ncov == 1L) "" else "s", d$numVars))
  }
  s <- x$mcmc.settings
  cat(sprintf("\nMCMC: %d iterations, %d burn-in, thin %d, %d chain%s (%d draws kept per chain)\n",
              s$n.iter, s$n.burnin, s$n.thin, s$n.chains,
              if (s$n.chains > 1L) "s" else "", s$n.kept))
  if (is.null(x$samples)) {
    cat("\nNot yet run. The model, configuration and MCMC objects are available for editing.\n")
  } else {
    cat(sprintf("Run time: %.1f minutes\n", as.numeric(x$run.time)))
    cat("\nPosterior inclusion probabilities:\n")
    ip <- inclusionProbs(x)
    if (!is.list(ip)) ip <- list(ip)
    for (nm in names(ip)) {
      v <- ip[[nm]]
      if (is.matrix(v)) {
        cat(sprintf("  %s: species-specific, use inclusionProbs(fit, \"%s\")\n", nm, nm))
        cat(sprintf("     mean across species: %s\n",
                    paste(sprintf("%s=%.2f", colnames(v), colMeans(v, na.rm = TRUE)),
                          collapse = "  ")))
      } else {
        cat(sprintf("  %-4s %s\n", nm,
                    paste(sprintf("%s=%.2f", names(v), v), collapse = "  ")))
      }
    }
    cat("\nUse summary() for parameter estimates.\n")
  }
  invisible(x)
}

#' Summarise a fitted model
#'
#' @param object A fitted model from \code{\link{bvsSSOM}},
#'   \code{\link{bvsSDOM}} or \code{\link{bvsMSOM}}.
#' @param pars Optional regular expression; only parameters whose names match
#'   are returned. For example \code{"^beta_psi"}.
#' @param digits Number of digits to round to.
#' @param ... Ignored.
#'
#' @return A data frame with the posterior mean, standard deviation, 2.5th,
#'   50th and 97.5th percentiles, the potential scale reduction factor and the
#'   effective sample size.
#'
#' @section Reading the output:
#' Coefficients are sampled under variable selection, so the posterior of a
#' coefficient is a mixture over the models that contain it and those that do
#' not. When a covariate is excluded its coefficient is left at its last value
#' and does not enter the likelihood, so a coefficient summary is only
#' interpretable together with the inclusion probability from
#' \code{\link{inclusionProbs}}. A scale reduction factor of \code{NA} is
#' normal for an indicator that never moves, such as the intercept.
#'
#' @export
summary.bvsOccupancy <- function(object, pars = NULL, digits = 3, ...) {
  if (is.null(object$samples)) {
    stop("this model has not been run, so there are no posterior samples.", call. = FALSE)
  }
  s <- object$summary
  if (is.null(s)) s <- .mcmcSummary(object$samples)
  if (!is.null(pars)) {
    keep <- grepl(pars, rownames(s))
    if (!any(keep)) {
      stop(sprintf("no monitored parameter matches %s.", sQuote(pars)), call. = FALSE)
    }
    s <- s[keep, , drop = FALSE]
  }
  num <- vapply(s, is.numeric, logical(1))
  s[num] <- lapply(s[num], round, digits = digits)
  s
}

#' Posterior samples as a coda object
#'
#' @param x A fitted model from \code{\link{bvsSSOM}}, \code{\link{bvsSDOM}} or
#'   \code{\link{bvsMSOM}}.
#' @param ... Ignored.
#'
#' @return An \code{\link[coda]{mcmc.list}}, ready for \pkg{coda} or
#'   \pkg{MCMCvis}.
#'
#' @importFrom coda as.mcmc.list
#' @method as.mcmc.list bvsOccupancy
#' @export
as.mcmc.list.bvsOccupancy <- function(x, ...) {
  if (is.null(x$samples)) {
    stop("this model has not been run, so there are no posterior samples.", call. = FALSE)
  }
  x$samples
}

#' Median probability model
#'
#' Reports the covariates whose posterior inclusion probability exceeds a
#' threshold, the usual summary of a variable selection run.
#'
#' @param object A fitted model.
#' @param threshold Inclusion probability above which a covariate is reported.
#' @param component Optional component to restrict to; see
#'   \code{\link{inclusionProbs}}.
#'
#' @return A named list of character vectors, one per component. For
#'   multi-species models, a list of lists, one per species.
#'
#' @export
medianModel <- function(object, threshold = 0.5, component = NULL) {
  ip <- inclusionProbs(object, component)
  if (!is.list(ip)) ip <- stats::setNames(list(ip), component %||% "component")
  lapply(ip, function(v) {
    if (is.matrix(v)) {
      out <- apply(v, 1L, function(r) colnames(v)[which(r > threshold)])
      if (!is.list(out)) out <- as.list(as.data.frame(out, stringsAsFactors = FALSE))
      out
    } else {
      names(v)[which(v > threshold)]
    }
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a
