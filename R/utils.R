## Internal helpers.
##
## The samplers work on flattened vectors: `values(model, ynodes)` is compared
## with `values(model, znodes)`, and rows of the flattened design matrix must
## line up with both. Rather than assume nimble expands node names in a
## particular order, everything here is derived from the node names the model
## actually reports, so the alignment holds whatever nimble does internally.

## "y[3, 2, 7]" -> c(3, 2, 7), returned as an integer matrix with one row per node.
.parseNodeIndices <- function(nodes) {
  if (!length(nodes)) {
    stop("no nodes to parse; the model does not contain the expected variable.", call. = FALSE)
  }
  inside <- sub("\\]\\s*$", "", sub("^[^\\[]*\\[", "", nodes))
  parts <- strsplit(inside, ",", fixed = TRUE)
  k <- length(parts[[1L]])
  if (any(lengths(parts) != k)) {
    stop("nodes have inconsistent numbers of indices; expandNodeNames() returned something unexpected.",
         call. = FALSE)
  }
  out <- matrix(suppressWarnings(as.integer(trimws(unlist(parts, use.names = FALSE)))),
                ncol = k, byrow = TRUE)
  if (anyNA(out)) {
    stop("could not parse node indices; use returnScalarComponents = TRUE when expanding nodes.",
         call. = FALSE)
  }
  out
}

## Position of each row of `idx` within expandNodeNames(var). This is exactly
## what the samplers' `zindex` and `omega_pos` control arguments need.
.linIndex <- function(model, var, idx, what = var) {
  nodes <- model$expandNodeNames(var, returnScalarComponents = TRUE)
  nix <- .parseNodeIndices(nodes)
  if (ncol(nix) != ncol(idx)) {
    stop(sprintf("`%s` has %d indices but %d were supplied.", what, ncol(nix), ncol(idx)),
         call. = FALSE)
  }
  key <- function(m) do.call(paste, c(split(m, col(m)), list(sep = ",")))
  pos <- match(key(idx), key(nix))
  if (anyNA(pos)) {
    stop(sprintf("some observations have no matching node in `%s`; the model and the data disagree.",
                 what), call. = FALSE)
  }
  as.numeric(pos)
}

## The community samplers reach into a flattened vector of species-level nodes
## using the arithmetic (column - 1) * nspecies + species, which is only valid
## if nimble lists the nodes of a matrix variable column by column. Rather than
## take that on trust, confirm it once when the model is built.
.assertColumnMajor <- function(model, var, nr) {
  nix <- .parseNodeIndices(model$expandNodeNames(var, returnScalarComponents = TRUE))
  if (ncol(nix) != 2L) {
    stop(sprintf("`%s` should be a matrix-valued variable.", var), call. = FALSE)
  }
  expected <- (nix[, 2L] - 1L) * nr + nix[, 1L]
  if (!identical(as.integer(expected), seq_len(nrow(nix)))) {
    stop(sprintf("the nodes of `%s` are not listed column by column, which the community ",
                 var),
         "samplers assume. This indicates an unexpected node ordering in nimble; ",
         "please report it as a bug.", call. = FALSE)
  }
  invisible(TRUE)
}

## Check a detection array and report its non-missing cells.
.checkY <- function(y, ndim, argname = "y") {
  if (is.null(dim(y)) || length(dim(y)) != ndim) {
    stop(sprintf("`%s` must be a %d-dimensional array.", argname, ndim), call. = FALSE)
  }
  ok <- is.na(y) | y == 0 | y == 1
  if (!all(ok)) {
    stop(sprintf("`%s` must contain only 0, 1 and NA.", argname), call. = FALSE)
  }
  storage.mode(y) <- "double"
  y
}

## Move the surveyed visits of each site to the front, permuting the detection
## covariates the same way. Visits within a site are exchangeable given their
## covariates, so this changes nothing about the model; it just lets the
## likelihood be declared as `for (j in 1:Jindex[s])`, which keeps unsurveyed
## cells out of the model graph entirely.
.leftPack <- function(y, X, visitDim = 2L) {
  d <- dim(y)
  changed <- FALSE
  if (length(d) == 2L) {
    for (i in seq_len(d[1L])) {
      o <- order(is.na(y[i, ]))
      if (!identical(o, seq_len(d[2L]))) {
        changed <- TRUE
        y[i, ] <- y[i, o]
        if (!is.null(X)) X[i, , ] <- X[i, , o, drop = FALSE]
      }
    }
  } else if (length(d) == 3L) {
    for (i in seq_len(d[1L])) {
      for (t in seq_len(d[3L])) {
        o <- order(is.na(y[i, , t]))
        if (!identical(o, seq_len(d[2L]))) {
          changed <- TRUE
          y[i, , t] <- y[i, o, t]
          if (!is.null(X)) X[i, , , t] <- X[i, , o, t, drop = FALSE]
        }
      }
    }
  } else {
    stop("`.leftPack` handles two- and three-dimensional detection arrays only.", call. = FALSE)
  }
  list(y = y, X = X, changed = changed)
}

## Draw from PG(1, 0), used only to initialise the latent variables. The
## representation is an infinite weighted sum of exponentials; truncating it is
## more than accurate enough for a starting value, and avoids a dependency.
.rpg0 <- function(n, terms = 500L) {
  d <- (seq_len(terms) - 0.5)^2
  vapply(seq_len(n), function(i) sum(stats::rexp(terms) / d) / (2 * pi^2), numeric(1))
}

## Flatten a design array so that row r corresponds to node r of `nodes`.
## `idx` gives, for each node, the index into each dimension of the array;
## `map` says which columns of `idx` map to which array dimension.
.flattenDesign <- function(X, idx, map) {
  V <- dim(X)[2L]
  n <- nrow(idx)
  out <- matrix(0, n, V)
  for (v in seq_len(V)) {
    sel <- cbind(idx[, map[1L]], rep.int(v, n))
    if (length(map) > 1L) {
      for (m in map[-1L]) sel <- cbind(sel, idx[, m])
    }
    out[, v] <- X[sel]
  }
  out
}

## Split a design array into a single slice, keeping it a matrix.
.slice <- function(X, ...) {
  idx <- list(...)
  args <- c(list(X, TRUE, TRUE), idx, list(drop = FALSE))
  s <- do.call("[", args)
  matrix(as.vector(s), nrow = dim(X)[1L], ncol = dim(X)[2L])
}

## ------------------------------------------------------------- summaries ---

## Split-chain potential scale reduction factor. Returns NA for parameters that
## never move, which is the normal state of an intercept indicator fixed at 1.
.rhat <- function(x) {
  m <- ncol(x)
  n <- nrow(x)
  if (m < 2L || n < 2L) return(NA_real_)
  W <- mean(apply(x, 2L, stats::var))
  if (!is.finite(W) || W <= .Machine$double.eps) return(NA_real_)
  B <- n * stats::var(colMeans(x))
  sqrt(((n - 1) / n * W + B / n) / W)
}

.mcmcSummary <- function(samples) {
  if (inherits(samples, "mcmc.list")) {
    chains <- lapply(samples, as.matrix)
  } else {
    chains <- list(as.matrix(samples))
  }
  all <- do.call(rbind, chains)
  pars <- colnames(all)
  ess <- tryCatch(as.numeric(coda::effectiveSize(samples)), error = function(e) rep(NA_real_, length(pars)))
  rh <- vapply(seq_along(pars), function(j) {
    .rhat(vapply(chains, function(ch) ch[, j], numeric(nrow(chains[[1L]]))))
  }, numeric(1))
  q <- t(apply(all, 2L, stats::quantile, probs = c(0.025, 0.5, 0.975), names = FALSE))
  out <- data.frame(
    mean = colMeans(all),
    sd = apply(all, 2L, stats::sd),
    `2.5%` = q[, 1L],
    `50%` = q[, 2L],
    `97.5%` = q[, 3L],
    Rhat = round(rh, 3),
    n.eff = round(ess),
    check.names = FALSE,
    row.names = pars
  )
  out
}

## Per-chain initial values. Chains after the first are jittered so that
## convergence diagnostics have something to work with.
.makeInits <- function(base, n.chains, user = NULL, jitter = 0.25) {
  if (!is.null(user)) {
    if (is.list(user) && length(user) && is.list(user[[1L]])) {
      if (length(user) != n.chains) {
        stop("`inits` given as a list of lists must have one element per chain.", call. = FALSE)
      }
      return(lapply(user, function(u) utils::modifyList(base, u)))
    }
    base <- utils::modifyList(base, user)
  }
  out <- vector("list", n.chains)
  out[[1L]] <- base
  if (n.chains > 1L) {
    for (k in 2:n.chains) {
      ini <- base
      for (nm in names(ini)) {
        if (grepl("^beta", nm) && is.numeric(ini[[nm]])) {
          ini[[nm]] <- ini[[nm]] + stats::rnorm(length(ini[[nm]]), 0, jitter)
        }
        if (grepl("^Omega", nm) && is.numeric(ini[[nm]])) {
          ini[[nm]][] <- .rpg0(length(ini[[nm]]))
        }
      }
      out[[k]] <- ini
    }
  }
  out
}

.tick <- function(verbose, ...) {
  if (isTRUE(verbose)) message(...)
  invisible(NULL)
}
