## Design matrices and the covariate grouping that variable selection acts on.

#' Build a design matrix for variable selection
#'
#' Creates the design matrix for one component of an occupancy model, together
#' with the grouping vector that tells the samplers which columns belong to
#' which covariate. A factor contributes several dummy columns, and those
#' columns are selected or dropped as a single block so that the answer does
#' not depend on which level happens to be the reference.
#'
#' The grouping is read straight off the \code{"assign"} attribute of
#' \code{\link[stats]{model.matrix}}, so interactions and polynomial terms are
#' grouped correctly as well.
#'
#' @param formula A one-sided model formula, for example \code{~ elev + habitat}.
#'   An intercept is required: it is always retained by the sampler and acts as
#'   the reference model.
#' @param data A data frame containing the variables in \code{formula}. One row
#'   per site.
#' @param contrasts Optional list of contrasts, passed to
#'   \code{\link[stats]{model.matrix}}. Treatment contrasts are recommended,
#'   since the prior is built for dummy coding.
#'
#' @return An object of class \code{"bvsDesign"}: a list with components
#'   \describe{
#'     \item{X}{The design matrix, with the intercept in the first column.}
#'     \item{indexes_covariates}{Integer vector, one entry per column of
#'       \code{X}, giving the covariate group. The intercept is group 1.}
#'     \item{ncov}{Number of candidate covariates, excluding the intercept.}
#'     \item{numVars}{Number of columns of \code{X}.}
#'     \item{labels}{Covariate labels, one per group.}
#'   }
#'
#' @seealso \code{\link{bvsDesignArray}} for covariates that vary within a
#'   site, and \code{\link{bvsPrior}} for the matching prior.
#'
#' @examples
#' site <- data.frame(elev = rnorm(50), habitat = factor(sample(1:3, 50, TRUE)))
#' d <- bvsDesign(~ elev + habitat, site)
#' d$indexes_covariates   # 1 2 3 3: habitat's two dummies form one group
#'
#' @export
bvsDesign <- function(formula, data, contrasts = NULL) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula, for example ~ elev + habitat.", call. = FALSE)
  }
  if (length(formula) == 3L) formula <- stats::delete.response(stats::terms(formula))
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  tt <- attr(mf, "terms")
  X <- stats::model.matrix(tt, mf, contrasts.arg = contrasts)
  .finishDesign(X, tt, call = match.call())
}

#' Build an array design matrix for varying covariates
#'
#' Detection covariates usually vary between visits to the same site, and the
#' persistence and colonisation covariates of a dynamic model vary between
#' seasons. This function builds the corresponding design array from data in
#' long format.
#'
#' \code{data} must be sorted so that the first dimension varies fastest: site
#' within visit within season, the ordering produced by
#' \code{expand.grid(site = 1:M, visit = 1:J, season = 1:T)}. Building the
#' whole array from a single call to \code{\link[stats]{model.matrix}}
#' guarantees that every slice gets the same columns even when a factor level
#' is absent from some visits.
#'
#' @inheritParams bvsDesign
#' @param data A data frame with \code{prod(dims)} rows, in the order described
#'   above.
#' @param dims Integer vector of dimensions, excluding the covariate dimension.
#'   Use \code{c(nsites, nvisits)} for a single-season detection design and
#'   \code{c(nsites, nvisits, nseasons)} for a dynamic one.
#'
#' @return An object of class \code{"bvsDesign"} whose \code{X} component is an
#'   array of dimension \code{c(dims[1], numVars, dims[-1])}, matching the
#'   layout the model code expects.
#'
#' @examples
#' M <- 20; J <- 3
#' visits <- data.frame(
#'   wind = rnorm(M * J),
#'   obs  = factor(sample(c("a", "b"), M * J, TRUE))
#' )
#' d <- bvsDesignArray(~ wind + obs, visits, dims = c(M, J))
#' dim(d$X)   # 20 sites x 3 columns x 3 visits
#'
#' @export
bvsDesignArray <- function(formula, data, dims, contrasts = NULL) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  dims <- as.integer(dims)
  if (length(dims) < 2L || anyNA(dims) || any(dims < 1L)) {
    stop("`dims` must be two or more positive integers.", call. = FALSE)
  }
  if (nrow(data) != prod(dims)) {
    stop(sprintf(
      "`data` has %d rows but `dims` implies %d. Rows must be ordered with the first dimension varying fastest.",
      nrow(data), prod(dims)
    ), call. = FALSE)
  }
  if (length(formula) == 3L) formula <- stats::delete.response(stats::terms(formula))
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  tt <- attr(mf, "terms")
  mm <- stats::model.matrix(tt, mf, contrasts.arg = contrasts)
  V <- ncol(mm)

  ## as.vector(mm) runs down the rows first, so reshaping to c(dims, V) puts
  ## the covariate index last; aperm then moves it into position 2.
  arr <- array(as.vector(mm), dim = c(dims, V))
  perm <- c(1L, length(dims) + 1L, seq.int(2L, length(dims)))
  X <- aperm(arr, perm)

  out <- .finishDesign(mm, tt, call = match.call())
  out$X <- X
  out$dims <- dim(X)
  out
}

#' Wrap an existing design matrix
#'
#' Use this when the design matrix was built elsewhere, for instance by code
#' that accompanies a published analysis.
#'
#' @param X A design matrix or design array. The intercept must be the first
#'   column, and for an array the covariate index must be the second dimension.
#' @param indexes_covariates Integer vector, one entry per column of \code{X},
#'   giving the covariate group. Must start at 1 for the intercept and increase
#'   in steps of one, with the columns of a group contiguous.
#' @param labels Optional covariate labels, one per group.
#'
#' @return An object of class \code{"bvsDesign"}.
#'
#' @examples
#' X <- cbind(1, matrix(rnorm(50 * 3), 50, 3))
#' as.bvsDesign(X, c(1, 2, 3, 4))
#'
#' @export
as.bvsDesign <- function(X, indexes_covariates, labels = NULL) {
  if (inherits(X, "bvsDesign")) return(X)
  d <- dim(X)
  if (is.null(d) || length(d) < 2L) {
    stop("`X` must be a matrix or an array.", call. = FALSE)
  }
  V <- d[2L]
  ic <- .checkIndexes(indexes_covariates, V)
  ncov <- as.integer(max(ic)) - 1L
  if (is.null(labels)) {
    labels <- c("(Intercept)", paste0("x", seq_len(ncov)))
  }
  if (length(labels) != ncov + 1L) {
    stop(sprintf("`labels` must have %d entries, one per covariate group including the intercept.",
                 ncov + 1L), call. = FALSE)
  }
  structure(
    list(X = X, indexes_covariates = ic, ncov = ncov, numVars = V,
         labels = labels, dims = d, terms = NULL, call = match.call()),
    class = "bvsDesign"
  )
}

#' Prior for the regression coefficients
#'
#' Builds the multivariate normal prior \eqn{\beta \sim N(b, B)} used by every
#' component of the model. The intercept gets a diffuse variance of its own.
#' The remaining columns get a block diagonal covariance: a continuous
#' covariate contributes a single variance, while the dummy columns of a
#' categorical covariate share an exchangeable block with correlation one half.
#' That block structure is what makes the fitted effects invariant to the
#' choice of reference level.
#'
#' @param design A \code{"bvsDesign"} object, or an integer vector of covariate
#'   group indices.
#' @param mu Prior mean for the intercept, on the logit scale. The default of
#'   0.5 corresponds to an occupancy or detection probability slightly above
#'   one half.
#' @param phi.mu Prior variance for the intercept.
#' @param phi.beta Scale of the prior covariance for the remaining
#'   coefficients. Smaller values shrink harder and favour smaller models.
#'
#' @return A list with the prior mean \code{b}, the prior covariance \code{B},
#'   and the unscaled block correlation matrix \code{C}.
#'
#' @examples
#' site <- data.frame(elev = rnorm(50), habitat = factor(sample(1:3, 50, TRUE)))
#' d <- bvsDesign(~ elev + habitat, site)
#' p <- bvsPrior(d)
#' round(p$B, 3)
#'
#' @export
bvsPrior <- function(design, mu = 0.5, phi.mu = 4, phi.beta = 1 / 4) {
  ic <- if (inherits(design, "bvsDesign")) design$indexes_covariates else
    .checkIndexes(design, length(design))
  V <- length(ic)
  if (!is.numeric(mu) || length(mu) != 1L) stop("`mu` must be a single number.", call. = FALSE)
  if (phi.mu <= 0 || phi.beta <= 0) stop("`phi.mu` and `phi.beta` must be positive.", call. = FALSE)

  C <- matrix(0, nrow = max(V - 1L, 0L), ncol = max(V - 1L, 0L))
  if (V > 1L) {
    l <- 0L
    for (g in 2:max(ic)) {
      L <- sum(ic == g)
      if (L == 1L) {
        C[l + 1L, l + 1L] <- 1
      } else {
        C[l + seq_len(L), l + seq_len(L)] <- 0.5 * (diag(L) + matrix(1, L, L))
      }
      l <- l + L
    }
  }

  b <- c(mu, rep(0, V - 1L))
  B <- matrix(0, V, V)
  B[1L, 1L] <- phi.mu
  if (V > 1L) B[-1L, -1L] <- C * phi.beta

  list(b = b, B = B, C = C, mu = mu, phi.mu = phi.mu, phi.beta = phi.beta)
}

#' @export
print.bvsDesign <- function(x, ...) {
  cat("<bvsDesign>\n")
  cat("  dimensions        :", paste(x$dims, collapse = " x "), "\n")
  cat("  columns           :", x$numVars, "\n")
  cat("  candidate covars  :", x$ncov, "\n")
  grp <- split(seq_len(x$numVars), x$indexes_covariates)
  lab <- x$labels
  cat("  groups            :\n")
  for (g in seq_along(grp)) {
    cat(sprintf("    %2d  %-24s columns %s\n", g,
                if (g <= length(lab)) lab[g] else paste0("group", g),
                paste(grp[[g]], collapse = ",")))
  }
  invisible(x)
}

## ---------------------------------------------------------------- internal --

## Validate a user supplied prior against the design it belongs to. Returns the
## prior unchanged so it can be used inline.
.checkPrior <- function(prior, design, name) {
  if (!is.list(prior) || !all(c("b", "B") %in% names(prior))) {
    stop(sprintf("`%s` must be a list with components `b` and `B`, as returned by bvsPrior().",
                 name), call. = FALSE)
  }
  V <- design$numVars
  if (!is.numeric(prior$b) || length(prior$b) != V) {
    stop(sprintf("`%s` gives a prior mean of length %d but the design has %d columns.",
                 name, length(prior$b), V), call. = FALSE)
  }
  if (!is.matrix(prior$B) || nrow(prior$B) != V || ncol(prior$B) != V) {
    stop(sprintf("`%s` gives a %s prior covariance but the design has %d columns.",
                 name, paste(dim(as.matrix(prior$B)), collapse = " x "), V), call. = FALSE)
  }
  if (!isTRUE(all.equal(unname(prior$B), unname(t(prior$B))))) {
    stop(sprintf("the prior covariance in `%s` is not symmetric.", name), call. = FALSE)
  }
  ev <- tryCatch(min(eigen(prior$B, symmetric = TRUE, only.values = TRUE)$values),
                 error = function(e) NA_real_)
  if (is.na(ev) || ev <= 0) {
    stop(sprintf("the prior covariance in `%s` is not positive definite.", name),
         call. = FALSE)
  }
  prior
}


.finishDesign <- function(X, tt, call) {
  asn <- attr(X, "assign")
  if (is.null(asn)) stop("could not read the term assignment from the design matrix.", call. = FALSE)
  if (!any(asn == 0L)) {
    stop("the formula must include an intercept: it is always retained and defines the null model.",
         call. = FALSE)
  }
  if (anyNA(X)) {
    stop("the design matrix contains missing values. Covariates must be complete; ",
         "impute them or drop the affected rows before fitting.", call. = FALSE)
  }
  ic <- as.integer(asn) + 1L
  ic <- .checkIndexes(ic, ncol(X))
  labels <- c("(Intercept)", attr(tt, "term.labels"))
  structure(
    list(X = X, indexes_covariates = ic, ncov = as.integer(max(ic)) - 1L, numVars = ncol(X),
         labels = labels, dims = dim(X), terms = tt, call = call),
    class = "bvsDesign"
  )
}

.checkIndexes <- function(ic, V) {
  ic <- as.integer(ic)
  if (length(ic) != V) {
    stop(sprintf("`indexes_covariates` has length %d but the design matrix has %d columns.",
                 length(ic), V), call. = FALSE)
  }
  if (anyNA(ic) || ic[1L] != 1L) {
    stop("`indexes_covariates` must start at 1, the group of the intercept.", call. = FALSE)
  }
  if (is.unsorted(ic)) {
    stop("`indexes_covariates` must be non-decreasing: the columns of a covariate must be contiguous.",
         call. = FALSE)
  }
  if (!identical(sort(unique(ic)), seq_len(max(ic)))) {
    stop("`indexes_covariates` must use every group from 1 to its maximum, with no gaps.",
         call. = FALSE)
  }
  if (sum(ic == 1L) != 1L) {
    stop("group 1 must contain only the intercept column.", call. = FALSE)
  }
  ## Returned as double, not integer, on purpose. This vector is handed to
  ## nimbleFunctions whose run signatures declare `indexes_covariates` (and
  ## `group_id`) as double(1). An integer vector compiles to NimArr<1,int> and
  ## the generated C++ then fails to match the double signature, so the model
  ## builds in R and only breaks at compileNimble(). Validation above is done
  ## on the integer copy; only the storage type changes here.
  as.numeric(ic)
}

## Coerce a user-supplied design and check its shape.
.asDesign <- function(design, argname, ndim, indexes = NULL) {
  if (!inherits(design, "bvsDesign")) {
    if (is.null(indexes)) {
      stop(sprintf("`%s` must be a bvsDesign object; build one with bvsDesign() or as.bvsDesign().",
                   argname), call. = FALSE)
    }
    design <- as.bvsDesign(design, indexes)
  }
  if (length(dim(design$X)) != ndim) {
    stop(sprintf("`%s` must have %d dimensions but has %d. %s",
                 argname, ndim, length(dim(design$X)),
                 if (ndim > 2L) "Use bvsDesignArray() for covariates that vary within a site."
                 else "Use bvsDesign() for site-level covariates."),
         call. = FALSE)
  }
  if (design$ncov < 1L) {
    stop(sprintf("`%s` has no candidate covariates, so there is nothing to select.", argname),
         call. = FALSE)
  }
  if (anyNA(design$X)) {
    stop(sprintf("`%s` contains missing values; covariates must be complete.", argname),
         call. = FALSE)
  }
  design
}
