#' bvsOccupancy: Bayesian Variable Selection for Occupancy Models
#'
#' Fits single-season, dynamic (multi-season) and multi-species occupancy
#' models with Bayesian variable selection on every component of the model.
#' Occupancy, detection, persistence and colonisation each carry their own
#' vector of binary inclusion indicators, and the sampler visits models of
#' different dimension by adding, deleting and swapping covariates.
#'
#' @section Why Polya-Gamma augmentation:
#' Writing the Bernoulli likelihoods in Polya-Gamma augmented form makes the
#' regression coefficients conditionally Gaussian. Two things follow. The
#' coefficients can be drawn in one conjugate block rather than tuned, and,
#' more importantly, they can be integrated out analytically when the inclusion
#' vector is updated. The add / delete / swap move therefore compares models on
#' their marginal likelihood without any reversible-jump machinery, which is
#' what keeps mixing usable once the number of candidate covariates grows.
#'
#' @section Categorical covariates:
#' A factor enters the design matrix as several dummy columns, and selecting
#' those columns independently would make the answer depend on which level
#' happens to be the reference. Throughout the package, columns are grouped by
#' \code{indexes_covariates} and a covariate is included or excluded as a whole
#' block. The prior covariance gives the dummy columns of a group an
#' exchangeable correlation of one half, which makes the species-level and
#' community-level effects invariant to the reference level.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{bvsSSOM}}}{Single-season, single-species occupancy.}
#'   \item{\code{\link{bvsSDOM}}}{Single-species dynamic occupancy, with
#'     persistence and colonisation.}
#'   \item{\code{\link{bvsMSOM}}}{Multi-species occupancy with community-level
#'     means and variances.}
#'   \item{\code{\link{bvsDesign}}, \code{\link{bvsDesignArray}}}{Build design
#'     matrices and the covariate grouping from a model formula.}
#'   \item{\code{\link{bvsPrior}}}{Construct the prior on the regression
#'     coefficients.}
#'   \item{\code{\link{inclusionProbs}}}{Posterior inclusion probabilities.}
#'   \item{\code{\link{simSSOM}}, \code{\link{simSDOM}},
#'     \code{\link{simMSOM}}}{Simulate data from each model.}
#' }
#'
#' @references
#' Polson, N. G., Scott, J. G. and Windle, J. (2013) Bayesian inference for
#' logistic models using Polya-Gamma latent variables.
#' \emph{Journal of the American Statistical Association} \strong{108}, 1339--1349.
#'
#' MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle, J. A.
#' and Langtimm, C. A. (2002) Estimating site occupancy rates when detection
#' probabilities are less than one. \emph{Ecology} \strong{83}, 2248--2255.
#'
#' MacKenzie, D. I., Nichols, J. D., Hines, J. E., Knutson, M. G. and Franklin,
#' A. B. (2003) Estimating site occupancy, colonization, and local extinction
#' when a species is detected imperfectly. \emph{Ecology} \strong{84}, 2200--2207.
#'
#' @keywords internal
"_PACKAGE"

#' @import nimble
#' @importFrom stats model.matrix model.frame terms na.pass cov var sd quantile
#'   median rbinom rnorm runif rexp plogis setNames delete.response reformulate
#' @importFrom utils modifyList
#' @importFrom coda mcmc mcmc.list effectiveSize
NULL
