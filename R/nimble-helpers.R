## Shared nimbleFunction helpers: linear algebra, the marginal likelihood of an
## inclusion vector, and the linear predictor under variable selection.


#' Right-multiply a matrix by a diagonal matrix
#'
#' Computes \code{X \%*\% diag(D)} without forming the diagonal matrix. Used
#' internally when building the Polya-Gamma weighted cross-product
#' \eqn{X' \Omega X}.
#'
#' @param X A numeric matrix.
#' @param D A numeric vector holding the diagonal entries.
#'
#' @return A numeric matrix with the same dimensions as \code{X}.
#' @keywords internal
#' @export
diagMatrixProd <- nimbleFunction(
  run = function( X = double(2), D= double(1) ){
    returnType(double(2))
    
    result <- matrix(0, dim(X)[1], length(D))
    
    for (i in 1:dim(result)[1]) {
      for (j in 1:dim(result)[2]) {
        result[i, j] = X[i,j] * D[j]
      }
      
    }
    
    return(result)
  }
)


#' Marginal log-likelihood of an inclusion vector
#'
#' Evaluates the Polya-Gamma augmented marginal likelihood of a binary
#' inclusion vector \code{gamma} with the regression coefficients integrated
#' out analytically. This is the quantity that drives the Metropolis-Hastings
#' acceptance ratio in the add / delete / swap samplers
#' (\code{\link{gamma_sampler_psi}} and friends).
#'
#' Columns of the design matrix are grouped by \code{indexes_covariates}, so a
#' categorical covariate is either retained or dropped as a block, together
#' with all of its dummy columns.
#'
#' @param x Binary inclusion vector, one entry per covariate group (the
#'   intercept occupies the first position and is always retained).
#' @param X Design matrix restricted to the rows that contribute to the
#'   likelihood.
#' @param indexes_covariates Integer vector mapping each column of \code{X} to
#'   its covariate group; see \code{\link{bvsDesign}}.
#' @param b Prior mean vector for the regression coefficients.
#' @param B Prior covariance matrix for the regression coefficients.
#' @param Omega Current Polya-Gamma latent variables, one per row of \code{X}.
#' @param k Working response, \eqn{y - 1/2}.
#' @param log Logical (as an integer); return the log-likelihood. Defaults to
#'   \code{TRUE}.
#'
#' @return A scalar (log-)likelihood.
#' @references
#' Polson, N. G., Scott, J. G. and Windle, J. (2013) Bayesian inference for
#' logistic models using Polya-Gamma latent variables.
#' \emph{Journal of the American Statistical Association} \strong{108}, 1339--1349.
#' @keywords internal
#' @export
dLgamma <- nimbleFunction(
  run= function( x = double(1) , X= double(2),  indexes_covariates = double(1), b = double(1),  B= double(2),  Omega= double(1),  k= double(1), log = integer(0, default= 1) ){
    returnType(double(0))
    
    gamma <- x
    m <- length(indexes_covariates)
    index_present <-  numeric(m)
    l = 0
    
    #gamma = rep(1, numVars) # data to test likelihood works correctly
    #  gamma = c(1,1,1,0,0,0)
    
    for (i in 1:m ) {
      if(gamma[indexes_covariates[i]] == 1 ){
        index_present[l+1] = i
        l = l+1
      }
    }
    
    
    X_gamma <- matrix(0, dim(X)[1], l)
    b_gamma <- numeric(l)
    B_gamma <- matrix(0, l, l)
    
    
    for (i in 1:l) {
      X_gamma[,i] = X[, index_present[i]]
      b_gamma[i] = b[index_present[i]]
      for (j in 1:l) {
        B_gamma[i,j] <- B[index_present[i], index_present[j]]
        
      }
    }
    
    tX = t(X_gamma)
    tXOmega = diagMatrixProd(tX, Omega)
    cholXgOmX = chol(tXOmega  %*% X_gamma + inverse(B_gamma) )
    
    firstTerm = (.5) * logdet(inverse(B_gamma)) - logdet(cholXgOmX)
    tXKbplusBb = t(X_gamma) %*% k + inverse(B_gamma) %*% b_gamma  
    v = solve(t(cholXgOmX), tXKbplusBb) # problem with trimatl and trimatl and so I've taken it out 
    vtv = t(v) %*% v
    
    secondTerm = - .5 * ( (t(b_gamma) %*% inverse(B_gamma) %*% b_gamma) - vtv)
    
    loglikelihood <- firstTerm + secondTerm[1,1]
    
    if(log) return(loglikelihood)
    else return(exp(loglikelihood))
    
  }
)


#' Right-multiply a matrix by a diagonal matrix
#'
#' Identical to \code{\link{diagMatrixProd}}; retained as a separate
#' \code{nimbleFunction} so that the coefficient and inclusion samplers can be
#' compiled independently.
#'
#' @inheritParams diagMatrixProd
#' @return A numeric matrix with the same dimensions as \code{X}.
#' @keywords internal
#' @export
dMatrixProd <- nimbleFunction(
  run = function( X = double(2), D= double(1) ){
    returnType(double(2))
    
    result <- matrix(0, dim(X)[1], length(D))
    
    for (i in 1:dim(result)[1]) {
      for (j in 1:dim(result)[2]) {
        result[i, j] = X[i,j] * D[j]
      }
      
    }
    
    return(result)
  }
)


#' Fast multivariate normal draw from a Cholesky factor
#'
#' @param mu Mean vector.
#' @param cholsigma Cholesky factor of the covariance matrix.
#'
#' @return A single draw as a numeric vector.
#' @keywords internal
#' @export
mvrnormQuick <- nimbleFunction(
  run = function(mu = double(1), cholsigma = double(2)){
    returnType(double(1))
    ncols = dim(cholsigma)[2] 
    Y =  rnorm(ncols, 0,1)  
    return (c(mu + cholsigma %*% Y)) # shouldn' this be a matrix?, no it returns a vector in c++
    #return ( Y) # shouldn' this be a matrix? 
  }
)



#' Conjugate update of the regression coefficients
#'
#' Draws the active regression coefficients from their full conditional
#' Gaussian distribution given the Polya-Gamma latent variables.
#'
#' @param X Design matrix restricted to the active columns and to the rows that
#'   contribute to the likelihood.
#' @param B Prior covariance matrix of the active coefficients.
#' @param b Prior mean of the active coefficients.
#' @param Omega Polya-Gamma latent variables, one per row of \code{X}.
#' @param k Working response, \eqn{y - 1/2}.
#'
#' @return A vector of sampled coefficients, one per active column.
#' @keywords internal
#' @export
sample_beta <- nimbleFunction(
  run = function(X= double(2), B= double(2), b = double(1), Omega = double(1) ,k= double(1) ){
    returnType(double(1))
    
    tX <- t(X)
    tXOmega = dMatrixProd(tX, Omega)
    
    L = t(chol(tXOmega %*% X + inverse(B)) )
    
    tmp = solve(L, tX %*% k + inverse(B) %*% b)
    alpha = solve(t(L),tmp) # up to here compiles, # here alpha is a matrix 
    
    result = mvrnormQuick(c(alpha), t(inverse(L))) # here is the problem as alpha is take to be a vector 
    
    
    return(result) 
    
  }
  
)


#' Placeholder density for Polya-Gamma latent variables
#'
#' The Polya-Gamma latent variables are updated in closed form by dedicated
#' samplers (\code{\link{PG_sampler_psi}} and friends), so the model only
#' needs a distribution that contributes a constant to the log posterior.
#' \code{dPG} returns zero on the log scale, and \code{rPG} deliberately
#' stops if it is ever called.
#'
#' The distribution is registered with \pkg{nimble} when \pkg{bvsOccupancy}
#' is loaded, which is what allows model code to contain
#' \code{Omega[1:n] ~ dPG()}.
#'
#' @param x Vector of latent variables.
#' @param n Number of draws; required by \pkg{nimble} but unused.
#' @param log Logical (as an integer); return the value on the log scale.
#'
#' @return \code{dPG} returns \code{0} (or \code{1} when \code{log = FALSE}).
#' @name dPG
#' @aliases rPG
#' @export
dPG <- nimbleFunction(
  run = function(x = double(1), log=integer(0, default=1)){
    returnType(double(0))
    ans <- 0
    if(log) return(ans)
    else return(exp(ans))
  }
)

# some non-significant dist


#' @rdname dPG
#' @export
rPG <- nimbleFunction(
  run = function(n = integer(0)) {
    stop('rPG should never be run')
    returnType(double(1))
    return(0)
  })



#' Linear predictor under the current inclusion vector
#'
#' Forms \eqn{X_\gamma \beta_\gamma}, that is, the linear predictor built
#' from only those covariate groups currently flagged as active. Coefficients
#' of excluded groups are simply not used, so they remain at their last sampled
#' value without influencing the likelihood.
#'
#' @param X Full design matrix.
#' @param beta Full coefficient vector.
#' @param indexes_covariates Integer vector mapping columns of \code{X} to
#'   covariate groups.
#' @param gamma Binary inclusion vector, one entry per covariate group.
#'
#' @return A one-column matrix holding the linear predictor.
#' @export
compute_predictor = nimbleFunction(
  run = function(X = double(2), beta= double(1),  indexes_covariates= double(1), gamma = double(1)){
    returnType(double(2))
    
    L <- length(indexes_covariates)
    index_present <- logical(L)
    
    for (i in 1:L) {
      if(any(indexes_covariates[i] ==which(gamma==1))){
        index_present[i] <- TRUE
      } else index_present[i] <- FALSE
      
    }
    
    X_gamma <- X[,index_present,drop = FALSE]
    beta_gamma <- beta[index_present]
    S <- dim(X_gamma)[1]
    linpred <- X_gamma %*% beta_gamma#[1:S,1]  
    
    return(linpred)
    
  }
)


#' Block community covariance matrix
#'
#' Builds the community-level covariance matrix used by the multi-species
#' model. Each covariate group receives its own variance, and the dummy columns
#' of a categorical covariate share that variance with an exchangeable
#' correlation of one half. This keeps the species-level effects invariant to
#' the choice of reference level.
#'
#' @param tau Vector of group-level variances, one per covariate group.
#' @param group_id Integer vector mapping columns to covariate groups.
#' @param numVars Total number of columns in the design matrix.
#'
#' @return A \code{numVars} by \code{numVars} covariance matrix.
#' @export
build_block_cov <- nimbleFunction(
  run = function(tau = double(1), group_id = double(1), numVars = double(0)) {
    returnType(double(2))
    S <- matrix(0, nrow = numVars, ncol = numVars)
    for (a in 1:numVars) {
      for (b in 1:numVars) {
        if (group_id[a] == group_id[b]) {
          if (a == b) S[a, b] <- tau[group_id[a]]
          else        S[a, b] <- tau[group_id[a]] * 0.5
        }
      }
    }
    return(S)
  }
)
