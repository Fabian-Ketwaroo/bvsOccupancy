## Variable selection samplers for the occupancy (psi) component.


#' Variable selection sampler for occupancy
#'
#' Add / delete / swap Metropolis-Hastings update of the inclusion vector for the occupancy component, with the coefficients integrated out. The intercept is always retained, and the dummy columns of a categorical covariate move as a block.
#'
#' @section Use:
#' This sampler is assigned automatically by \code{\link{bvsSSOM}},
#' \code{\link{bvsSDOM}} and \code{\link{bvsMSOM}}. It is exported so that
#' it can be added by hand to a custom \pkg{nimble} model configuration with
#' \code{conf$addSampler()}.
#'
#' @param model The model object.
#' @param mvSaved A \code{modelValues} object used to restore the model state.
#' @param target The node being sampled.
#' @param control A named list of control parameters:
#'   \describe{
#'   \item{ncov}{Number of covariates excluding the intercept.}
#'   \item{indexes_covariates}{Integer vector mapping columns of \code{designMatrix} to covariate groups.}
#'   \item{designMatrix}{Design matrix, with rows aligned to the nodes named below.}
#'   \item{PG}{Name of the Polya-Gamma node.}
#'   \item{fixedEffects}{Name of the regression coefficient node.}
#'   }
#'
#' @return A \pkg{nimble} sampler object.
#' @export
gamma_sampler_psi <- nimbleFunction(
  name = 'gamma_sampler_psi',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    X <- control$designMatrix
    N <- dim(X)[1]; M <- dim(X)[2]
    fixedIndexes <- 1
    ncov <-  if(!is.null(control$ncov)) control$ncov else 7
    d_bar <- ifelse(ncov <= 2, 1, 2)#2#ncov  
    d_bar <- ifelse(d_bar <= ncov, d_bar, 2)
    if (is.null(control$indexes_covariates)) {
      stop("`indexes_covariates` must be supplied in `control` for gamma_sampler_psi.", call. = FALSE)
    }
    indexes_covariates <- as.numeric(control$indexes_covariates)
    znodes <- model$getDependencies(target, stochOnly = TRUE, self = FALSE)
    PG_nodes <- model$expandNodeNames(control$PG, returnScalarComponents = TRUE)
    fixedeffects <- control$fixedEffects
  },
  run = function() {
    
    ## Get current values.
    # Needs to be update everything as it is latent
    b <- model$getParam(fixedeffects, 'mean') #model$beta.comm
    B <- model$getParam(fixedeffects , 'cov') #model$var.comm
    
    z <- values(model, znodes)
    Kst <- rep(1, N)  
    k <- z[1:N] - (Kst[1:N]/2) 
    
    Omega <- values(model, PG_nodes)#model$Omega
    gamma <- values(model, target) #model$gamma
    gamma_star <- gamma
    h_ratio <- 1
    M1 <- ncov+1
    
    if( runif(1,0, 1) < 0.33333) { # add
      
      if( (sum(gamma[1:M1]) - fixedIndexes) != ncov ){
        
        
        # Find zero covariates
        numOfZeroCov = ncov - ( sum(gamma[1:M1])  - fixedIndexes ) # correct to here
        zeroCov = nimNumeric(numOfZeroCov)
        
        i <- 1
        for (l in (fixedIndexes+1):length(gamma)) {
          if(gamma[l] == 0){
            zeroCov[i] = l
            i = i+ 1
          }
          
        }
        
        covariate_to_update = zeroCov[rcat(1, rep(1/numOfZeroCov, numOfZeroCov))]
        
        
        # covariate_to_update =  sample_int(zeroCov[1:numOfZeroCov])
        
        #
        gamma_star[covariate_to_update] = 1
        h_ratio = (ncov -(sum(gamma[1:M1]) - fixedIndexes)) / (ncov - (sum(gamma[1:M1]) - fixedIndexes) - 1 + ((ncov - d_bar) / (d_bar)) )
      }
      
    } else if( runif(1,0, 1) < .5){ # delete
      #
      if( ( sum(gamma[1:M1]) - fixedIndexes) != 0){
        
        # Find non zero covariates
        numOfNonZeroCov = sum(gamma[1:M1]) - fixedIndexes
        nonZeroCov = integer(numOfNonZeroCov)
        
        
        i = 1
        for (l in (fixedIndexes+1):length(gamma) ) {
          if(gamma[l] == 1) {
            nonZeroCov[i] = l
            i = i + 1
          }
          
        }
        
        covariate_to_update = nonZeroCov[rcat(1, rep(1/numOfNonZeroCov, numOfNonZeroCov))]#sample_int(nonZeroCov)
        
        gamma_star[covariate_to_update] = 0
        
        h_ratio = (ncov - (sum(gamma[1:M1]) - fixedIndexes) + ((ncov - d_bar) / (d_bar)) ) / (ncov - (sum(gamma[1:M1]) - fixedIndexes) + 1)
      }
      #
    }  else { # swap
      # 
      if( (sum(gamma[1:M1]) - fixedIndexes) != 0 & (sum(gamma[1:M1]) - fixedIndexes) != ncov ) {
        
        # Find zero covariates
        numOfZeroCov = ncov - (sum(gamma[1:M1]) - fixedIndexes)
        zeroCov = integer(numOfZeroCov)
        
        i = 1
        for (l in (fixedIndexes+1) : length(gamma)) {
          if(gamma[l] == 0){
            zeroCov[i] = l
            i = i+ 1
          }
          
        }
        
        covariates2_to_swap = zeroCov[rcat(1, rep(1/numOfZeroCov, numOfZeroCov))] #sample_int(zeroCov)
        
        # Find non zero covariates
        numOfNonZeroCov = sum(gamma[1:M1]) - fixedIndexes
        nonZeroCov = integer(numOfNonZeroCov)
        
        i = 1
        for (l in (fixedIndexes+1) : length(gamma) ) {
          if(gamma[l] == 1) {
            nonZeroCov[i] = l
            i = i + 1
          }
          
        }
        
        covariates1_to_swap = nonZeroCov[rcat(1, rep(1/numOfNonZeroCov, numOfNonZeroCov))]#sample_int(nonZeroCov)
        
        gamma_star[covariates1_to_swap] = 0
        gamma_star[covariates2_to_swap] = 1
        
        h_ratio = 1
        
      }
      # 
    }
    
    model[[target]][1:M1] <<- gamma_star[1:M1] # update proposal
    
    # Compute log likelihood of proposed and current value
    L_gamma_star = dLgamma(gamma_star[1:M1], X[1:N,1:M], indexes_covariates[1:M], b[1:M], B[1:M,1:M], Omega[1:N], k[1:N])
    
    L_gamma = dLgamma(gamma[1:M1], X[1:N,1:M], indexes_covariates[1:M], b[1:M], B[1:M,1:M], Omega[1:N], k[1:N])
    
    h_ratio = h_ratio * exp(L_gamma_star - L_gamma) # MH ratio
    
    if( runif(1,0, 1) < h_ratio ){
      jump <- TRUE
    }  else jump <- FALSE
    
    # keep the model and mvSaved objects consistent
    if(jump) copy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
    else copy(from = mvSaved, to = model, row = 1,  nodes = target, logProb = TRUE)
    
  },
  methods = list(
    reset = function() { }
  )
)


#' Polya-Gamma update for occupancy
#'
#' Draws the Polya-Gamma latent variables of the occupancy component from their full conditional distribution, one per site.
#'
#' @section Use:
#' This sampler is assigned automatically by \code{\link{bvsSSOM}},
#' \code{\link{bvsSDOM}} and \code{\link{bvsMSOM}}. It is exported so that
#' it can be added by hand to a custom \pkg{nimble} model configuration with
#' \code{conf$addSampler()}.
#'
#' @param model The model object.
#' @param mvSaved A \code{modelValues} object used to restore the model state.
#' @param target The node being sampled.
#' @param control A named list of control parameters:
#'   \describe{
#'   \item{indexes_covariates}{Integer vector mapping columns of \code{designMatrix} to covariate groups.}
#'   \item{designMatrix}{Design matrix, with rows aligned to the nodes named below.}
#'   \item{fixedEffects}{Name of the regression coefficient node.}
#'   \item{Gamma}{Name of the inclusion indicator node.}
#'   }
#'
#' @return A \pkg{nimble} sampler object.
#' @export
PG_sampler_psi <- nimbleFunction(
  name = 'PG_sampler_psi',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    X <- control$designMatrix
    nsize = dim(X)[1] # assume that all observations have the same number of samples
    m <- dim(X)[2]
    indexes_covariates <- as.numeric(if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:m)
    gamma_nodes <- model$expandNodeNames(control$Gamma, returnScalarComponents = TRUE) 
    beta_nodes <- model$expandNodeNames(control$fixedEffects, returnScalarComponents = TRUE) 
  },
  run = function() {
    
    n <- rep(1, nsize)  
    gamma <- values(model, gamma_nodes)#model$gamma
    beta <- values(model, beta_nodes)#model$beta
    
    # resize beta and X 
    l <- 0
    index_present <- integer(length(gamma))
    
    for (i in 1:length(gamma)) {
      if(gamma[indexes_covariates[i]] == 1){
        index_present[l+1] <- i
        l <- l + 1
      }
    }
    
    X_gamma <- matrix(0, nsize,l)
    beta_gamma <- numeric(l)
    
    for (i in 1:l) {
      X_gamma[,i] <- X[, index_present[i]]
      beta_gamma[i] <- beta[index_present[i]]
    }
    
    model[[target]][1:nsize] <<- sample_Omega(X_gamma[1:nsize, 1:l], beta_gamma[1:l], n[1:nsize]) # this matches the result obtain from Alex's function
    
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
  },
  methods = list(
    reset = function() { }
  )
)



#' Conjugate coefficient update for occupancy
#'
#' Draws the active occupancy coefficients from their Gaussian full conditional distribution given the Polya-Gamma latent variables. Coefficients of excluded covariates are left untouched.
#'
#' @section Use:
#' This sampler is assigned automatically by \code{\link{bvsSSOM}},
#' \code{\link{bvsSDOM}} and \code{\link{bvsMSOM}}. It is exported so that
#' it can be added by hand to a custom \pkg{nimble} model configuration with
#' \code{conf$addSampler()}.
#'
#' @param model The model object.
#' @param mvSaved A \code{modelValues} object used to restore the model state.
#' @param target The node being sampled.
#' @param control A named list of control parameters:
#'   \describe{
#'   \item{indexes_covariates}{Integer vector mapping columns of \code{designMatrix} to covariate groups.}
#'   \item{designMatrix}{Design matrix, with rows aligned to the nodes named below.}
#'   \item{Gamma}{Name of the inclusion indicator node.}
#'   \item{PG}{Name of the Polya-Gamma node.}
#'   }
#'
#' @return A \pkg{nimble} sampler object.
#' @export
beta_sampler_psi <- nimbleFunction(
  name = 'beta_sampler_psi',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    #calcNodes <- model$getDependencies(target) # all the betas
    X <- control$designMatrix
    N <- dim(X)[1]; M <- dim(X)[2]
    indexes_covariates <- as.numeric(if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:M)
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    znodes <- model$getDependencies(target, stochOnly = TRUE, self = FALSE)
    
    calcNodes <- model$getDependencies(target)
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    gamma_nodes <- model$expandNodeNames(control$Gamma, returnScalarComponents = TRUE) 
    PG_nodes <- model$expandNodeNames(control$PG, returnScalarComponents = TRUE)
  },
  run = function() {
    
    b <- model$getParam(target, 'mean') #model$beta.comm
    B <- model$getParam(target, 'cov') #model$var.comm
    
    z <- values(model, znodes)
    Kst <- rep(1, N)  #rep(model$K, N)
    k <- z[1:N] - (Kst[1:N]/2) # works with K fixed
    #  k <- model$y - (21/2) # works with K fixed
    gamma <- values(model, gamma_nodes) #model$gamma
    Omega <- values(model, PG_nodes)
    
    # resize X, b and B
    l = 0
    index_present = integer(length(indexes_covariates))
    
    # gamma <- c(1,1,1,0,0)
    # indexes_covariates <- 1:5
    for (i in 1:length(indexes_covariates)) {
      
      if(gamma[indexes_covariates[i]] == 1){
        index_present[l+1] = i
        l = l + 1
      }
    }
    
    X_gamma2 = matrix(0, dim(X)[1], l)
    b_gamma = numeric(l)
    B_gamma = matrix(0, l,l)
    
    for (i in 1:l) {
      X_gamma2[,i] = X[, index_present[i]]
      b_gamma[i] = b[index_present[i]]
      for (j in 1:l) {
        B_gamma[i,j] = B[index_present[i], index_present[j]]
      }
    }
    
    # sample beta
    beta_new  <- sample_beta(X_gamma2[1:N,1:l], B_gamma[1:l,1:l], b_gamma[1:l], Omega[1:N], k[1:N])
    #model[[target]][1:l] <<- sample_beta(X_gamma2[1:N,1:l], B_gamma[1:l,1:l], b_gamma[1:l], Omega[1:N], k[1:N]) # when l not equal to length of gamma, I'm updating the wrong gamma
    
    for (i in 1:l) {
      # beta[index_present[i]] = beta_new[i]
      model[[target]][index_present[i]] <<- beta_new[i]  
    }
    
    update_other_nodes <- model$calculateDiff(calcNodesDepStage)
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
    
  },
  methods = list(
    reset = function() { }
  )
)
