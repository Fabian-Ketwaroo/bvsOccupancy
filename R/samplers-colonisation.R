## Variable selection samplers for the colonisation (eta) component.


#' Variable selection sampler for colonisation
#'
#' Add / delete / swap Metropolis-Hastings update of the inclusion vector for the colonisation component of a dynamic occupancy model. Only sites that were unoccupied in the current season contribute.
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
#'   \item{znodes}{Name of the latent occupancy node.}
#'   \item{ynodes}{Name of the occupancy node one season ahead. Optional; when omitted the stochastic dependencies of \code{target} are used instead.}
#'   }
#'
#' @return A \pkg{nimble} sampler object.
#' @export
gamma_sampler_eta <- nimbleFunction(
  name = 'gamma_sampler_eta',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    X_phi <- control$designMatrix
    S <- dim(X_phi)[1]; M <- dim(X_phi)[2]#; Trep <- dim(X_p)[3]
    fixedIndexes <- 1
    ncov <-  if(!is.null(control$ncov)) control$ncov else 7
    d_bar <- ifelse(ncov <= 2, 1, 2)
    d_bar <- ifelse(d_bar <= ncov, d_bar, 2)
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:ncov
    ## `ynodes` must be aligned element-for-element with `znodes`. Model
    ## declaration order need not match `expandNodeNames()` order, so the
    ## fitting functions supply the nodes explicitly; falling back on the
    ## dependencies of `target` preserves the original behaviour.
    ynodes <- if (!is.null(control$ynodes)) {
      model$expandNodeNames(control$ynodes, returnScalarComponents = TRUE)
    } else {
      model$getDependencies(target, stochOnly = TRUE, self = FALSE)
    }
    PG_nodes <- model$expandNodeNames(control$PG, returnScalarComponents = TRUE)
    fixedeffects <- control$fixedEffects
    znodes <- model$expandNodeNames(control$znodes, returnScalarComponents = TRUE)
  },
  run = function() {
    
    # Get current values 
    b <- model$getParam(fixedeffects, 'mean') #model$beta.comm
    B <- model$getParam(fixedeffects , 'cov') #model$var.comm
    
    
    y <- values(model, ynodes) 
    z <- values(model, znodes)
    
    k <- y[z == 0] - 0.5
    N <- length(k)
    
    # Select covariates at samples site and z = 0
    
    X <- X_phi[z ==0, ]
    
    sel   <- which(z == 0)
    Omega <- values(model, PG_nodes)[sel]
    
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


# sampler for detection probability (p)

#' Polya-Gamma update for colonisation
#'
#' Draws the Polya-Gamma latent variables of the colonisation component for a single season, at the sites that were unoccupied in that season.
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
#'   \item{znodes}{Name of the latent occupancy node.}
#'   }
#'
#' @return A \pkg{nimble} sampler object.
#' @export
PG_sampler_eta <- nimbleFunction(
  name = 'PG_sampler_eta',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    X_p <- control$designMatrix
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:m
    gamma_nodes <- model$expandNodeNames(control$Gamma, returnScalarComponents = TRUE) 
    beta_nodes <- model$expandNodeNames(control$fixedEffects, returnScalarComponents = TRUE) 
    znodes <- model$expandNodeNames(control$znodes, returnScalarComponents = TRUE)
  },
  run = function() {
    
    # Get current values
    z <- values(model, znodes)
    length(z)
    
    xid <- which(z ==0) 
    X = X_p[xid, ] 
    nsize <- dim(X)[1]
    m <- dim(X)[2]
    
    n <- rep(1, nsize)  #
    
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
    
    
    Omega_new <- sample_Omega(X_gamma[1:nsize, 1:l], beta_gamma[1:l], n[1:nsize]) # this matches the result obtain from Alex's function
    
    for (i in 1:nsize) {
      model[[target]][xid[i]] <<- Omega_new[i]  # need to test here if betas that are not updated remain the same or what should the new value of beta be
    }
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
  },
  methods = list(
    reset = function() { }
  )
)



#' Conjugate coefficient update for colonisation
#'
#' Draws the active colonisation coefficients from their Gaussian full conditional distribution given the Polya-Gamma latent variables.
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
#'   \item{znodes}{Name of the latent occupancy node.}
#'   \item{ynodes}{Name of the occupancy node one season ahead. Optional; when omitted the stochastic dependencies of \code{target} are used instead.}
#'   }
#'
#' @return A \pkg{nimble} sampler object.
#' @export
beta_sampler_eta <- nimbleFunction(
  name = 'beta_sampler_eta',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    X_p <- control$designMatrix
    S <- dim(X_p)[1]; M <- dim(X_p)[2]#; Trep <- dim(X_p)[3]
    indexes_covariates <- if(!is.null(control$indexes_covariates)) control$indexes_covariates else 1:M
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    ## `ynodes` must be aligned element-for-element with `znodes`. Model
    ## declaration order need not match `expandNodeNames()` order, so the
    ## fitting functions supply the nodes explicitly; falling back on the
    ## dependencies of `target` preserves the original behaviour.
    ynodes <- if (!is.null(control$ynodes)) {
      model$expandNodeNames(control$ynodes, returnScalarComponents = TRUE)
    } else {
      model$getDependencies(target, stochOnly = TRUE, self = FALSE)
    }
    znodes <- model$expandNodeNames(control$znodes, returnScalarComponents = TRUE)
    
    calcNodes <- model$getDependencies(target)
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    gamma_nodes <- model$expandNodeNames(control$Gamma, returnScalarComponents = TRUE) 
    PG_nodes <- model$expandNodeNames(control$PG, returnScalarComponents = TRUE)
  },
  run = function() {
    
    
    # Get current values
    b <- model$getParam(target, 'mean') #model$beta.comm
    B <- model$getParam(target, 'cov') #model$var.comm
    y <- values(model, ynodes)
    z <- values(model, znodes)
    
    
    k <- y[z == 0] - 0.5
    N <- length(k)
    
    # Select covariates at samples site and z = 0
    #z[zindex]==1 # acroos all replicated 
    #dim(X_p)
    X <- X_p[z==0, ]
    #dim(X)
    
    gamma <- values(model, gamma_nodes) #model$gamma
    
    sel   <- which(z == 0)
    Omega <- values(model, PG_nodes)[sel]
    
    
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
    
    for (i in 1:l) {
      model[[target]][index_present[i]] <<- beta_new[i]  # need to test here if betas that are not updated remain the same or what should the new value of beta be
    }
    
    update_other_nodes <- model$calculateDiff(calcNodesDepStage)
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
    
  },
  methods = list(
    reset = function() { }
  )
)
