## Conjugate samplers for the community-level parameters of the multi-species model.


#' Community mean update for a single coefficient
#'
#' Conjugate normal update of one community-level mean in a multi-species model. Only species for which the covariate is currently included contribute; if no species includes it, the mean is drawn from its prior. Use \code{\link{Comm_multivariate_mean_sampler}} when categorical covariates are present.
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
#'   \item{beta.species}{Species-level coefficient nodes for this covariate.}
#'   \item{gamma.species}{Species-level inclusion indicators for this covariate.}
#'   \item{sig2.comm}{Community variance node for this covariate.}
#'   }
#'
#' @return A \pkg{nimble} sampler object.
#' @export
Comm_mean_sampler <- nimbleFunction(
  name = 'Comm_mean_sampler',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    # Get hyperprior mean and standard deviation from the target node
    m0 <- model$getParam(target, 'mean')
    v0 <- model$getParam(target, 'var')
    prec0 <- 1 / v0                  # Prior precision
    
    # Expand and store node names
    beta.species  <- model$expandNodeNames(control$beta.species, returnScalarComponents = TRUE)
    gamma.species <- model$expandNodeNames(control$gamma.species, returnScalarComponents = TRUE)
    sig2.comm     <- model$expandNodeNames(control$sig2.comm, returnScalarComponents = TRUE)
    
    S <- length(gamma.species)
  },
  run = function() {
    
    gamma <- values(model, gamma.species)
    beta  <- values(model, beta.species)
    sig2  <- values(model, sig2.comm)[1] # Current community variance
    
    N <- 0
    sum_beta <- 0
    
    # Fast C++ scalar loop: count active species and sum active beta values
    for (s in 1:S) {
      if (gamma[s] == 1) {
        N <- N + 1
        sum_beta <- sum_beta + beta[s]
      }
    }
    
    # Conjugate Normal draw if active species exist; otherwise sample from prior
    if (N > 0) {
      prec_data <- N / sig2
      prec_post <- prec0 + prec_data
      var_post  <- 1 / prec_post
      
      mean_post <- var_post * ((m0 * prec0) + (sum_beta / sig2))
      
      model[[target]] <<- rnorm(n = 1, mean = mean_post, sd = sqrt(var_post))
    } else {
      model[[target]] <<- rnorm(n = 1, mean = m0, sd = sqrt(v0))
    }
    
    # Recalculate target log probability and update MCMC chain trace
    calculate(model, target)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
  },
  methods = list(
    reset = function() { }
  )
)




#' Community variance update
#'
#' Conjugate inverse-gamma update of the community-level variance of one covariate group, using the exchangeable block structure that keeps categorical effects invariant to the reference level. Only species that currently include the covariate contribute.
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
#'   \item{beta.species}{Species-level coefficient nodes for the columns in this group.}
#'   \item{gamma.species}{Species-level inclusion indicators for this group.}
#'   \item{beta.comm}{Community mean nodes for the columns in this group.}
#'   }
#'
#' @return A \pkg{nimble} sampler object.
#' @export
Comm_var_sampler <- nimbleFunction(
  name = 'Comm_var_sampler',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    a <- model$getParam(target, 'shape')
    b <- model$getParam(target, 'scale')
    
    beta.comm     <- model$expandNodeNames(control$beta.comm,     returnScalarComponents = TRUE)
    gamma.species <- model$expandNodeNames(control$gamma.species, returnScalarComponents = TRUE)
    
    L <- length(beta.comm)
    S <- length(gamma.species)
    
    # flatten to a single vector, ordered dummy-within-species:
    #   entry (s-1)*L + l  is species s, dummy l
    beta.flat <- character(L * S)
    for (l in 1:L) {
      nodes_l <- model$expandNodeNames(control$beta.species[l], returnScalarComponents = TRUE)
      for (s in 1:S) beta.flat[(s - 1) * L + l] <- nodes_l[s]
    }
    
    Cinv <- 2 * (diag(L) - matrix(1, L, L) / (L + 1))
  },
  run = function() {
    
    gamma <- values(model, gamma.species)
    mu    <- values(model, beta.comm)
    bet   <- values(model, beta.flat)        # length L*S, resolved at setup
    
    Nact <- 0
    quad <- 0
    
    for (s in 1:S) {
      if (gamma[s] == 1) {
        Nact <- Nact + 1
        d <- nimNumeric(L)
        for (l in 1:L) {
          d[l] <- bet[(s - 1) * L + l] - mu[l]
        }
        quad <- quad + inprod(d[1:L], (Cinv[1:L, 1:L] %*% d[1:L])[1:L, 1])
      }
    }
    
    if (Nact > 0) {
      model[[target]] <<- rinvgamma(n = 1, shape = a + L * Nact / 2, scale = b + quad / 2)
    } else {
      model[[target]] <<- rinvgamma(n = 1, shape = a, scale = b)
    }
    
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
  },
  methods = list(reset = function() { })
)


# Community mean sampler with one variance per group not dummy level. This ensure invarance in the species-level effects

#' Joint community mean update
#'
#' Conjugate multivariate normal update of the full vector of community-level means. Each covariate group contributes through its own variance and its exchangeable block correlation, so categorical covariates are handled correctly. Species that currently exclude a group contribute nothing to that block.
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
#'   \item{beta.species}{Matrix of species-level coefficient nodes.}
#'   \item{gamma.species}{Matrix of species-level inclusion indicators.}
#'   \item{sig2.comm}{Vector of community variance nodes, one per covariate group.}
#'   \item{indexes_covariates}{Integer vector mapping columns of \code{designMatrix} to covariate groups.}
#'   }
#'
#' @return A \pkg{nimble} sampler object.
#' @export
Comm_multivariate_mean_sampler <- nimbleFunction(
  name = 'Comm_multivariate_mean_sampler',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    
    calcNodes <- model$getDependencies(target)
    
    mu0  <- model$getParam(target, 'mean')
    cov0 <- model$getParam(target, 'cov')
    D_beta <- length(mu0)
    
    P0     <- inverse(cov0)
    P0_mu0 <- (P0 %*% mu0)[1:D_beta, 1]
    
    beta.nodes  <- model$expandNodeNames(control$beta.species,  returnScalarComponents = TRUE)
    gamma.nodes <- model$expandNodeNames(control$gamma.species, returnScalarComponents = TRUE)
    sig2.nodes  <- model$expandNodeNames(control$sig2.comm)
    
    indexes_covariates <- as.numeric(control$indexes_covariates)
    S <- length(beta.nodes) / D_beta
    
    n_groups <- max(indexes_covariates)
    
    # for each group: its parameter columns, padded, plus size and C^{-1}
    grp_cols <- matrix(0, nrow = n_groups, ncol = D_beta)
    grp_size <- numeric(n_groups)
    for (g in 1:n_groups) {
      cols <- which(indexes_covariates == g)
      grp_size[g] <- length(cols)
      grp_cols[g, 1:length(cols)] <- cols
    }
    Lmax <- max(grp_size)
    
    # C^{-1} for each group, padded to Lmax x Lmax
    Cinv_all <- array(0, dim = c(n_groups, Lmax, Lmax))
    for (g in 1:n_groups) {
      L <- grp_size[g]
      Cinv_all[g, 1:L, 1:L] <- 2 * (diag(L) - matrix(1, L, L) / (L + 1))
    }
  },
  
  run = function() {
    
    gamma_val <- values(model, gamma.nodes)
    beta_val  <- values(model, beta.nodes)
    sig2_val  <- values(model, sig2.nodes)     # length n_groups
    
    v_data <- nimNumeric(D_beta, value = 0)
    P_post <- P0
    
    for (g in 1:n_groups) {
      L  <<- grp_size[g]
      s2 <- sig2_val[g]
      if (s2 > 0) {
        
        N_g   <- 0
        sumb  <- nimNumeric(L, value = 0)
        
        for (s in 1:S) {
          if (gamma_val[(g - 1) * S + s] == 1) {
            N_g <- N_g + 1
            for (l in 1:L) {
              sumb[l] <- sumb[l] + beta_val[(grp_cols[g, l] - 1) * S + s]
            }
          }
        }
        
        if (N_g > 0) {
          # precision block: N_g * C^{-1} / s2 ; information: C^{-1} %*% sumb / s2
          contrib <- (Cinv_all[g, 1:L, 1:L] %*% sumb[1:L])[1:L, 1] / s2
          for (l in 1:L) {
            v_data[grp_cols[g, l]] <- v_data[grp_cols[g, l]] + contrib[l]
            for (m in 1:L) {
              P_post[grp_cols[g, l], grp_cols[g, m]] <-
                P_post[grp_cols[g, l], grp_cols[g, m]] +
                N_g * Cinv_all[g, l, m] / s2
            }
          }
        }
      }
    }
    
    M_post <- solve(P_post, P0_mu0 + v_data)[1:D_beta]
    
    model[[target]] <<- rmnorm_chol(n = 1, mean = M_post,
                                    cholesky = chol(P_post), prec_param = TRUE)
    
    calculate(model, calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  
  methods = list(
    reset = function() { }
  )
)