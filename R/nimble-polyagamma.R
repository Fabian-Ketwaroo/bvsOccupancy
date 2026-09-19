## Exact Polya-Gamma sampler of Polson, Scott and Windle (2013).


#' Coefficient of the alternating series used in the Polya-Gamma sampler
#'
#' Internal building block of the exact Polya-Gamma sampler of Polson, Scott
#' and Windle (2013); see \code{\link{samplepg}}.
#'
#' @param n,x,t,mu,z Arguments of the underlying algorithm.
#'
#' @return A scalar draw.
#' @keywords internal
#' @export
aterm <- nimbleFunction(
  run = function(n = integer(), x= double(0), t = double(0)){
    returnType(double(0))
    f=0
    if(x <= t) {
      f = log(pi) + log(n + 0.5) + 1.5*(-log(pi/2) - log(x)) - 2*(n + 0.5)*(n + 0.5)/x
    }
    else {
      f = log(pi) + log(n + 0.5) - x * pi^2/2 * (n + 0.5)*(n + 0.5)
    }    
    return(exp(f))
    
  }
)


#' Draw from an exponential distribution
#'
#' Internal building block of the exact Polya-Gamma sampler of Polson, Scott
#' and Windle (2013); see \code{\link{samplepg}}.
#'
#' @param n,x,t,mu,z Arguments of the underlying algorithm.
#'
#' @return A scalar draw.
#' @keywords internal
#' @export
exprnd <- nimbleFunction(
  run = function(mu = double(0)){
    returnType(double(0))
    return(-mu*log(1.0 - runif(1,0.0,1.0) ))
  }
)



#' Draw from the truncated gamma proposal of Windle (2013)
#'
#' Internal building block of the exact Polya-Gamma sampler of Polson, Scott
#' and Windle (2013); see \code{\link{samplepg}}.
#'
#' @param n,x,t,mu,z Arguments of the underlying algorithm.
#'
#' @return A scalar draw.
#' @keywords internal
#' @export
truncgamma <- nimbleFunction(
  run=function( ){
    returnType(double(0))
    
    c <- pi/2
    done <- FALSE
    while(!done){
      X = exprnd(1.0) * 2.0 + c
      gX = sqrt(pi/2) / sqrt(X)
      
      if(runif(1,0,1) <= gX ) done = TRUE
    }
    
    return(X)
  }
)




#' Draw from an inverse Gaussian distribution
#'
#' Internal building block of the exact Polya-Gamma sampler of Polson, Scott
#' and Windle (2013); see \code{\link{samplepg}}.
#'
#' @param n,x,t,mu,z Arguments of the underlying algorithm.
#'
#' @return A scalar draw.
#' @keywords internal
#' @export
randinvg <- nimbleFunction(
  run=function(mu = double(0)){
    returnType(double(0))
    u = rnorm(1,0,1)
    V = u*u
    out = mu + 0.5*mu *( mu*V - sqrt(4.0*mu*V + mu*mu * V*V) )
    
    if(runif(1,0,1) > mu /(mu+out) )  out = mu*mu / out
    
    return(out)
  }
)


#' Draw from a truncated inverse Gaussian distribution
#'
#' Internal building block of the exact Polya-Gamma sampler of Polson, Scott
#' and Windle (2013); see \code{\link{samplepg}}.
#'
#' @param n,x,t,mu,z Arguments of the underlying algorithm.
#'
#' @return A scalar draw.
#' @keywords internal
#' @export
tinvgauss <- nimbleFunction(
  run= function(z= double(0), t = double(0)) {
    returnType(double(0))
    
    mu <- 1.0/z
    
    done <- FALSE
    
    # Pick sampler
    if(mu > t) {
      # Sampler based on truncated gamma 
      # Algorithm 3 in the Windle (2013) PhD thesis, page 128
      while(!done) {
        u = runif(1,0.0, 1.0)
        X = 1.0 / truncgamma()
        
        if ( log(u) < (-z*z*0.5*X) ) {
          #nimbreak()
          done = TRUE
          #break
        }
      }
    }  
    else {
      # Rejection sampler
      X = t + 1.0
      while(X >= t) {
        X = randinvg(mu)
      }
    } 
    
    return(X)
    
  }
)



#' Exact draw from a Polya-Gamma PG(1, z) distribution
#'
#' Implements the accept-reject sampler of Polson, Scott and Windle (2013),
#' following the algorithm on page 130 of Windle's thesis.
#'
#' @param z Tilting parameter, typically the current linear predictor.
#'
#' @return A single Polya-Gamma draw.
#' @references
#' Polson, N. G., Scott, J. G. and Windle, J. (2013) Bayesian inference for
#' logistic models using Polya-Gamma latent variables.
#' \emph{Journal of the American Statistical Association} \strong{108}, 1339--1349.
#' @keywords internal
#' @export
samplepg <- nimbleFunction(
  run= function( z= double(0) ){
    returnType(double(0))
    
    # PG(b, z) = 0.25 * J*(b, z/2)
    z = abs(z)*0.5
    # Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
    MATH_2_PI <- 0.636619772367581343075535053490057448137838582961825794990
    t = MATH_2_PI
    
    done <- breakWhile <- FALSE
    
    # Compute p, q and the ratio q / (q + p)
    # (derived from scratch; derivation is not in the original paper)
    K = z*z/2.0 + pi^2/8.0
    logA = log(4.0) - log(pi) - z
    logK = log(K)
    Kt = K * t
    w = sqrt(pi/2)
    
    logf1 = logA + pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt
    logf2 = logA + 2*z + pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt
    p_over_q = exp(logf1) + exp(logf2)
    ratio = 1.0 / (1.0 + p_over_q)
    
    #breakWhile <- FALSE
    
    # Main sampling loop; page 130 of the Windle PhD thesis
    while(!breakWhile) {
      # Step 1: Sample X ? g(x|z)
      u = runif(1,0.0,1.0)
      if(u < ratio) {
        # truncated exponential
        X = t + exprnd(1.0)/K
      }
      else{
        # truncated Inverse Gaussian
        X = tinvgauss(z, t)
      }
      
      # Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
      i = 1
      Sn = aterm(0, X, t)
      U = runif(1,0.0,1.0) * Sn
      asgn = -1
      even = FALSE
      
      while(!done) {
        
        Sn = Sn + asgn * aterm(i, X, t)
        
        # Accept if n is odd
        if(!even & (U <= Sn)) {
          X = X * 0.25
          return(X)
        }
        
        # Return to step 1 if n is even
        if(even & (U > Sn) ) {
          # nimbreak()
          breakWhile = done = TRUE 
          # break # this is where the problem is 
        }
        
        even = !even
        asgn = -asgn
        i = i+1
        
      }
      
    }
    return(X) 
  }
)




#' Draw from a Polya-Gamma PG(n, z) distribution
#'
#' Sums \code{n} independent PG(1, z) draws.
#'
#' @param n Integer shape parameter, the number of Bernoulli trials.
#' @param z Tilting parameter.
#'
#' @return A single Polya-Gamma draw.
#' @keywords internal
#' @export
pgr <- nimbleFunction(
  run = function( n = integer(), z = double(0) ){
    returnType(double(0))
    
    x = 0
    for (i in 1:n) {
      x =  x + samplepg(z)
    }
    
    return(x)
    
  }
)





#' Vectorised Polya-Gamma update
#'
#' Draws one Polya-Gamma latent variable per row of the design matrix, tilted
#' by that row's linear predictor.
#'
#' @param X Design matrix restricted to the active columns.
#' @param beta Coefficients of the active columns.
#' @param n Vector of trial counts, one per row of \code{X}.
#'
#' @return A vector of Polya-Gamma draws, one per row of \code{X}.
#' @keywords internal
#' @export
sample_Omega = nimbleFunction(
  run = function( X = double(2), beta = double(1), n = double(1) ){
    returnType(double(1))
    
    nsize = length(n) # n is the number of samples for each observation
    Omega_vec = numeric(nsize)
    m <- dim(X)[2]
    
    for (i in 1:nsize) {
      b = X[i, 1:m ] %*%  beta[1:m] 
      Omega_vec[i] = pgr(n[i], b[1,1]) 
    }
    
    return(Omega_vec)
    
  }
)
