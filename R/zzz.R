## Registration of the placeholder Polya-Gamma distribution.
##
## Model code in this package contains `Omega[1:n] ~ dPG()`. nimble can
## auto-register a user distribution when it meets one in model code, but doing
## it once at load time is cheaper and keeps the message out of every call to
## nimbleModel(). Registration is wrapped in a `try` because re-registering an
## existing distribution is harmless but noisy.

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  try(
    suppressMessages(
      nimble::registerDistributions(
        list(
          dPG = list(
            BUGSdist = "dPG()",
            types    = c("value = double(1)"),
            pqAvail  = FALSE,
            discrete = FALSE
          )
        ),
        userEnv = ns,
        verbose = FALSE
      )
    ),
    silent = TRUE
  )
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion(pkgname)
  packageStartupMessage(
    "bvsOccupancy ", v,
    ": Bayesian variable selection for occupancy models.\n",
    "Models are compiled by nimble on first use, which can take a few minutes."
  )
}

## Names that appear inside nimbleCode() blocks. They are BUGS variables, not R
## objects, so R's code checker reports them as undefined. Declaring them here
## keeps `R CMD check` quiet without changing any behaviour.
utils::globalVariables(c(
  ## data and constants
  "nsites", "nvisits", "nyears", "nspecies",
  "X_psi", "X_phi", "X_eta", "X_p", "X_z",
  "ic_psi", "ic_phi", "ic_eta", "ic_p", "ic_z",
  ## parameters
  "beta_psi", "beta_phi", "beta_eta", "beta_p", "beta_z",
  "gamma_psi", "gamma_phi", "gamma_eta", "gamma_p", "gamma_z",
  "mu_beta_z", "mu_beta_p", "tau_beta_z", "tau_beta_p",
  "TBz", "TBp",
  ## BUGS language functions
  "logit<-"
))
