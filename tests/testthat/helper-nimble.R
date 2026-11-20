## Anything that touches a nimble model is skipped unless nimble is installed
## and can actually build one. Checking that a trivial model builds is stricter
## than checking that the package is on the library path, and it keeps the
## suite green on machines where nimble is present but not working.

has_working_nimble <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    ok <- FALSE
    if (requireNamespace("nimble", quietly = TRUE)) {
      ok <- tryCatch({
        code <- nimble::nimbleCode({ x ~ dnorm(0, 1) })
        m <- nimble::nimbleModel(code, inits = list(x = 0), calculate = FALSE)
        inherits(m, "modelBaseClass") || !is.null(m)
      }, error = function(e) FALSE)
    }
    cached <<- isTRUE(ok)
    cached
  }
})

skip_without_nimble <- function() {
  testthat::skip_if_not(has_working_nimble(), "nimble is not available")
}
