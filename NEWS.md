# bvsOccupancy 0.3.0

## Bug fixes

Two corrections to the samplers, both supplied by the author from upstream work
on the original sampler code. Both change results.

* The four Polya-Gamma samplers (`PG_sampler_psi`, `PG_sampler_p`,
  `PG_sampler_phi`, `PG_sampler_eta`) built the index of included columns with
  `length(gamma)`, the number of covariate GROUPS, where it should have been
  `length(indexes_covariates)`, the number of design COLUMNS. Whenever a
  categorical covariate was in the model the design matrix had more columns than
  groups, and the trailing columns were silently left out of the linear
  predictor used to draw Omega. Since the coefficient and indicator updates both
  condition on Omega, every part of the fit was affected. Models with only
  continuous covariates were unaffected, because there columns and groups
  coincide.

* `samplepg()` left the rejection branch of its alternating-series sampler
  through `breakWhile = done = TRUE`, which exited both loops and returned the
  current proposal instead of rejecting it. It now sets only `done`, leaves the
  inner loop and draws a fresh proposal, as the algorithm requires. The inner
  flag is also reset for each new proposal. Omega draws were previously not from
  the Polya-Gamma distribution.

## New features

* The target model size behind the variable-selection moves is now an argument.
  `bvsSSOM()` and `bvsMSOM()` take `d.bar.occ` and `d.bar.det`; `bvsSDOM()` and
  `bvsMSDOM()` additionally take `d.bar.per` and `d.bar.col`. Each sets the
  number of active covariates the prior expects for that component, so the
  prior inclusion probability per covariate is `d.bar / ncov`. It was
  previously fixed at 2 and could not be changed without editing the samplers,
  which made a sensitivity analysis on the prior impossible. Leaving the
  arguments unset reproduces the old behaviour exactly, and the resolved values
  are stored on the fitted object as `fit$d.bar`.

* The four gamma samplers accept `d.bar` in their `control` list, for anyone
  assembling a model by hand rather than through the model functions.

## Bug fixes (continued)

* `.buildAndRun()` assembled the fitted object by naming each field of `extra`
  explicitly, so a field added there and not named was dropped without warning.
  Fields it does not name are now carried through, and `fit$d.bar` is populated.


* With exactly one candidate covariate the old rule set the target model size
  equal to the number of candidates, which makes the implied Beta(1, 0) prior
  degenerate: the add move always accepted and the delete move never did, so
  that covariate was forced into every model with inclusion probability 1. The
  default is now 0.5 in that case, giving a prior inclusion probability of 0.5.
  Designs with two or more candidate covariates are unaffected.

## Documentation

* The vignette gains a section on what the selector does and does not promise,
  with operating characteristics at a representative design: active covariates
  recovered in every replicate, null covariates with a median inclusion
  probability of 0.08, and about a 2.5% chance per data set that one of five
  nulls crosses 0.5.

* The vignette gains a section on the prior over model size. The add, delete
  and swap moves implement a Beta-Binomial prior whose inclusion probability is
  `d_bar / ncov` with `d_bar` fixed at 2, so the prior expects about two active
  covariates however many are offered and the prior inclusion probability per
  covariate falls as candidates are added: 0.5 with four candidates, 0.2 with
  ten, 0.1 with twenty. With few candidates the prior is close to uniform and a
  null covariate cannot be driven near zero, which is the prior behaving as
  designed rather than the selector failing. Inclusion probabilities should be
  read against that prior.

Anyone who ran a model containing a factor, or who is relying on earlier Omega
draws, should re-run. Results from models with only continuous covariates are
affected by the `samplepg()` fix alone.

# bvsOccupancy 0.2.0

## New features

* `bvsMSDOM()` fits a multi-species dynamic occupancy model: several species at
  the same sites over several seasons, with variable selection on initial
  occupancy, persistence, colonisation and detection, separately for every
  species, and community-level means and variances tying the species together.
  It is the combination of `bvsSDOM()` and `bvsMSOM()` and uses the same
  samplers; nothing in the method is new.

* `simMSDOM()` simulates from that model, with optional unequal replication
  across site-seasons.

# bvsOccupancy 0.1.0

First release, accompanying the manuscript.

## Models

* `bvsSSOM()` fits a single-season single-species occupancy model with
  variable selection on occupancy and detection.
* `bvsSDOM()` fits a dynamic occupancy model with variable selection on
  initial occupancy, persistence, colonisation and detection.
* `bvsMSOM()` fits a multi-species occupancy model with species-specific
  variable selection and conjugate community-level parameters.

## Supporting functions

* `bvsDesign()` and `bvsDesignArray()` build design matrices and derive the
  covariate grouping from a model formula, so that categorical covariates are
  selected as a block without the grouping having to be written out by hand.
* `bvsPrior()` constructs the prior on the regression coefficients.
* `inclusionProbs()` and `medianModel()` summarise the selection results.
* `simSSOM()`, `simSDOM()` and `simMSOM()` simulate data from each model.

## Changes to the reference implementation

Four changes were made while packaging the manuscript code. The first two
are corrections to the samplers themselves and affect results; the last two
are defensive checks in the R layer.

* The persistence and colonisation samplers now take the node holding the
  next season's occupancy state from their `control` list. Previously it came
  from `model$getDependencies()`, which returns nodes in model declaration
  order, while the current season's state came from `expandNodeNames()`. For a
  matrix-valued state declared inside nested loops these two orderings differ,
  so the two vectors could be compared out of step. The fitting functions now
  pass both node vectors explicitly and check that they align.

* Design matrices are flattened in the order of the node names the model
  actually reports, rather than assuming a column-major layout, and the row
  counts are checked against the number of nodes. In the reference code the
  flattened persistence design covered all seasons while the state vector
  covered all but the last, so the two had different lengths and R's recycling
  rules silently selected the wrong rows.

* Priors supplied through `prior.occ`, `prior.per`, `prior.col` and
  `prior.det` are now validated against the design they belong to. The length
  of the prior mean, the dimensions and symmetry of the prior covariance, and
  its positive definiteness are all checked before the model is built, so a
  mismatched prior produces a clear message in R rather than an opaque failure
  inside nimble.

* `simSDOM()` compared `dim(Jindex)` against a double vector, so the
  comparison was never true and every valid `Jindex` was rejected. Simulating
  unequal numbers of visits across sites and seasons now works.

* `indexes_covariates` is now stored as a double vector rather than an integer
  one. The nimbleFunctions that consume it (`dLgamma`, `compute_predictor` and
  `build_block_cov`) declare the argument as `double(1)`, so an integer vector
  compiled to `NimArr<1, int>` and the generated C++ failed to match the
  double signature. The symptom was a model that built and configured without
  complaint in R and then failed at `compileNimble()` with "no matching
  function for call to rcFun_...". The sampler setup code coerces the vector
  as well, so a hand-built `control` list cannot reintroduce the mismatch.

* `gamma_sampler_psi()` now stops with a clear message when
  `indexes_covariates` is absent from `control`. It previously fell through to
  `print()`, which returned a character vector and pushed the failure into the
  C++ stage.
