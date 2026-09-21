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

