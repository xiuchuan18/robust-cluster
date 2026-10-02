# Experimental protocol

## Initialization budget

`R_RESTART = 1` means one full EM/ECM optimization. For every non-GMM method,
that fit begins from exactly one k-means call with `nstart = 5`. These are two
different layers of the budget: five Lloyd starts choose one partition, and
that one partition seeds one complete model fit.

GMM is the sole exception: `mclust::Mclust(..., initialization = NULL)` retains
the package's deterministic model-based hierarchical initialization.

## Package-specific implementation

- tMM, MSN, MST (`EMMIXskew` 1.0.3): one
  `stats::kmeans(..., nstart = 5)` partition is passed to the package's
  `initEmmix`, then to `EmSkew(init = ...)`. This preserves the native
  initializer but blocks `init.mix`'s hidden ten-random-start fallback.
- MSGHD (`MixGHD` 2.3.7): `gpar0` is constructed using the exact internal
  `rmgparMS` sequence: k-means centers and labels, `combinewk`, `rgparMS`, and
  ten `EMgrstepMS` updates. Only k-means `nstart` changes from 1 to 5.
- CGMM (`ContaminatedMixt` 1.3.8): `initialization = "kmeans"`; its native
  branch calls k-means once with `nstart = 5`.
- MixSIT: `init_methods = "kmeans"`, `n_short = 1`, and the estimator's k-means
  call uses `nstart = 5`.

The preflight test proves that the MSGHD wrapper at `nstart = 1` is identical
to the package-default `method = "km"` path in both map and log-likelihood
trajectory.

## Parallel and timing policy

The machine has 8 physical and 16 logical cores. Worker counts vary with experiment memory use (3--14); the full runner records
the count for each module. BLAS/OpenMP thread counts are fixed to one per job to avoid nested
oversubscription. The Section 6 scalability grid runs alone at the same fixed
eight-process load for all cells, so its reported wall times remain comparable.
Every driver uses the same checked process queue: a non-zero child exit, timeout,
or missing/invalid result stops the driver and preserves the child log.

## Cache and provenance policy

Per-job caches are written to a temporary file and renamed only after successful
serialization. A cache is reusable only when its embedded run identifier is
`mixsit-final-nstart5-eta5-2026-09-25` and its RNG scheme matches the source.
The tracked aggregate results are not used as completion markers. Package
versions are declared in `software-versions.csv`; `tests/preflight.R` enforces
them before fitting, and `docs/sessionInfo.txt` records the realized
environment.

Cache provenance includes model and information-criterion revision identifiers.
The publication repository contains aggregate outputs rather than per-job caches;
a clean clone recomputes the jobs when the full runner is invoked.

## Information criteria

All primary scores use the raw observed log-likelihood and are minimized:
`AIC = -2*ell + 2*p`, `BIC = -2*ell + p*log(N)`, and
`ICL = BIC + 2*H`, where `H` is soft posterior classification entropy.
Posterior rows are normalized, and probabilities are floored at `1e-10`
for the numerical entropy evaluation. One-component entropy is exactly zero;
missing or invalid posteriors produce a missing ICL, never zero entropy.
MSGHD uses
`p = (K-1) + K*(5*m + m*(m-1)/2)`; its five fitted vectors are mu, alpha,
phi, omega and lambda. Package and former over-count conventions are retained
only as labeled phase-map sensitivity scores.

MixGHD stores its posterior matrix in an S4 vector slot. The shared extractor
restores its column-major N-by-K dimensions before entropy is computed; passing
the vector directly to `as.matrix` would incorrectly treat it as one column.
The regression test compares the restored matrix with native `weightsMS`.

Both the production and installed-penalty nu updates compare candidates with
the actual incoming nu, including the appropriate penalty, and preserve the
incoming value exactly when a candidate worsens the conditional objective.
The full AECM observed-objective guard retains its documented numerical tolerance.

## Data scope

- Section 4: six scenarios, 100 replicates each.
- Section 5, primary arm: original Biopsy, full 11-variable AIS, and Wine data;
  no deduplication, jitter, feature selection, or artificial corruption.
- Section 5, semi-real arm: 50 independently perturbed copies of each primary
  dataset. After the clean variables have been standardized, each cell is
  independently selected with probability `alpha = 0.05` and receives additive
  `N(0, 5^2)` noise. The perturbed matrices are not standardized again. Class
  labels and the prespecified number of components are unchanged and used only
  for evaluation.
- Section 6: recovery, lambda sensitivity, main K-selection,
  and scalability. The angle-information numerical companion is not part of
  the production rerun; direction estimation is discussed in the manuscript
  using established identification theory and the delta method.

## Parameter-recovery matching

Recovery errors are evaluated after exact one-to-one assignments. Generating
components are matched to fitted components by maximizing total classification
overlap. Within each matched component, the observed-space latent directions
(the normalized columns of `t(U) %*% t(O)`) are matched by minimizing the total
sign-invariant angular loss. Both assignments use the deterministic exact
dynamic-programming solver in `Stage2/Section6_TruthLib.R`; no greedy matching
or external assignment package is used. The fitted component and axis selected
for every target are stored in the final recovery and lambda-sensitivity RDS
files together with the matching-rule identifiers.

## Section 5 random-number and failure policy

For semi-real replicate `r`, `mk_seed("Section5", dataset, "cellwise", 0.05,
5, "contam", r)` determines both the Bernoulli mask and additive Gaussian
draws. Each method then receives its own fitting stream through the existing
`fit_multistart` semantic tag `Section5CellWise05/<dataset>` and replicate
number. Consequently, all methods see the identical perturbed copy within a
replicate, while contamination, initialization, and different methods never
share an advancing random stream.

A numerical failure is recorded as a failed return and is never replaced by a
different initialization. ARI summaries are conditional on returned
partitions and must be reported with the return count. Formal package
convergence flags are reported separately because a returned classification
and a certified converged optimum are not the same event.

## Safeguards and timing scope

The recorded non-increasing-step flags and trajectories describe the main AECM
loop. Flip refinement is accepted by its final objective comparison, and its
internal iterations are not merged into those main-loop diagnostics. The observed
objective acceptance guard permits numerical changes within its 1e-10 tolerance;
the output is not a proof of exact arithmetic monotonicity of every substep.

Axis sign/reciprocal canonicalization preserves density. For EE, a desired sign
transformation outside the fixed Cayley chart can be skipped, so the stored
representative need not satisfy every preferred sign convention.

Scalability uses eight concurrent single-threaded fits. Supplied VV/EV timings
come from the final safeguarded rerun; unaffected GMM, MSGHD, and EE timings were
retained from their original runs under the same stated load. They are not a
simultaneous paired timing rerun. `t_base / iters` describes the main-loop time
per recorded iteration; `t_tot` also includes initialization and flip overhead.
