# Final numerical record

Aggregate RDS files are the scientific record. CSV files are numerical exports
or summaries of that record. Internal codes `MixSIT_None`, `MixSIT_U`, and
`MixSIT_Full` mean VV, EV, and EE. Per-job caches are excluded.

| File stem under Stage2/results | Content |
|---|---|
| Section4_Simulation_Results | Six lists of 100 repetitions with method endpoints |
| Section4_BayesCeiling | Oracle generating-model MAP ARI reference, CSV/RDS |
| Section4_Wbar_S7_CellWise10 | 100 replicate diagnostic records with contamination masks, matched latent-scale weights, and axis maps |
| Section4_Wbar_Summary | Contamination diagnostic numerical summary |
| Section5_RealData | 27 dataset-method endpoints, CSV/RDS |
| Section5_CellWise05 | 1,350 contaminated dataset-method endpoints, CSV/RDS |
| Section6_Recovery | 22,500 parameter-axis rows; 750 fits |
| Section6_LambdaSensitivity | 8,400 parameter-axis rows; 280 fits |
| Section6_KSelect | 7,200 candidate-K fits |
| Section6_Scale | 1,800 timing records; intersecting sweeps share the anchor configuration |
| Section6_Scaling_Points | Median-time plotting data |
| Section6_Recovery_Summary | Per-variant/N recovery metrics |
| Section6_Recovery_Slopes | Upper-half sample-size log-log recovery slopes |
| Section6_Lambda_Summary | Paired-replicate lambda sensitivity metrics |
| Section6_KSelect_Summary | Cellwise completion, hit rate, signed error, and equal-cell averages for AIC/BIC/ICL |
| Section6_Scale_Summary | Mean-time empirical exponents, iteration counts, cap rates, anchor and m=20 times |

Under `calibrated-flexibility/`, `CalFlex_Banknote_Clean` and
`CalFlex_Banknote_CellWise05` hold endpoints, `CalFlex_Banknote_Summary` summarizes
them, `CalFlex_BayesCeiling` holds the same generating-MAP reference used for
the controlled cells, and `CalFlex_CeilingAttainment_ByCell`/`CalFlex_CeilingAttainment_Profile` contain
reference-attainment summaries. `CalFlex_PhaseMap` holds 10,000 candidate fits;
`CalFlex_PhaseMap_Paths` and `CalFlex_PhaseMap_Summary` describe complete-path
selection and incomplete paths. Primary MSGHD scores use the actual five-vector
parameter count. Alternative counting conventions are explicitly labeled
sensitivity fields, not primary criteria.

## Denominators and conventions

Clustering accuracy is conditional on a returned partition; failures and return
counts remain visible. A returned partition does not imply a convergence flag.
K selection requires every candidate on the prespecified K grid to return with a
finite criterion. Missing candidates are not removed to create a shorter grid.
Hit rates and signed errors are conditional on these complete paths. Main-table
averages give every represented scenario equal weight.

Recovery uses exact global component-overlap matching and sign-invariant
latent-axis angular matching. Fit-level quantities are counted once per fit,
component-level location errors once per component, and axis errors once per
axis. Reported recovery slopes use N at or above the median of the five sample
sizes. Lambda summaries use repetitions available at every lambda within N and
average per-fit axis RMSEs; this differs from pooling all squared errors.

Scaling exponents regress log mean total time on log swept size. The plotted
curves use median total time. `t_base/iters` covers recorded main-loop iterations,
while total time includes additional overhead. Cap fractions use each method's
cap (100 for MSGHD, 500 for MixSIT); GMM iteration fields are missing. See the
timing comparability limits in `PROTOCOL.md`.

The release retains all existing formal numerical and figure files byte-for-byte.
Five summary CSVs are extracted from the formal aggregates for this code release.
`MANIFEST.csv` covers publication files only and excludes itself. Running the
full workflow regenerates it after replacing outputs.

See the [four-competitor revision](../Stage2/CompetitorExtension/README.md) for initialization, scope, return counts, and the NPICA scoring and timing protocols.
