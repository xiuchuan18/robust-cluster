# Formal experiments

Run drivers from `Stage2/`, or use `scripts/run_all.R` from the repository root.
The `Section*` names are historical file identifiers; consult the scientific
description below when mapping them to the current manuscript.

| Experiment | Design | Driver |
|---|---|---|
| Controlled simulation | Six configured scenarios, 100 repetitions each; nine methods | `Section4_Run.R` |
| Generating-model MAP reference | Same six cells and data seeds | `Section4_BayesCeiling.R`, `CalFlex_BayesProfile.R` |
| Latent-scale contamination diagnostic | 100 repetitions of S7_CellWise10; three MixSIT variants; exact observed-direction axis matching | `Section4_Run.R`, `Section4_WbarAnalysis.R` |
| Clean empirical benchmarks | Biopsy, AIS, Wine; nine methods | `Section5_RealData.R` |
| Semi-real contamination | 50 copies per dataset; cell probability .05, additive Gaussian SD 5 after standardization | `Section5_CellWise05.R` |
| Parameter recovery | D3, m=10, K=3, separation=6; N=200,400,800,1600,3200; 50 repetitions; VV/EV/EE | `Section6_Recovery.R` |
| Penalty sensitivity | Same D3 family; N=200,400,800,1600; lambda=0,.5,1,2,5,10,25; 10 repetitions; EV | `Section6_LambdaSensitivity.R` |
| Main K selection | Eight configured generating cells; K=1,...,6; five methods; 30 repetitions | `Section6_KSelect.R` |
| Scalability | Sweep N=200,400,800,1600; m=3,5,10,15,20; K=2,3,4,5,6; anchor (800,10,3); 30 repetitions; five methods | `Section6_Scale.R` |
| Banknote | Clean and 50 contaminated copies; same nine methods and contamination rule | `CalFlex_Banknote.R` |
| K-selection phase map | Eight prespecified conditions; K=1,...,5; five methods; 50 repetitions (10,000 candidate fits) | `CalFlex_PhaseMap.R` |

Configuration and truth definitions live in `00_CONFIG.R`, `Section4_Lib.R`,
`Section6_TruthLib.R`, `Section6_KSelectLib.R`, `Section6_ScaleLib.R`, and
`CalFlex_FormalLib.R`. Semantic RNG streams preserve each job's data and fitting
draws across different scheduling orders. Run identifiers and RNG namespaces are
preserved in source and RDS provenance.

The production library and sensitivity study use the same actual-incoming-nu
baseline safeguard. The published sensitivity study changes only the strength
of the log-nu heterogeneity penalty. The release includes the active D3 truth
and log penalty; unused exploratory alternatives are excluded.

`Final_Summaries.R` reads aggregate RDS files and writes five numerical CSV
summaries. It never reads per-job caches. `Section6_ScalingFigure.R` reproduces
the formal median-time plot without writing manuscript fragments or running
unreported selection grids.

The standalone calibrated-flexibility runner assumes the main simulation and
three-dataset aggregate results already exist, as they do in this release.
For complete regeneration use `scripts/run_all.R`, which orders those dependencies.
