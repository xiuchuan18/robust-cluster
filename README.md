# robust-cluster

MixSIT: model-based mixture clustering with independent split-t latent coordinates
and a hierarchy of structural constraints. This repository contains the estimator,
formal experiment drivers, final numerical results, and six vector PDF figures.
The manuscript source project is distributed separately.

## Structural variants

| Published name | Internal method code | Structure |
|---|---|---|
| MixSIT-VV | `MixSIT_None` | Component-specific triangular and orthogonal factors |
| MixSIT-EV | `MixSIT_U` | Shared triangular factor, component-specific orthogonal factors |
| MixSIT-EE | `MixSIT_Full` | Shared triangular and orthogonal factors |

Tail indices and split parameters remain component-specific. The default log-nu
heterogeneity penalty has strength 5; the formal sensitivity experiment varies it.

## Repository layout

```text
Stage1/              estimator, flip moves, C++ backend, build script
Stage2/              formal experiments and shared configuration
  results/           aggregate RDS files and CSV summaries
    calibrated-flexibility/   oracle profile, Banknote, K-selection phase map
  figures/           six final one-page vector PDFs
scripts/             environment setup and workflow runners
tests/               numerical, scoring, RNG, and artifact checks
docs/                protocol, experiment map, data, results, and figure guide
vendor/EMMIXskew/     pinned dependency source and upstream license
software-versions.csv
MANIFEST.csv         SHA-256 checksums of publication files (excluding itself)
LICENSE              GPL-3
```

Historical `Section4`, `Section5`, and `Section6` filenames are stable execution
identifiers. Following the manuscript reorganization, they correspond to its
controlled-evidence, empirical-benchmark, and recovery/selection/scalability
sections. They do not imply the current section numbering.

## Environment and reproduction

The recorded environment is Windows with **R 4.4.1**. Direct package versions
are pinned in `software-versions.csv`; `docs/sessionInfo.txt` records the environment
used for the supplied results. Use the matching R compilation toolchain
(Rtools44 on Windows). The setup script downloads pinned dependencies and builds
the local `.r-lib/`; a working Internet connection and compilers are required.

From the repository root, with Rscript on PATH:

```powershell
Rscript scripts/setup_runtime.R
Rscript tests/preflight.R
Rscript tests/method_revision.R
Rscript tests/nu_actual_baseline.R
Rscript tests/information_criteria.R
Rscript tests/wbar_axis_matching.R
Rscript tests/seed_audit.R
Rscript tests/calflex_seed_audit.R
Rscript tests/result_integrity.R
Rscript tests/calflex_integrity.R
Rscript tests/figure_integrity.R
```

To rebuild numerical summaries from the supplied aggregate results, without
fitting the experiments again:

```powershell
Push-Location Stage2
Rscript Final_Summaries.R
Pop-Location
```

The complete experiment workflow is computationally expensive:

```powershell
Rscript scripts/run_all.R
```

A fresh clone contains no per-job caches: the full runner fits the experiments
again, then rebuilds summaries and figures. It writes into `Stage2/results/` and
`Stage2/figures/`. Use a separate clone if you want to preserve the supplied
numerical record while rerunning. Completed jobs are resumable only when their
run identifier and RNG provenance match. Cache files, logs, installed packages,
and compiled binaries are excluded from Git.

## Interpretation and documentation

- [Protocol](docs/PROTOCOL.md): initialization, objective safeguards, scoring,
  failure handling, RNG provenance, and timing limits.
- [Experiments](docs/EXPERIMENTS.md): formal experiment designs and entry points.
- [Data](docs/DATA.md): benchmark sources and preprocessing.
- [Results](docs/RESULTS.md): aggregate files and summary denominators.
- [Figures](docs/FIGURES.md): figure-to-script mapping.

Oracle MAP ARI is a reference attained by the generating-model classifier; it is
not an upper bound on ARI. Returned partitions and formal convergence are reported
separately. Primary AIC/BIC/ICL use raw observed likelihood and a common soft
posterior entropy definition. Timings are empirical measurements from the recorded
machine and fitting protocol, not machine-independent complexity bounds.

The source is GPL-3. Vendored EMMIXskew retains its upstream copyright and license;
its compatibility changes are described in `vendor/EMMIXskew/PATCHES.md`.

## Four additional competitors

[Experiment guide](Stage2/CompetitorExtension/README.md) and [reported results](Stage2/results/competitor-extension). Rebuild Figures 1, 3, 5 and S1 from the supplied CSVs with Rscript scripts/render_competitor_figures.R.
