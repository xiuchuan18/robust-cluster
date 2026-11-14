# Additional competitor experiments
Run the commands below from Stage2/CompetitorExtension. The reported results are already included.
The four displayed names map to the frozen implementation IDs as follows: MMST = FW2014_MSt,
NPICA = ZH2019_ICA, MMSCN = T2025_Manly, PMSCN = T2025_Power.

| File in Stage2/results/competitor-extension/ | Inputs | Attempted fits |
|---|---:|---:|
| simulation.csv | 600, S1-S6 x 100 | 2400 |
| section5.csv | 4 clean + 200 contaminated | 816 |
| section6.csv | 250 recovery + 360 timing | 610 NPICA |

Each archive contains the exact saved R inputs in data/ and original fit records in fits/.
The final Section 6 direction scores are in section6.csv. The initial angle fields in the
raw fit records predate the corrected matrix reshaping. code/score_section6.R reproduces
the final scores, checks latent signals, and compares them to the original Stage2 results.
Section 6.3 has ten repeated native process exits. Failed attempts count in the denominator;
timing medians are over returned fits. Section 6.2 is excluded. Convergence flags are
distinct from returned fits throughout. Historical method rows in the two plot-data files
are the reference methods required to reproduce the published figures.

Use R 4.4.1 with the original repository library plus MSclust 1.0.4 and icamix 1.0.6.
Run Rscript scripts/setup_runtime.R from the repository root to install pinned packages.
Other dependencies include mclust, digest, clue, and ggplot2. The numerical implementation
is from v2.0.1-direct-residual-jacobi. Original method seeds, input hashes, versions and
raw objects are preserved. The archive portable entry points use relative repository paths.

Examples (replace Rscript by the full executable path if necessary):

    Rscript --vanilla code/figures.R
    Rscript --vanilla code/unpack.R section6
    Rscript --vanilla code/score_section6.R
    Rscript --vanilla code/unpack.R section4
    Rscript --vanilla code/refit.R section4 1

The refit command accepts section4, section5 or section6 and a one-based row number in
the corresponding table. It loads the original input, checks the semantic hash and
uses the recorded seed. New outputs go into refits/. For native crash isolation,
launch each requested row in its own R process. Section 6.3 original benchmarks used
eight concurrent fits; individual replay times are environment dependent.
The scoring script writes derived tables under replay/section6. Figure regeneration
writes the four PDFs in Stage2/figures/.

MMST uses the MSclust core with a vectorized equivalent direction objective; NPICA uses
the icamix native core. MMSCN and PMSCN are independent implementations of Tortora et al.
All four use one fit after a five-start k-means initialization. Tortora et al. used PAM.
See validation/ for the parity and numerical checks for these four methods.

Download section4.zip, section5.zip or section6.zip and SHA256SUMS.txt from the [revision release](https://github.com/xiuchuan18/robust-cluster/releases/tag/four-competitors-v1). Verify the digest and place the ZIP under Stage2/CompetitorExtension/archives/ before unpacking. Public records omit local machine provenance paths.
