# Local compatibility patch

This directory contains the GPL-licensed CRAN source of EMMIXskew 1.0.3.
Three C source files (`src/density.c`, `src/distance.c`, and `src/module.c`)
were changed only to supply the hidden Fortran character-length arguments
required by current R toolchains. The statistical algorithms, defaults, and R
interface are unchanged. `MD5` records the checksums of the patched sources.

The unmodified release is archived at
<https://cran.r-project.org/src/contrib/Archive/EMMIXskew/>.
