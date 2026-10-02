# Final figure set

`Stage2/figures/` contains only the six final one-page vector PDFs intended
for manuscript or supplementary-material use. All are regenerated from the
tracked aggregate results; no values are edited at the graphics stage.

## Shared design

- Maximum width: 7.2 inches (full manuscript text width).
- Typeface: device-independent sans serif through Cairo PDF.
- Method colours, point shapes, and line types: `Stage2/PlotStyle.R`.
- Panel strips, grid weight, margins, legends, and base type size are shared.
- Colour is never the only distinction in line plots; shapes and line types are
  also encoded. Boxplots retain the shared colour mapping and fixed method
  order.
- Figure titles and prose captions are intentionally excluded from the PDFs;
  those belong in the manuscript caption and cannot become inconsistent with
  the paper.
- Sparse designs are not forced into line-plot panels: the phase map uses a
  condition-by-method matrix.
- Banknote is integrated with the other three datasets in the common
  contamination figure rather than presented as an artificial standalone case.

## Generating scripts

| Final PDF | Generating script |
|---|---|
| `Figure1_ARI_Boxplots.pdf` | `Stage2/Section4_Figures.R` |
| `Figure2_BIC_Boxplots.pdf` | `Stage2/Section4_Figures.R` |
| `Section4_Wbar_Density.pdf` | `Stage2/Section4_WbarAnalysis.R` |
| `Section5_CellWise05_ARI.pdf` | `Stage2/Section5_CellWise05_Figure.R` |
| `Section6_Scaling.pdf` | `Stage2/Section6_ScalingFigure.R` |
| `calibrated-flexibility/CalFlex_PhaseMap.pdf` | `Stage2/CalFlex_PhaseMapAnalysis.R` |

Run `tests/figure_integrity.R` after rendering. It rejects extra files,
missing files, invalid PDF signatures, implausibly small outputs, and multi-page
figures when Poppler's `pdfinfo` is available.
