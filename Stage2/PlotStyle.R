# PlotStyle.R
# Shared journal-ready styling for every manuscript figure.

JOURNAL_METHOD_COLORS <- c(
  "GMM"       = "#4E79A7",
  "tMM"       = "#76B7B2",
  "MSN"       = "#F28E2B",
  "MST"       = "#E15759",
  "MSGHD"     = "#B07AA1",
  "CGMM"      = "#9C755F",
  "MixSIT-VV" = "#009E73",
  "MixSIT-EV" = "#0072B2",
  "MixSIT-EE" = "#C89B00"
)

JOURNAL_METHOD_SHAPES <- c(
  "GMM" = 16, "tMM" = 17, "MSN" = 15, "MST" = 18, "MSGHD" = 8,
  "CGMM" = 3, "MixSIT-VV" = 0, "MixSIT-EV" = 1, "MixSIT-EE" = 2
)

JOURNAL_METHOD_LINETYPES <- c(
  "GMM" = "solid", "tMM" = "dashed", "MSN" = "dotdash",
  "MST" = "longdash", "MSGHD" = "twodash", "CGMM" = "dotted",
  "MixSIT-VV" = "dotted", "MixSIT-EV" = "dashed", "MixSIT-EE" = "solid"
)

theme_journal <- function(base_size = 8.5, base_family = "sans") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      line = ggplot2::element_line(linewidth = 0.35, colour = "#303030"),
      rect = ggplot2::element_rect(linewidth = 0.35, colour = "#303030"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "#E6E6E6"),
      strip.background = ggplot2::element_rect(fill = "#F0F0F0", colour = "#303030"),
      strip.text = ggplot2::element_text(face = "bold", margin = ggplot2::margin(3, 3, 3, 3)),
      axis.title = ggplot2::element_text(colour = "#202020"),
      axis.text = ggplot2::element_text(colour = "#303030"),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(1, 1, 1, 1),
      legend.box.spacing = grid::unit(2, "pt"),
      plot.margin = ggplot2::margin(5, 6, 5, 5)
    )
}

scale_colour_methods <- function(...) {
  ggplot2::scale_colour_manual(values = JOURNAL_METHOD_COLORS, ...)
}

scale_fill_methods <- function(...) {
  ggplot2::scale_fill_manual(values = JOURNAL_METHOD_COLORS, ...)
}

scale_shape_methods <- function(...) {
  ggplot2::scale_shape_manual(values = JOURNAL_METHOD_SHAPES, ...)
}

scale_linetype_methods <- function(...) {
  ggplot2::scale_linetype_manual(values = JOURNAL_METHOD_LINETYPES, ...)
}

journal_save <- function(filename, plot, width = 7.2, height = 4.8) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename, plot, width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, bg = "white"
  )
}
