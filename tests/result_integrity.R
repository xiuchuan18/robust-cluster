args0 <- commandArgs(FALSE)
this <- sub("^--file=", "", args0[grep("^--file=", args0)])
root <- normalizePath(file.path(dirname(this), ".."), winslash = "/")
resdir <- file.path(root, "Stage2", "results")
read_result <- function(name) {
  f <- file.path(resdir, name)
  if (!file.exists(f)) stop("missing final result: ", name)
  readRDS(f)
}
unique_rows <- function(x, cols) {
  if (anyDuplicated(x[cols])) stop("duplicate keys: ", paste(cols, collapse = ", "))
}

# The final recovery scorer uses an exact assignment, not a greedy heuristic.
# This matrix is a counterexample on which row-wise greedy matching scores 16
# while the global assignment scores 22.
truth_env <- new.env(parent = globalenv())
sys.source(file.path(root, "Stage2", "Section6_TruthLib.R"), envir = truth_env)
benefit <- matrix(c(9, 8, 0,
                    8, 0, 0,
                    0, 7, 6), 3, byrow = TRUE)
stopifnot(identical(truth_env$solve_square_assignment(-benefit), c(2L, 1L, 3L)))

s4 <- read_result("Section4_Simulation_Results.rds")
stopifnot(length(s4) == 6L, all(vapply(s4, nrow, integer(1)) == 100L))
for (x in s4) unique_rows(x, "sim_id")

# The S5 latent-scale diagnostic must be invariant to latent-axis relabelling.
# Every component carries a one-to-one map from observed coordinates to fitted
# latent axes, and the published summary must use the matched weights.
wbar <- read_result("Section4_Wbar_S7_CellWise10.rds")
wbar_version <- "global-observed-direction-v1"
stopifnot(length(wbar) == 100L,
          all(vapply(wbar, function(x)
            identical(x$axis_match_version, wbar_version), logical(1))))
for (rec in wbar) for (v in c("MixSIT_None", "MixSIT_U", "MixSIT_Full")) {
  z <- rec$wbar[[v]]
  stopifnot(is.list(z), identical(z$matching, wbar_version),
            identical(dim(z$weights), c(700L, 10L)),
            identical(dim(z$axis_map), c(3L, 10L)),
            identical(dim(z$axis_abs_cosine), c(3L, 10L)),
            all(is.finite(z$weights)), all(is.finite(z$axis_abs_cosine)),
            all(z$axis_abs_cosine >= 0 & z$axis_abs_cosine <= 1))
  if (!all(apply(z$axis_map, 1, function(p)
      identical(sort(as.integer(p)), seq_len(10L)))))
    stop("Wbar contains a non-permutation axis assignment")
}
wbar_summary <- read.csv(file.path(resdir, "Section4_Wbar_Summary.csv"))
stopifnot(nrow(wbar_summary) == 3L, all(wbar_summary$n_rep == 100L),
          all(wbar_summary$auc > 0.8),
          all(wbar_summary$mean_abs_cosine > 0.9))

s5 <- read_result("Section5_RealData.rds")
stopifnot(nrow(s5) == 27L, setequal(unique(s5$dataset), c("Biopsy", "AIS", "Wine")))
unique_rows(s5, c("dataset", "method"))

s5c <- read_result("Section5_CellWise05.rds")
stopifnot(nrow(s5c) == 1350L, all(s5c$alpha == 0.05), all(s5c$eta == 5))
unique_rows(s5c, c("dataset", "rep", "method"))

checks <- list(
  Section6_Recovery.rds = c(22500L, "N", "variant", "rep", "k", "axis"),
  Section6_LambdaSensitivity.rds = c(8400L, "lam", "N", "rep", "k", "axis"),
  Section6_KSelect.rds = c(7200L, "cell", "method", "rep", "K"),
  Section6_Scale.rds = c(1800L, "N", "m", "K", "method", "rep"))
for (nm in names(checks)) {
  x <- read_result(nm); spec <- checks[[nm]]
  if (nrow(x) != as.integer(spec[1])) stop(nm, " has ", nrow(x), " rows")
  unique_rows(x, spec[-1])
}

for (nm in c("Section6_Recovery.rds", "Section6_LambdaSensitivity.rds")) {
  x <- read_result(nm)
  required <- c("matched_component", "matched_axis", "component_matching",
                "axis_matching")
  if (!all(required %in% names(x))) stop(nm, " lacks global-matching audit fields")
  stopifnot(all(x$component_matching == "global_maximum_overlap"),
            all(x$axis_matching == "global_minimum_total_angle"))
  fit_keys <- if ("variant" %in% names(x)) c("design", "N", "rep", "variant")
              else c("design", "N", "rep", "lam")
  component_rows <- unique(x[c(fit_keys, "k", "matched_component")])
  component_groups <- split(component_rows$matched_component,
                            interaction(component_rows[fit_keys], drop = TRUE))
  if (!all(vapply(component_groups, function(z)
    identical(sort(as.integer(z)), seq_along(z)), logical(1))))
    stop(nm, " contains a non-permutation component assignment")
  axis_groups <- split(x$matched_axis,
                       interaction(x[c(fit_keys, "k")], drop = TRUE))
  if (!all(vapply(axis_groups, function(z)
    identical(sort(as.integer(z)), seq_along(z)), logical(1))))
    stop(nm, " contains a non-permutation axis assignment")
}

cat("final result integrity passed\n")
