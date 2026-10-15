# Plot the summary.csv produced by test_mmseqs_params.py as a single
# size-vs-purity scatter
#
# Usage:
#   Rscript plot_mmseqs_results.R \
#     --summary results/summary.csv \
#     --outdir results/plots \
#     --helpers analysis/publication_plot_helpers.R
#
# Produces one figure, tradeoff.png: mean_cluster_size (x) vs
# weighted_purity (y), one point per parameter combination.
#   color = min_seq_id       (main lever for cluster size)
#   shape = cluster_mode x cov_mode combination (4 shapes)
#   size  = coverage_threshold (2 sizes)
# Combos you want sit toward the top-right: large AND pure.
#
# Requirements: ggplot2. An optional helper file can provide
# theme_publication().

# Kea palette, reordered for contrast: the two browns (#7B5C34, #AB7C47)
# sit next to the green (#6C803A) in luminance and read as near-identical
# at small point sizes, so they're pushed to the end. The first 5 —
# orange, charcoal, green, yellow, light grey — are the ones actually
# used for the 5 min_seq_id levels and are chosen to be maximally
# distinguishable from each other.
kea_palette <- c("#D73202", "#272318", "#6C803A", "#CCAE42", "#D3CDBF", "#7B5C34", "#AB7C47")

parse_args <- function(args) {
  defaults <- list(
    summary = NULL,
    outdir = NULL,
    helpers = NULL
  )
  for (i in seq_along(args)) {
    if (args[i] == "--summary") defaults$summary <- args[i + 1]
    if (args[i] == "--outdir") defaults$outdir <- args[i + 1]
    if (args[i] == "--helpers") defaults$helpers <- args[i + 1]
  }
  if (is.null(defaults$summary) || is.null(defaults$outdir)) {
    stop("Usage: Rscript plot_mmseqs_results.R --summary <path> --outdir <path> [--helpers <path>]")
  }
  defaults
}

load_plot_theme <- function(helper_path) {
  if (!is.null(helper_path) && file.exists(helper_path)) {
    source(helper_path)
  }
  if (!exists("theme_publication", mode = "function")) {
    theme_publication <<- function(legend_position = "right") {
      ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(legend.position = legend_position)
    }
  }
}

#' Single scatter of mean_cluster_size vs weighted_purity across the full
#' parameter sweep. min_seq_id -> color, (cluster_mode x cov_mode) -> shape,
#' coverage_threshold -> point size — no faceting, everything on one axis.
plot_size_purity_tradeoff <- function(df,
                                       size_col = "mean_cluster_size",
                                       purity_col = "mean_purity",
                                       title = "Cluster size vs. purity trade-off",
                                       save_path = NULL,
                                       width = 7.5,
                                       height = 5.5,
                                       dpi = 300,
                                       color_palette = kea_palette,
                                       point_sizes = c(2.4, 3.8, 5.2),
                                       theme = theme_publication(legend_position = "right")) {
  plot_df <- df[!is.na(df[[size_col]]) & df[[size_col]] != "" &
                  !is.na(df[[purity_col]]) & df[[purity_col]] != "", ]
  plot_df[[size_col]] <- suppressWarnings(as.numeric(plot_df[[size_col]]))
  plot_df[[purity_col]] <- suppressWarnings(as.numeric(plot_df[[purity_col]]))
  plot_df <- plot_df[!is.na(plot_df[[size_col]]) &
                       !is.na(plot_df[[purity_col]]), ]
  if (nrow(plot_df) == 0) {
    stop("No rows with numeric mean_cluster_size and weighted_purity values.")
  }

  plot_df$min_seq_id <- factor(plot_df$min_seq_id)
  plot_df$coverage_threshold <- factor(plot_df$coverage_threshold)
  plot_df$mode_combo <- interaction(plot_df$cluster_mode, plot_df$cov_mode,
                                     sep = " / ")

  n_colors <- nlevels(plot_df$min_seq_id)
  if (n_colors > length(color_palette)) {
    stop("More min_seq_id levels than palette colors — extend color_palette.")
  }
  n_shapes <- nlevels(plot_df$mode_combo)
  shape_values <- c(16, 17, 15, 18, 3, 8)
  if (n_shapes > length(shape_values)) {
    stop("More cluster_mode/cov_mode combinations than defined shape_values — extend shape_values.")
  }
  n_sizes <- nlevels(plot_df$coverage_threshold)
  if (n_sizes > length(point_sizes)) {
    stop("More coverage_threshold levels than defined point_sizes — extend point_sizes.")
  }

  p <- ggplot(plot_df, aes(x = .data[[size_col]], y = .data[[purity_col]],
                            color = min_seq_id,
                            shape = mode_combo,
                            size = coverage_threshold)) +
    geom_point(alpha = 0.85, stroke = 1.1) +
    scale_color_manual(values = color_palette[seq_len(n_colors)], name = "min_seq_id") +
    scale_shape_manual(values = shape_values[seq_len(n_shapes)], name = "cluster_mode / cov_mode") +
    scale_size_manual(values = point_sizes[seq_len(n_sizes)], name = "coverage_threshold") +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "mean cluster size", y = "mean purity", title = title) +
    theme

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(save_path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  }

  p
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  suppressPackageStartupMessages(library(ggplot2))
  load_plot_theme(args$helpers)

  df <- read.csv(args$summary, stringsAsFactors = FALSE)
  cat(sprintf("Read %d rows from %s\n", nrow(df), args$summary))

  plot_size_purity_tradeoff(
    df,
    save_path = file.path(args$outdir, "tradeoff.pdf")
  )

  cat(sprintf("Wrote %s\n", file.path(args$outdir, "tradeoff.pdf")))
}

main()