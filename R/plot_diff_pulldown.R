## plot_diff_pulldown.R
## Visualisation functions for DifferentialPulldownAnalysis results.
##
## Implements all gene-level, isoform-level, component-level, and integrative
## plots for exploring how compared conditions differ at the isoform interaction
## level.
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: plot_cross_gene_overview, plot_cross_gene_scatter,
##           plot_cross_gene_boxplot_overview, plot_cross_gene_boxplot

#' @include diff_pulldown_analysis_class.R
NULL
#' @include diff_pulldown_gene_class.R
NULL
#' @include component_analysis.R
NULL

utils::globalVariables(c("group", "delta__", "label__", "x__", "y__"))


# ==============================================================================
#  Internal helpers
# ==============================================================================

#' Check that ggplot2 and required extension packages are available
#' @noRd
.check_diff_plot_pkgs <- function(extra = character()) {
  required <- c("ggplot2", extra)
  for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf(
        "Package '%s' is required. Install with: install.packages('%s').",
        pkg, pkg
      ), call. = FALSE)
  }
}

#' Shared ggplot2 theme for differential pulldown plots
#' @noRd
.diff_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 13, face = "bold",
                                              margin = ggplot2::margin(b = 4)),
      plot.subtitle   = ggplot2::element_text(size = 9, colour = "grey40",
                                              margin = ggplot2::margin(b = 8)),
      plot.caption    = ggplot2::element_text(size = 8, colour = "grey55",
                                              hjust = 0),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right",
      legend.title     = ggplot2::element_text(size = 9, face = "bold"),
      legend.text      = ggplot2::element_text(size = 8),
      strip.text       = ggplot2::element_text(size = 10, face = "bold")
    )
}

#' CVD-safe direction colours (Okabe-Ito)
.col_up   <- "#E69F00"   # orange  - up condition
.col_down <- "#56B4E9"   # sky blue - down condition

#' Validate that a DifferentialPulldownAnalysis object has usable data
#' @noRd
.check_diff_analysis <- function(obj) {
  if (!is(obj, "DifferentialPulldownAnalysis"))
    stop("`obj` must be a DifferentialPulldownAnalysis object.", call. = FALSE)
  if (nrow(obj@cross_gene_summary) == 0L)
    stop("cross_gene_summary is empty - no genes to plot.", call. = FALSE)
  invisible(TRUE)
}


# ==============================================================================
#  Exported: plot_cross_gene_scatter, plot_cross_gene_boxplot
# ==============================================================================

#' Internal: select genes to label from a data frame
#'
#' Collects candidates from each requested strategy, then trims to the top
#' \code{label_n} overall ranked by absolute delta.  This guarantees a hard
#' maximum of \code{label_n} labels regardless of how many strategies are
#' requested.
#'
#' Scatter strategies: \code{"delta"} (furthest from diagonal),
#' \code{"x"} (highest x value), \code{"y"} (highest y value).
#'
#' Boxplot strategies: \code{"delta"} (largest absolute delta),
#' \code{"up"} (largest positive delta - most orange),
#' \code{"down"} (largest negative delta - most blue).
#'
#' @param df       Data frame with columns \code{gene_id}, \code{x__},
#'   \code{y__}, \code{delta__}.
#' @param label    Character vector of strategy names.
#' @param label_n  Integer. Maximum total genes to label.
#' @return Character vector of gene IDs to label.
#' @noRd
.select_labels <- function(df, label, label_n) {
  if (label_n <= 0L || length(label) == 0L) return(character(0))

  # Collect candidate gene IDs from each strategy using a generous pool
  # (up to 2x label_n per strategy) then trim the union to label_n overall
  pool_n <- max(label_n * 2L, 10L)
  ids <- character(0)

  # Scatter strategies
  if ("delta" %in% label)
    ids <- union(ids, df$gene_id[
      rank(-abs(df$delta__), ties.method = "first") <= pool_n])
  if ("x" %in% label)
    ids <- union(ids, df$gene_id[
      rank(-df$x__, ties.method = "first") <= pool_n])
  if ("y" %in% label)
    ids <- union(ids, df$gene_id[
      rank(-df$y__, ties.method = "first") <= pool_n])

  # Boxplot strategies
  if ("up" %in% label)
    ids <- union(ids, df$gene_id[
      rank(-df$delta__, ties.method = "first") <= pool_n])
  if ("down" %in% label)
    ids <- union(ids, df$gene_id[
      rank(df$delta__, ties.method = "first") <= pool_n])

  if (length(ids) == 0L) return(character(0))

  # Trim to label_n by ranking candidates on absolute delta as tiebreaker
  candidates <- df[df$gene_id %in% ids, ]
  top <- candidates$gene_id[
    rank(-abs(candidates$delta__), ties.method = "first") <= label_n]
  top
}


#' Internal: build one cross-gene scatter ggplot without printing
#'
#' Shared core used by both \code{plot_cross_gene_scatter} and
#' \code{plot_cross_gene_overview}.  Returns the ggplot object silently;
#' the exported wrapper handles printing.
#'
#' @noRd
.build_cross_gene_scatter <- function(cgs, analysis,
                                      x_col, y_col,
                                      x_axis_log, y_axis_log,
                                      x_log_offset = 1,
                                      y_log_offset = 1,
                                      label_n, label,
                                      title,
                                      panel_theme = FALSE) {

  for (col in c(x_col, y_col)) {
    if (!col %in% names(cgs))
      stop(sprintf(
        "Column '%s' not found in cross_gene_summary.\nAvailable: %s",
        col, paste(names(cgs), collapse = ", ")
      ), call. = FALSE)
    if (!is.numeric(cgs[[col]]))
      stop(sprintf("Column '%s' must be numeric.", col), call. = FALSE)
  }

  df <- cgs[!is.na(cgs[[x_col]]) & !is.na(cgs[[y_col]]), ]
  if (nrow(df) == 0L) return(NULL)

  x_vals <- df[[x_col]]
  y_vals <- df[[y_col]]
  # offset=1: log10(x+1) keeps zeros at origin; offset=0: plain log10 (for
  # columns like nodes where minimum is 1 and zeros never occur)
  if (x_axis_log) {
    x_vals[x_vals + x_log_offset <= 0] <- NA
    x_vals <- log10(x_vals + x_log_offset)
  }
  if (y_axis_log) {
    y_vals[y_vals + y_log_offset <= 0] <- NA
    y_vals <- log10(y_vals + y_log_offset)
  }
  df$x__    <- x_vals
  df$y__    <- y_vals
  df$delta__ <- df$x__ - df$y__

  df <- df[!is.na(df$x__) & !is.na(df$y__), ]
  if (nrow(df) == 0L) return(NULL)

  label_n   <- as.integer(label_n)
  label_ids <- .select_labels(df, label, label_n)
  df$label__ <- ifelse(df$gene_id %in% label_ids, df$gene_id, NA_character_)

  all_vals <- c(df$x__, df$y__)
  ax_max   <- max(all_vals, na.rm = TRUE) * 1.05
  ax_min   <- min(all_vals, na.rm = TRUE)
  ax_min   <- if (ax_min > 0) 0 else ax_min * 1.05

  # Short readable axis labels - replace underscores, append scaling suffix
  .short_ax <- function(col, is_log, offset) {
    readable <- list(
      up_nodes       = "up nodes",       down_nodes       = "down nodes",
      up_edges       = "up edges",       down_edges       = "down edges",
      up_density     = "up density",     down_density     = "down density",
      up_avg_degree  = "up avg degree",  down_avg_degree  = "down avg degree",
      up_clustering  = "up clustering",  down_clustering  = "down clustering",
      n_edges_gained = "edges gained",   n_edges_lost     = "edges lost"
    )
    label <- readable[[col]] %||% gsub("_", " ", col)
    if (is_log) {
      suffix <- if (offset == 0) " (log)" else " (log+1)"
      paste0(label, suffix)
    } else label
  }
  x_lab <- .short_ax(x_col, x_axis_log, x_log_offset)
  y_lab <- .short_ax(y_col, y_axis_log, y_log_offset)

  auto_title <- if (!is.null(title)) title else {
    # Short human-readable titles for the overview panels
    pair_labels <- list(
      "up_nodes-down_nodes"           = "Number of nodes",
      "up_edges-down_edges"           = "Number of edges",
      "up_density-down_density"       = "Density",
      "up_avg_degree-down_avg_degree" = "Avg degree",
      "up_clustering-down_clustering" = "Clustering",
      "n_edges_gained-n_edges_lost"   = "Edges gained vs lost"
    )
    key <- paste0(x_col, "-", y_col)
    pair_labels[[key]] %||% sprintf("%s vs %s", x_col, y_col)
  }

  base_theme <- if (panel_theme) {
    # Lighter theme for small overview panels
    ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(size = 9, face = "bold",
                                                 margin = ggplot2::margin(b = 2)),
        plot.caption     = ggplot2::element_text(size = 7, colour = "grey55",
                                                 hjust = 0),
        panel.grid.minor = ggplot2::element_blank(),
        axis.title       = ggplot2::element_text(size = 8),
        axis.text        = ggplot2::element_text(size = 7)
      )
  } else {
    .diff_plot_theme()
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x__, y = y__, label = label__)) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", colour = "grey70", linewidth = 0.4
    ) +
    ggplot2::geom_point(colour = "grey30", size = 2, alpha = 0.75) +
    ggrepel::geom_text_repel(
      na.rm = TRUE, size = 2.4, colour = "grey20",
      max.overlaps = 15, box.padding = 0.25
    ) +
    ggplot2::coord_fixed(xlim = c(ax_min, ax_max), ylim = c(ax_min, ax_max)) +
    ggplot2::labs(
      title   = auto_title,
      x       = x_lab,
      y       = y_lab,
      caption = paste0("n = ", nrow(df))
    ) +
    base_theme

  p
}


#' Plot cross-gene summary: scatter of two columns per gene
#'
#' Plots any two numeric columns from \code{cross_gene_summary} as a scatter,
#' one point per gene.  A diagonal identity line makes condition asymmetry
#' immediately visible.  For a quick five-panel overview across all standard
#' metrics, see \code{\link{plot_cross_gene_overview}}.
#'
#' @param obj       An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis  Character scalar.  Label of the
#'   \linkS4class{DifferentialPulldownAnalysis}; use \code{list_analyses(obj)}
#'   to see available labels.
#' @param x_col     Character scalar.  Column for the X-axis (must be numeric
#'   in \code{cross_gene_summary}).
#' @param y_col     Character scalar.  Column for the Y-axis (must be numeric
#'   in \code{cross_gene_summary}).
#' @param x_axis_log Logical.  Apply \code{log10(x + 1)} scaling to the
#'   X-axis?  Default \code{FALSE}.  Zeros map to 0 and are retained;
#'   the axis label will read \code{log10(<col> + 1)}.
#' @param y_axis_log Logical.  Apply \code{log10(x + 1)} scaling to the
#'   Y-axis?  Default \code{FALSE}.
#' @param x_log_offset Numeric scalar. Offset added before \code{log10()} on
#'   the x-axis when \code{x_axis_log = TRUE}. Default \code{1}, giving
#'   \code{log10(x + 1)}.
#' @param y_log_offset Numeric scalar. Offset added before \code{log10()} on
#'   the y-axis when \code{y_axis_log = TRUE}. Default \code{1}.
#' @param label_n   Integer.  Maximum total genes to label across all strategies.
#'   Default \code{10}.  Set \code{0} to suppress all labels.
#' @param label     Character vector.  One or more of \code{"delta"}
#'   (genes furthest from the identity diagonal), \code{"x"} (highest
#'   X values), \code{"y"} (highest Y values).  Strategies are unioned, so
#'   \code{c("delta", "x")} labels both.  Default \code{"delta"}.
#' @param title     Character scalar or \code{NULL} (auto-generated).
#'
#' @return A \code{ggplot} object, invisibly.
#' @export
plot_cross_gene_scatter <- function(obj,
                                    analysis,
                                    x_col,
                                    y_col,
                                    x_axis_log   = FALSE,
                                    y_axis_log   = FALSE,
                                    x_log_offset = 1,
                                    y_log_offset = 1,
                                    label_n      = 10L,
                                    label        = "delta",
                                    title        = NULL) {
  .check_diff_plot_pkgs(c("ggrepel"))
  diff_res <- .resolve_diff_analysis(obj, analysis)
  label    <- match.arg(label, c("delta", "x", "y"), several.ok = TRUE)

  p <- .build_cross_gene_scatter(
    cgs          = diff_res@cross_gene_summary,
    analysis     = analysis,
    x_col        = x_col,
    y_col        = y_col,
    x_axis_log   = x_axis_log,
    y_axis_log   = y_axis_log,
    x_log_offset = x_log_offset,
    y_log_offset = y_log_offset,
    label_n      = label_n,
    label        = label,
    title        = title,
    panel_theme  = FALSE
  )

  if (is.null(p))
    stop("No rows with non-NA values in both columns.", call. = FALSE)

  # Override auto-title with full analysis label when called directly
  if (is.null(title))
    p <- p + ggplot2::labs(
      title   = sprintf("%s vs %s  [%s]", x_col, y_col, analysis),
      caption = paste0("n = ", nrow(p$data), " genes  |  analysis: ", analysis)
    )

  p
}


# ==============================================================================
#  Overview plot - five-panel 2x3 grid
# ==============================================================================

#' Overview scatter: six standard up-vs-down comparisons in a 2x3 grid
#'
#' A first-look function that assembles six \code{plot_cross_gene_scatter}
#' panels in a 2x3 grid:
#' \enumerate{
#'   \item Nodes (up vs down)
#'   \item Edges (up vs down)
#'   \item Density (up vs down)
#'   \item Average degree (up vs down)
#'   \item Clustering coefficient (up vs down)
#'   \item Edges gained vs edges lost
#' }
#' Count columns (nodes, edges, edges gained/lost) use \code{log10(x + 1)}
#' scaling by default so that zero-interaction genes stay at the origin
#' rather than being dropped.  Labels are suppressed so the overview stays
#' clean; use \code{\link{plot_cross_gene_scatter}} directly on any panel
#' that catches your eye to add labels and explore further.
#'
#' Requires the \pkg{patchwork} package
#' (\code{install.packages("patchwork")}).
#'
#' @param obj        An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis   Character scalar.  Label of the
#'   \linkS4class{DifferentialPulldownAnalysis}; use \code{list_analyses(obj)}
#'   to see available labels.
#' @param log_counts Logical.  Apply \code{log10(x + 1)} scaling to the
#'   nodes, edges, and edges-gained/lost panels?  Default \code{TRUE}.
#'   Zeros map to 0 and are retained on the plot.
#' @param label_n    Integer.  Maximum total genes to label per panel.
#'   Default \code{0} (no labels - keeps the overview clean).  Set to a
#'   small value like \code{3} to flag the most extreme genes.
#' @param label      Character vector.  Labelling strategy passed to each
#'   panel.  One or more of \code{"delta"}, \code{"x"}, \code{"y"}.
#'   Only used when \code{label_n > 0}.  Default \code{"delta"}.
#' @param title      Character scalar or \code{NULL}.  Overall title for the
#'   combined figure.  Defaults to
#'   \code{"Cross-gene overview  [<analysis>]"}.
#'
#' @return A \code{patchwork} object, invisibly.  The combined figure is
#'   printed as a side-effect.
#' @export
plot_cross_gene_overview <- function(obj,
                                     analysis,
                                     log_counts = TRUE,
                                     label_n    = 0L,
                                     label      = "delta",
                                     title      = NULL) {
  .check_diff_plot_pkgs(c("ggrepel", "patchwork"))
  diff_res <- .resolve_diff_analysis(obj, analysis)
  cgs      <- diff_res@cross_gene_summary
  label    <- match.arg(label, c("delta", "x", "y"), several.ok = TRUE)

  # Panel definitions: list(x_col, y_col, x_log, y_log, x_offset, y_offset)
  # Row 1: count/magnitude metrics (Nodes, Edges, Avg degree)
  # Row 2: normalised/ratio metrics (Density, Clustering, Edges gained vs lost)
  # Nodes use offset=0 (plain log10, min value is 1 so zeros never occur)
  # Edges/degree/gained/lost use offset=1 (log10(x+1), zeros are valid)
  panels <- list(
    list(x = "up_nodes",       y = "down_nodes",     xl = log_counts, yl = log_counts, xo = 0, yo = 0),
    list(x = "up_edges",       y = "down_edges",     xl = log_counts, yl = log_counts, xo = 1, yo = 1),
    list(x = "up_avg_degree",  y = "down_avg_degree",xl = log_counts, yl = log_counts, xo = 1, yo = 1),
    list(x = "up_density",     y = "down_density",   xl = FALSE,      yl = FALSE,      xo = 1, yo = 1),
    list(x = "up_clustering",  y = "down_clustering",xl = FALSE,      yl = FALSE,      xo = 1, yo = 1),
    list(x = "n_edges_gained", y = "n_edges_lost",   xl = log_counts, yl = log_counts, xo = 1, yo = 1)
  )

  plots <- lapply(panels, function(pd) {
    .build_cross_gene_scatter(
      cgs          = cgs,
      analysis     = analysis,
      x_col        = pd$x,
      y_col        = pd$y,
      x_axis_log   = pd$xl,
      y_axis_log   = pd$yl,
      x_log_offset = pd$xo,
      y_log_offset = pd$yo,
      label_n      = label_n,
      label        = label,
      title        = NULL,
      panel_theme  = TRUE
    )
  })

  # Drop any NULL panels (columns absent from cgs)
  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0L)
    stop("No panels could be built - check that cross_gene_summary has the ",
         "expected up_*/down_* columns.", call. = FALSE)

  n_genes <- nrow(cgs[!is.na(cgs$up_nodes) & !is.na(cgs$down_nodes), ])

  combined <- patchwork::wrap_plots(plots, ncol = 3) +
    patchwork::plot_annotation(
      title    = title %||% sprintf("Cross-gene overview  [%s]", analysis),
      subtitle = paste0(
        "Each panel: up-seed neighbourhood (X) vs down-seed neighbourhood (Y).  ",
        "Points on the diagonal = no change between conditions.\n",
        if (log_counts) paste0(
          "Axes labelled (log) use log10(x) scaling.  ",
          "Axes labelled (log+1) use log10(x + 1) scaling.  "
        ) else "",
        "n = ", n_genes, " genes with up and down subgraph."
      ),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 13, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 8,  colour = "grey40",
                                              lineheight = 1.3)
      )
    )

  combined
}


#' Plot cross-gene summary: boxplot of the difference between two columns
#'
#' For each gene computes a per-gene delta (\code{x_col - y_col}, or the
#' log10 fold-change \code{log10(x_col / y_col)} when \code{log_ratio = TRUE})
#' and displays the distribution as a boxplot with individual gene points
#' overlaid as a jitter strip.  A horizontal reference line at zero marks no
#' difference between conditions.  Orange points have a positive delta
#' (x_col > y_col); blue have a negative delta.
#'
#' @section Log-ratio mode:
#' When \code{log_ratio = TRUE} the plotted value is
#' \code{log10((x_col + 1) / (y_col + 1))}.  Adding 1 before dividing means
#' zero-value genes map to 0 (no change) rather than being dropped, so all
#' genes are retained.  The zero reference line corresponds to a ratio of 1
#' (no change between conditions).
#'
#' @param obj       An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis  Character scalar.  Label of the
#'   \linkS4class{DifferentialPulldownAnalysis}; use \code{list_analyses(obj)}
#'   to see available labels.
#' @param x_col     Character scalar.  Column for the minuend / numerator
#'   (e.g. the up-seed neighbourhood metric). Must be numeric.
#' @param y_col     Character scalar.  Column for the subtrahend / denominator
#'   (e.g. the down-seed neighbourhood metric). Must be numeric.
#'   Delta = \code{x_col - y_col} (or \code{log10(x_col/y_col)} when
#'   \code{log_ratio = TRUE}).
#' @param log_ratio Logical.  When \code{TRUE}, plot
#'   \code{log10((x_col + 1) / (y_col + 1))} instead of
#'   \code{x_col - y_col}.  Zeros map to 0 and no genes are dropped.
#'   Default \code{FALSE}.
#' @param label_n   Integer.  Maximum total genes to label across all strategies.
#'   Default \code{10}.  Set \code{0} to suppress all labels.
#' @param label     Character vector.  One or more of \code{"delta"}
#'   (largest absolute delta), \code{"up"} (largest positive delta, orange
#'   points), \code{"down"} (largest negative delta,
#'   largest negative delta, blue points).  Default \code{"delta"}.
#' @param title     Character scalar or \code{NULL} (auto-generated).
#'
#' @return A \code{ggplot} object, invisibly.
#' @export
plot_cross_gene_boxplot <- function(obj,
                                    analysis,
                                    x_col,
                                    y_col,
                                    log_ratio = FALSE,
                                    label_n   = 10L,
                                    label     = "delta",
                                    title     = NULL) {
  .check_diff_plot_pkgs(c("ggrepel"))
  diff_res <- .resolve_diff_analysis(obj, analysis)
  cgs      <- diff_res@cross_gene_summary

  label <- match.arg(label, c("delta", "up", "down"), several.ok = TRUE)

  for (col in c(x_col, y_col)) {
    if (!col %in% names(cgs))
      stop(sprintf(
        "Column '%s' not found in cross_gene_summary.\nAvailable: %s",
        col, paste(names(cgs), collapse = ", ")
      ), call. = FALSE)
    if (!is.numeric(cgs[[col]]))
      stop(sprintf("Column '%s' must be numeric.", col), call. = FALSE)
  }

  df <- cgs[!is.na(cgs[[x_col]]) & !is.na(cgs[[y_col]]), ]
  if (nrow(df) == 0L)
    stop("No rows with non-NA values in both columns.", call. = FALSE)

  # Keep raw values on df for "x"/"y" label ranking (always on original scale)
  df$x__ <- df[[x_col]]
  df$y__ <- df[[y_col]]

  if (log_ratio) {
    df$delta__ <- log10((df$x__ + 1) / (df$y__ + 1))
    delta_lab  <- sprintf("log10((%s+1)/(%s+1))", x_col, y_col)
  } else {
    df$delta__ <- df$x__ - df$y__
    delta_lab  <- sprintf("%s \u2212 %s", x_col, y_col)
  }

  # Use short readable label for known metric pairs
  short_labels <- list(
    "up_nodes-down_nodes"           = list(lin = "\u0394 nodes",      log = "log ratio (nodes)"),
    "up_edges-down_edges"           = list(lin = "\u0394 edges",      log = "log ratio (edges)"),
    "up_density-down_density"       = list(lin = "\u0394 density",    log = "log ratio (density)"),
    "up_avg_degree-down_avg_degree" = list(lin = "\u0394 avg degree", log = "log ratio (avg degree)"),
    "up_clustering-down_clustering" = list(lin = "\u0394 clustering", log = "log ratio (clustering)")
  )
  key <- paste0(x_col, "-", y_col)
  if (!is.null(short_labels[[key]])) {
    delta_lab <- if (log_ratio) short_labels[[key]]$log else short_labels[[key]]$lin
  }

  label_n   <- as.integer(label_n)
  label_ids <- .select_labels(df, label, label_n)
  df$label__ <- ifelse(df$gene_id %in% label_ids, df$gene_id, NA_character_)

  df$group  <- ""
  pt_colour <- ifelse(df$delta__ > 0, .col_up, .col_down)

  zero_label <- if (log_ratio) "ratio = 1 (no change)" else "no difference"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = delta__)) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype   = "dashed",
      colour     = "grey60",
      linewidth  = 0.5
    ) +
    ggplot2::geom_boxplot(
      width = 0.35, outlier.shape = NA,
      fill  = "grey92", colour = "grey40", linewidth = 0.5
    ) +
    ggplot2::geom_jitter(
      width  = 0.12, height = 0,
      size   = 2,    alpha  = 0.65,
      colour = pt_colour
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label__),
      na.rm = TRUE, size = 2.8, colour = "grey15",
      max.overlaps = 20, box.padding = 0.35,
      direction = "y", nudge_x = 0.25
    ) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(
      title    = title %||% sprintf(
        "%s  [%s]", delta_lab, analysis),
      subtitle = paste0(
        "Orange = ", x_col, " > ", y_col, ".  ",
        "Blue = ", x_col, " < ", y_col, ".  ",
        "Dashed line = ", zero_label, "."
      ),
      x       = NULL,
      y       = delta_lab,
      caption = paste0("n = ", nrow(df), " genes  |  analysis: ", analysis)
    ) +
    .diff_plot_theme() +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )

  p
}


# ==============================================================================
#  Boxplot overview - six panels in a single row
# ==============================================================================

#' Internal: build one cross-gene boxplot panel without printing
#' @noRd
.build_cross_gene_boxplot <- function(cgs, analysis,
                                      x_col, y_col,
                                      log_ratio,
                                      label_n, label,
                                      title,
                                      panel_theme = FALSE) {

  for (col in c(x_col, y_col)) {
    if (!col %in% names(cgs)) return(NULL)
    if (!is.numeric(cgs[[col]])) return(NULL)
  }

  df <- cgs[!is.na(cgs[[x_col]]) & !is.na(cgs[[y_col]]), ]
  if (nrow(df) == 0L) return(NULL)

  df$x__ <- df[[x_col]]
  df$y__ <- df[[y_col]]

  if (log_ratio) {
    # log10((x+1)/(y+1)) - zeros map to 0 (no change) rather than being dropped
    df$delta__ <- log10((df$x__ + 1) / (df$y__ + 1))
    delta_lab  <- sprintf("log10((%s+1)/(%s+1))", x_col, y_col)
  } else {
    df$delta__ <- df$x__ - df$y__
    delta_lab  <- sprintf("%s \u2212 %s", x_col, y_col)
  }

  # Short readable y-axis labels for known metric pairs
  short_labels <- list(
    "up_nodes-down_nodes"           = list(lin = "\u0394 nodes",      log = "log ratio (nodes)"),
    "up_edges-down_edges"           = list(lin = "\u0394 edges",      log = "log ratio (edges)"),
    "up_density-down_density"       = list(lin = "\u0394 density",    log = "log ratio (density)"),
    "up_avg_degree-down_avg_degree" = list(lin = "\u0394 avg degree", log = "log ratio (avg degree)"),
    "up_clustering-down_clustering" = list(lin = "\u0394 clustering", log = "log ratio (clustering)")
  )
  key <- paste0(x_col, "-", y_col)
  if (!is.null(short_labels[[key]])) {
    delta_lab <- if (log_ratio) short_labels[[key]]$log else short_labels[[key]]$lin
  }

  label_n   <- as.integer(label_n)
  label_ids <- .select_labels(df, label, label_n)
  df$label__ <- ifelse(df$gene_id %in% label_ids, df$gene_id, NA_character_)

  df$group  <- ""
  pt_colour <- ifelse(df$delta__ > 0, .col_up, .col_down)

  # Short auto-titles for overview panels
  auto_title <- if (!is.null(title)) title else {
    pair_labels <- list(
      "up_nodes-down_nodes"           = "Nodes",
      "up_edges-down_edges"           = "Edges",
      "up_density-down_density"       = "Density",
      "up_avg_degree-down_avg_degree" = "Avg degree",
      "up_clustering-down_clustering" = "Clustering"
    )
    key <- paste0(x_col, "-", y_col)
    pair_labels[[key]] %||% sprintf("%s \u2212 %s", x_col, y_col)
  }

  base_theme <- if (panel_theme) {
    ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        plot.title         = ggplot2::element_text(size = 8, face = "bold",
                                                   margin = ggplot2::margin(b = 2)),
        plot.caption       = ggplot2::element_text(size = 6, colour = "grey55",
                                                   hjust = 0),
        panel.grid.minor   = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        axis.title.x       = ggplot2::element_blank(),
        axis.text.x        = ggplot2::element_blank(),
        axis.title.y       = ggplot2::element_text(size = 7),
        axis.text.y        = ggplot2::element_text(size = 6)
      )
  } else {
    .diff_plot_theme() +
      ggplot2::theme(
        axis.text.x        = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank()
      )
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = delta__)) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype   = "dashed",
      colour     = "grey60",
      linewidth  = 0.4
    ) +
    ggplot2::geom_boxplot(
      width = 0.25, outlier.shape = NA,
      fill  = "grey92", colour = "grey40", linewidth = 0.4
    ) +
    ggplot2::geom_jitter(
      width  = 0.08, height = 0,
      size   = if (panel_theme) 1.2 else 2,
      alpha  = 0.65,
      colour = pt_colour
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label__),
      na.rm = TRUE,
      size  = if (panel_theme) 2.0 else 2.8,
      colour = "grey15",
      max.overlaps = 15, box.padding = 0.3,
      direction = "y", nudge_x = 0.15
    ) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(
      title   = auto_title,
      y       = delta_lab,
      caption = paste0("n = ", nrow(df))
    ) +
    base_theme

  p
}


#' Overview boxplot: five delta distributions in a single row
#'
#' Assembles five \code{plot_cross_gene_boxplot} panels in one row, each
#' showing the distribution of per-gene deltas (up - down) for a
#' standard metric pair:
#' \enumerate{
#'   \item Nodes
#'   \item Edges
#'   \item Density
#'   \item Average degree
#'   \item Clustering coefficient
#' }
#' All panels include \code{expand_limits(y = 0)} so the zero reference
#' line is always in frame, making it easy to compare the direction of
#' shift across metrics.  Labels are suppressed by default for a clean
#' overview; use \code{\link{plot_cross_gene_boxplot}} to drill into any
#' panel with labels.
#'
#' When \code{log_ratio = TRUE}, the delta is computed as
#' \code{log10((x+1)/(y+1))} so that zero-value genes map to 0 rather
#' than being dropped.
#'
#' Requires the \pkg{patchwork} package.
#'
#' @param obj        An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis   Character scalar.  Label of the
#'   \linkS4class{DifferentialPulldownAnalysis}; use \code{list_analyses(obj)}
#'   to see available labels.
#' @param log_ratio  Logical.  Use \code{log10((x+1)/(y+1))} instead of
#'   \code{x - y}?  Default \code{FALSE}.  When \code{TRUE}, zeros
#'   map to 0 (no change) and no genes are dropped.
#' @param label_n    Integer.  Maximum total genes to label per panel.  Default \code{0}.
#' @param label      Character vector.  One or more of \code{"delta"}
#'   (largest absolute delta), \code{"up"} (largest positive delta, orange
#'   points), \code{"down"} (largest negative delta,
#'   largest negative delta, blue points).  Default \code{"delta"}.
#' @param title      Character scalar or \code{NULL} (auto-generated).
#'
#' @return A \code{patchwork} object, invisibly.
#' @export
plot_cross_gene_boxplot_overview <- function(obj,
                                             analysis,
                                             log_ratio = FALSE,
                                             label_n   = 0L,
                                             label     = "delta",
                                             title     = NULL) {
  .check_diff_plot_pkgs(c("ggrepel", "patchwork"))
  diff_res <- .resolve_diff_analysis(obj, analysis)
  cgs      <- diff_res@cross_gene_summary
  label    <- match.arg(label, c("delta", "up", "down"), several.ok = TRUE)

  panels <- list(
    list(x = "up_nodes",      y = "down_nodes"),
    list(x = "up_edges",      y = "down_edges"),
    list(x = "up_avg_degree", y = "down_avg_degree"),
    list(x = "up_density",    y = "down_density"),
    list(x = "up_clustering", y = "down_clustering")
  )

  plots <- lapply(panels, function(pd) {
    .build_cross_gene_boxplot(
      cgs        = cgs,
      analysis   = analysis,
      x_col      = pd$x,
      y_col      = pd$y,
      log_ratio  = log_ratio,
      label_n    = label_n,
      label      = label,
      title      = NULL,
      panel_theme = TRUE
    )
  })

  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0L)
    stop("No panels could be built - check cross_gene_summary columns.",
         call. = FALSE)

  combined <- patchwork::wrap_plots(plots, nrow = 1) +
    patchwork::plot_annotation(
      title    = title %||% sprintf("Cross-gene delta overview  [%s]", analysis),
      subtitle = paste0(
        "Each panel: distribution of per-gene deltas (up \u2212 down).  ",
        "Orange = up > down.  Blue = down > up.  ",
        "Dashed line = no difference."
      ),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 13, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 9,  colour = "grey40")
      )
    )

  combined
}


# ==============================================================================
#  Internal: qualitative colour palette
# ==============================================================================

#' Generate a qualitative CVD-safe palette for component colouring
#' @noRd
.qualitative_palette <- function(n) {
  okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000"
  )
  if (n <= length(okabe_ito)) return(okabe_ito[seq_len(n)])
  # Recycle with a warning for large n
  warning(sprintf(
    "More than %d components - recycling colour palette.", length(okabe_ito)
  ), call. = FALSE)
  rep_len(okabe_ito, n)
}


# ==============================================================================
#  Null-coalescing operator (internal)
# ==============================================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b
