## plot_network_diagnostics.R
## Plotting function for NetworkDiagnostics results.
##
## Produces a 2-row x 3-column panel figure suitable for thesis and
## publication use. No in-plot titles (Row/panel titles live in the LaTeX
## caption); panels are tagged A-F via patchwork. Panels are grouped and
## coloured by network level: isoform-level panels (top row) share one
## accent colour, gene-level panels (bottom row, left two) share a second
## accent colour, and the categorical connectivity-status panel (bottom
## right) uses a 3-colour categorical scheme.
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: plot_network_diagnostics

#' @include isonetanalyzer_class.R
NULL
#' @include network_diagnostics_class.R
NULL

# Suppress R CMD check NOTEs for ggplot2 column names used in aes()
utils::globalVariables(c(
  "degree",
  "clustering_coefficient",
  "n_nodes",
  "status",
  "n_genes",
  "n_isoforms_total",
  "log10_x",
  "n"
))

# Semantic accent colours (Okabe-Ito palette)
.ISOFORM_COLOUR <- "#0072B2"   # blue  - isoform-level panels
.GENE_COLOUR    <- "#E69F00"   # orange - gene-level panels


#' Bin a strictly-positive numeric vector in log10 space
#'
#' Standard fixed-width histogram binning done on log10(values), so the
#' result can be plotted on an ordinary linear x-axis (with custom power-of-
#' ten tick labels) using a single constant bin width. This avoids both the
#' zero-count log(0) problem (empty bins are simply dropped) and the
#' variable relative-width approach, which is fragile when combined with
#' scale_x_log10().
#' @noRd
.log_hist_table <- function(values, bins = 30L) {
  lv <- log10(values)
  rng <- range(lv)
  if (diff(rng) == 0) rng <- rng + c(-0.5, 0.5)   # guard: all-identical values
  breaks <- seq(rng[1], rng[2], length.out = bins + 1L)
  h <- graphics::hist(lv, breaks = breaks, plot = FALSE)
  out <- data.frame(
    log10_x    = h$mids,
    n          = h$counts,
    bin_width  = diff(h$breaks)[1]
  )
  out[out$n > 0, , drop = FALSE]   # drop empty bins; nothing left to log(0)
}


#' Build integer power-of-ten breaks/labels spanning a log10-scale range
#'
#' Labels are plain decimal numbers (1, 10, 100, ...), matching the default
#' label style ggplot2 produces for scale_x_log10()/scale_y_log10() in the
#' other panels, rather than exponent notation.
#' @noRd
.log10_breaks_labels <- function(log10_range) {
  brks <- seq(floor(log10_range[1]), ceiling(log10_range[2]))
  list(breaks = brks, labels = format(10^brks, scientific = FALSE))
}


# ==============================================================================
#  Internal panel helpers
# ==============================================================================

#' Shared minimal ggplot2 theme for all diagnostic panels
#' @noRd
.diag_theme <- function() {
  ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.title       = ggplot2::element_text(size = 9),
      axis.text        = ggplot2::element_text(size = 8),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      panel.grid.minor = ggplot2::element_blank()
    )
}


#' Panel A — Isoform degree distribution (isoform-level, top row)
#' @noRd
.panel_isoform_degree <- function(diag) {
  df <- .df(diag@isoform_degree_table)
  df <- df[df$degree > 0L, , drop = FALSE]   # exclude isolated

  ht <- .log_hist_table(df$degree, bins = 30L)
  axis_info <- .log10_breaks_labels(range(ht$log10_x))
  thr_log   <- log10(diag@sticky_isoform_threshold)

  ggplot2::ggplot(ht, ggplot2::aes(x = log10_x, y = n)) +
    ggplot2::geom_col(
      fill   = .ISOFORM_COLOUR,
      colour = "white",
      linewidth = 0.2,
      width  = ht$bin_width
    ) +
    ggplot2::scale_x_continuous(breaks = axis_info$breaks, labels = axis_info$labels) +
    ggplot2::scale_y_log10() +
    ggplot2::geom_vline(
      xintercept = thr_log,
      linetype   = "dashed",
      colour     = "grey20",
      linewidth  = 0.6
    ) +
    ggplot2::annotate(
      "text",
      x      = thr_log,
      y      = Inf,
      hjust  = -0.15,
      vjust  = 1.4,
      size   = 2.8,
      colour = "grey20",
      label  = "sticky\nthreshold"
    ) +
    ggplot2::labs(
      x = "Degree (log scale)",
      y = "Number of isoforms (log scale)"
    ) +
    .diag_theme()
}


#' Panel D — Gene degree distribution (gene-level, bottom row)
#' @noRd
.panel_gene_degree <- function(diag) {
  df <- .df(diag@gene_degree_table)
  df <- df[df$degree > 0L, , drop = FALSE]   # exclude isolated, consistent with panel A

  ht <- .log_hist_table(df$degree, bins = 30L)
  axis_info <- .log10_breaks_labels(range(ht$log10_x))
  thr_log   <- log10(diag@sticky_gene_threshold)

  ggplot2::ggplot(ht, ggplot2::aes(x = log10_x, y = n)) +
    ggplot2::geom_col(
      fill      = .GENE_COLOUR,
      colour    = "white",
      linewidth = 0.2,
      width     = ht$bin_width
    ) +
    ggplot2::scale_x_continuous(breaks = axis_info$breaks, labels = axis_info$labels) +
    ggplot2::scale_y_log10() +
    ggplot2::geom_vline(
      xintercept = thr_log,
      linetype   = "dashed",
      colour     = "grey20",
      linewidth  = 0.6
    ) +
    ggplot2::annotate(
      "text",
      x      = thr_log,
      y      = Inf,
      hjust  = -0.15,
      vjust  = 1.4,
      size   = 2.8,
      colour = "grey20",
      label  = "sticky\nthreshold"
    ) +
    ggplot2::labs(
      x = "Degree (log scale)",
      y = "Number of genes (log scale)"
    ) +
    .diag_theme()
}


#' Panel B — Connected component sizes (isoform-level, top row)
#' @noRd
.panel_components <- function(diag) {
  df <- .df(diag@component_table)

  # Show top 20 components; collapse remainder into an "other" bar
  top_n  <- 20L
  n_comp <- nrow(df)

  if (n_comp > top_n) {
    top_df  <- df[seq_len(top_n), , drop = FALSE]
    rest_n  <- sum(df$n_nodes[(top_n + 1L):n_comp])
    other   <- data.frame(component_id = top_n + 1L,
                          n_nodes      = rest_n,
                          n_edges      = NA_integer_,
                          stringsAsFactors = FALSE)
    plot_df <- rbind(top_df, other)
    x_label <- c(as.character(seq_len(top_n)), sprintf("rest\n(%d)", n_comp - top_n))
  } else {
    plot_df <- df
    x_label <- as.character(df$component_id)
  }

  plot_df$x_label <- factor(x_label, levels = x_label)

  # Show labels at position 1, multiples of 5, and the final "rest" bar.
  # Position 20 is suppressed when a "rest" bar is present to avoid overlap.
  n_bars         <- nrow(plot_df)
  has_rest       <- n_comp > top_n
  display_labels <- ifelse(
    seq_len(n_bars) == 1L |
      (seq_len(n_bars) %% 5 == 0 & !(has_rest & seq_len(n_bars) == top_n)) |
      seq_len(n_bars) == n_bars,
    as.character(plot_df$x_label),
    ""
  )

  ggplot2::ggplot(plot_df,
                  ggplot2::aes(x = x_label, y = n_nodes)) +
    ggplot2::geom_col(fill = .ISOFORM_COLOUR, colour = "white", linewidth = 0.2) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_discrete(labels = display_labels) +
    ggplot2::labs(
      x = "Component rank",
      y = "Nodes (log scale)"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7)) +
    .diag_theme()
}


#' Panel E — Isoforms per gene histogram (gene-level, bottom row)
#' @noRd
.panel_isoforms_per_gene <- function(diag) {
  df <- .df(diag@isoform_per_gene_table)

  counts <- as.data.frame(table(n_isoforms_total = df$n_isoforms_total),
                          stringsAsFactors = FALSE)
  names(counts) <- c("n_isoforms_total", "n")
  counts$n_isoforms_total <- as.numeric(as.character(counts$n_isoforms_total))

  ggplot2::ggplot(counts, ggplot2::aes(x = n_isoforms_total, y = n)) +
    ggplot2::geom_col(
      fill      = .GENE_COLOUR,
      colour    = "white",
      linewidth = 0.2,
      width     = 0.8
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Number of isoforms",
      y = "Number of genes (log scale)"
    ) +
    .diag_theme()
}


#' Panel C — Clustering coefficient vs degree scatter (isoform-level, top row)
#' @noRd
.panel_clustering <- function(diag) {
  df <- .df(diag@clustering_table)
  df <- df[df$degree > 0L, , drop = FALSE]

  # Downsample for large networks to keep rendering fast.
  # Uses a fixed local seed so the figure is reproducible across calls;
  # the caller's global RNG state is saved and restored afterwards.
  max_pts <- 5000L
  if (nrow(df) > max_pts) {
    keep <- withr::with_seed(42L, sample.int(nrow(df), max_pts))
    df <- df[keep, , drop = FALSE]
  }

  ggplot2::ggplot(df,
                  ggplot2::aes(x = degree,
                               y = clustering_coefficient)) +
    ggplot2::geom_point(
      alpha  = 0.25,
      size   = 0.8,
      colour = .ISOFORM_COLOUR
    ) +
    ggplot2::geom_smooth(
      method    = "loess",
      formula   = y ~ x,
      se        = TRUE,
      colour    = "grey10",
      fill      = "grey70",
      linewidth = 0.7
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Degree (log scale)",
      y = "Local clustering coefficient"
    ) +
    .diag_theme()
}


#' Panel F — Gene connectivity status bar chart (categorical, bottom row)
#' @noRd
.panel_connectivity_status <- function(diag) {
  df <- .df(diag@isoform_per_gene_table)

  counts <- as.data.frame(table(df$connectivity_status),
                          stringsAsFactors = FALSE)
  names(counts) <- c("status", "n_genes")

  # Ensure all three categories are present even if count is zero
  all_statuses <- c("all_connected", "partial", "all_isolated")
  missing      <- setdiff(all_statuses, counts$status)
  if (length(missing) > 0L) {
    counts <- rbind(counts,
                    data.frame(status  = missing,
                               n_genes = 0L,
                               stringsAsFactors = FALSE))
  }

  counts$status <- factor(counts$status,
                          levels = c("all_connected", "partial", "all_isolated"),
                          labels = c("All connected", "Partial", "All isolated"))

  ggplot2::ggplot(counts,
                  ggplot2::aes(x = status, y = n_genes, fill = status)) +
    ggplot2::geom_col(colour = "white", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_genes),
      vjust  = -0.4,
      size   = 3,
      colour = "grey20"
    ) +
    # Headroom so the count above the tallest bar is not clipped
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.1))
    ) +
    ggplot2::scale_fill_manual(
      values = c("All connected" = "#009E73",
                 "Partial"        = "#D55E00",
                 "All isolated"   = "#CC79A7"),
      guide  = "none"
    ) +
    ggplot2::labs(
      x = "Status",
      y = "Number of genes"
    ) +
    .diag_theme()
}


# ==============================================================================
#  Exported function
# ==============================================================================

#' Plot network diagnostics
#'
#' Produces a 2-row x 3-column panel figure summarising the
#' \linkS4class{NetworkDiagnostics} result stored in an
#' \linkS4class{IsoNetAnalyzeR} session object.  Panels carry no in-plot
#' titles; panel content and the overall figure title are intended to live
#' in the surrounding caption/body text rather than the image itself.
#'
#' \strong{Panel layout (2 rows x 3 columns), grouped by network level:}
#' \itemize{
#'   \item Row 1 (isoform-level, blue): (A) Isoform degree distribution |
#'     (B) Connected component sizes | (C) Clustering coefficient vs degree
#'   \item Row 2 (gene-level, orange / categorical): (D) Gene degree
#'     distribution | (E) Isoforms per gene | (F) Gene connectivity status
#' }
#'
#' Sticky degree thresholds are shown as dashed vertical lines in the two
#' degree distribution panels (A, D).  A LOESS smoothing line is overlaid on
#' the clustering coefficient panel (C).  Degree-distribution panels (A, D)
#' and the isoforms-per-gene panel (E) use log-scaled counts to make the
#' heavy right tail readable; no power-law shape is implied by this choice.
#'
#' @param obj    An \linkS4class{IsoNetAnalyzeR} session object with a
#'   completed network diagnostics run (see \code{\link{run_network_diagnostics}}).
#' @param title  Character scalar or \code{NULL}.  Optional overall figure
#'   title, off by default since titles are intended to live in the
#'   caption rather than the figure itself.  Set to a string to override.
#'
#' @return A \code{patchwork} object (invisibly).  Print it to display, or
#'   save with \code{ggplot2::ggsave()}.
#'
#' @examples
#' # A small deterministic network keeps this example fast; on the bundled
#' # network the diagnostics take about half a minute.
#' g <- igraph::make_graph("Zachary")
#' n <- igraph::vcount(g)
#' igraph::V(g)$name <- paste0("tx", seq_len(n))
#' igraph::V(g)$gene <- paste0("gene", ceiling(seq_len(n) / 2))
#' igraph::V(g)$transcript_biotype <- "protein_coding"
#' igraph::E(g)$interact <- 1L
#' toy <- IsoNetAnalyzeR("toy", network = g,
#'                       metadata = list(network_cutoff = 0.25))
#' toy <- run_network_diagnostics(toy)
#' p <- plot_network_diagnostics(toy)
#'
#' @export
plot_network_diagnostics <- function(obj, title = NULL) {

  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (is.null(obj@network_diagnostics))
    stop(
      "No network diagnostics found. Run `run_network_diagnostics(obj)` first.",
      call. = FALSE
    )

  diag <- obj@network_diagnostics

  p_iso_degree   <- .panel_isoform_degree(diag)        # A
  p_components   <- .panel_components(diag)            # B
  p_clustering   <- .panel_clustering(diag)             # C
  p_gene_degree  <- .panel_gene_degree(diag)            # D
  p_iso_per_gene <- .panel_isoforms_per_gene(diag)       # E
  p_connectivity <- .panel_connectivity_status(diag)     # F

  fig <- (p_iso_degree | p_components | p_clustering) /
    (p_gene_degree | p_iso_per_gene | p_connectivity)

  ann_args <- list(tag_levels = "A")
  if (!is.null(title)) {
    ann_args$title <- title
    ann_args$theme <- ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12)
    )
  }
  fig <- fig + do.call(patchwork::plot_annotation, ann_args)

  fig
}
