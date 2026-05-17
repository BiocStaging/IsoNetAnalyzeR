## plot_network_diagnostics.R
## Plotting function for NetworkDiagnostics results.
##
## Produces a 2-row x 3-column panel figure suitable for thesis and
## publication use.  All panels use a minimal, greyscale-friendly aesthetic.
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
  "n_isoforms_total"
))


# ==============================================================================
#  Internal panel helpers
# ==============================================================================

#' Shared minimal ggplot2 theme for all diagnostic panels
#' @noRd
.diag_theme <- function() {
  ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 10),
      axis.title       = ggplot2::element_text(size = 9),
      axis.text        = ggplot2::element_text(size = 8),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      panel.grid.minor = ggplot2::element_blank()
    )
}


#' Panel 1 — Isoform degree distribution
#' @noRd
.panel_isoform_degree <- function(diag) {
  df <- diag@isoform_degree_table
  df <- df[df$degree > 0L, , drop = FALSE]   # exclude isolated

  ggplot2::ggplot(df, ggplot2::aes(x = degree)) +
    ggplot2::geom_histogram(
      bins   = 40,
      fill   = "grey40",
      colour = "white",
      linewidth = 0.2
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(
      xintercept = diag@sticky_isoform_threshold,
      linetype   = "dashed",
      colour     = "grey20",
      linewidth  = 0.6
    ) +
    ggplot2::annotate(
      "text",
      x      = diag@sticky_isoform_threshold,
      y      = Inf,
      hjust  = -0.15,
      vjust  = 1.4,
      size   = 2.8,
      colour = "grey20",
      label  = "sticky\nthreshold"
    ) +
    ggplot2::labs(
      title = "Isoform Degree Distribution",
      x     = "Degree (log scale)",
      y     = "Count"
    ) +
    .diag_theme()
}


#' Panel 2 — Gene degree distribution
#' @noRd
.panel_gene_degree <- function(diag) {
  df <- diag@gene_degree_table

  ggplot2::ggplot(df, ggplot2::aes(x = degree)) +
    ggplot2::geom_histogram(
      bins      = 40,
      fill      = "grey40",
      colour    = "white",
      linewidth = 0.2
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::geom_vline(
      xintercept = diag@sticky_gene_threshold,
      linetype   = "dashed",
      colour     = "grey20",
      linewidth  = 0.6
    ) +
    ggplot2::annotate(
      "text",
      x      = diag@sticky_gene_threshold,
      y      = Inf,
      hjust  = -0.15,
      vjust  = 1.4,
      size   = 2.8,
      colour = "grey20",
      label  = "sticky\nthreshold"
    ) +
    ggplot2::labs(
      title = "Gene Degree Distribution",
      x     = "Degree (log scale)",
      y     = "Count"
    ) +
    .diag_theme()
}


#' Panel 3 — Connected component sizes
#' @noRd
.panel_components <- function(diag) {
  df <- diag@component_table

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
    ggplot2::geom_col(fill = "grey40", colour = "white", linewidth = 0.2) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_discrete(labels = display_labels) +
    ggplot2::labs(
      title = "Connected Component Sizes",
      x     = "Component rank",
      y     = "Nodes (log scale)"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7)) +
    .diag_theme()
}


#' Panel 4 — Isoforms per gene histogram
#' @noRd
.panel_isoforms_per_gene <- function(diag) {
  df <- diag@isoform_per_gene_table

  ggplot2::ggplot(df, ggplot2::aes(x = n_isoforms_total)) +
    ggplot2::geom_histogram(
      binwidth  = 1,
      fill      = "grey40",
      colour    = "white",
      linewidth = 0.2
    ) +
    ggplot2::labs(
      title = "Isoforms per Gene",
      x     = "Number of isoforms",
      y     = "Number of genes"
    ) +
    .diag_theme()
}


#' Panel 5 — Clustering coefficient vs degree scatter
#' @noRd
.panel_clustering <- function(diag) {
  df <- diag@clustering_table
  df <- df[df$degree > 0L, , drop = FALSE]

  # Downsample for large networks to keep rendering fast
  max_pts <- 5000L
  if (nrow(df) > max_pts)
    df <- df[sample.int(nrow(df), max_pts), , drop = FALSE]

  ggplot2::ggplot(df,
                  ggplot2::aes(x = degree,
                               y = clustering_coefficient)) +
    ggplot2::geom_point(
      alpha  = 0.25,
      size   = 0.8,
      colour = "grey30"
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
      title = "Clustering Coefficient vs Degree",
      x     = "Degree (log scale)",
      y     = "Local clustering coefficient"
    ) +
    .diag_theme()
}


#' Panel 6 — Gene connectivity status bar chart
#' @noRd
.panel_connectivity_status <- function(diag) {
  df <- diag@isoform_per_gene_table

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
                  ggplot2::aes(x = status, y = n_genes)) +
    ggplot2::geom_col(fill = "grey40", colour = "white", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_genes),
      vjust  = -0.4,
      size   = 3,
      colour = "grey20"
    ) +
    ggplot2::labs(
      title = "Gene Connectivity Status",
      x     = "Status",
      y     = "Number of genes"
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
#' \linkS4class{IsoNetAnalyzeR} session object.
#'
#' \strong{Panel layout (2 rows x 3 columns):}
#' \itemize{
#'   \item Row 1: Isoform degree distribution | Gene degree distribution | Connected component sizes
#'   \item Row 2: Isoforms per gene | Clustering coefficient vs degree | Gene connectivity status
#' }
#'
#' Sticky degree thresholds are shown as dashed vertical lines in the two
#' degree distribution panels.  A LOESS smoothing line is overlaid on the
#' clustering coefficient panel.
#'
#' @param obj    An \linkS4class{IsoNetAnalyzeR} session object with a
#'   completed network diagnostics run (see \code{\link{run_network_diagnostics}}).
#' @param title  Character scalar or \code{NULL}.  Overall figure title.
#'   Defaults to \code{"Network Diagnostics"}.
#'
#' @return A \code{patchwork} object (invisibly).  Print it to display, or
#'   save with \code{ggplot2::ggsave()}.
#'
#' @examples
#' \dontrun{
#' obj <- run_network_diagnostics(obj)
#'
#' # Display
#' plot_network_diagnostics(obj)
#'
#' # Save to file
#' p <- plot_network_diagnostics(obj)
#' ggplot2::ggsave("network_diagnostics.pdf", p, width = 12, height = 7)
#' }
#'
#' @export
plot_network_diagnostics <- function(obj, title = "Network Diagnostics") {

  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (is.null(obj@network_diagnostics))
    stop(
      "No network diagnostics found. Run `run_network_diagnostics(obj)` first.",
      call. = FALSE
    )

  if (!requireNamespace("ggplot2",   quietly = TRUE))
    stop("Package 'ggplot2' is required. Install with install.packages('ggplot2').",
         call. = FALSE)

  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required. Install with install.packages('patchwork').",
         call. = FALSE)

  diag <- obj@network_diagnostics

  p1 <- .panel_isoform_degree(diag)
  p2 <- .panel_gene_degree(diag)
  p3 <- .panel_components(diag)
  p4 <- .panel_isoforms_per_gene(diag)
  p5 <- .panel_clustering(diag)
  p6 <- .panel_connectivity_status(diag)

  fig <- (p1 | p2 | p3) / (p4 | p5 | p6)

  if (!is.null(title)) {
    fig <- fig +
      patchwork::plot_annotation(
        title = title,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 12)
        )
      )
  }

  print(fig)
  invisible(fig)
}
