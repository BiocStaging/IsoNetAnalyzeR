## network_diagnostics_class.R
## S4 class definition for a network diagnostics result.
##
## Produced by run_network_diagnostics(). Stores structural and topological
## summary statistics of the loaded isoform interaction network, degree
## distributions at isoform and gene level, connected component structure,
## sticky node tables, and gene-level isolation statistics.
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: NetworkDiagnostics class, show method

#' @include isonetanalyzer_class.R
NULL


# ==============================================================================
#  Class definition
# ==============================================================================

#' NetworkDiagnostics - network characterisation result
#'
#' Top-level container produced by \code{\link{run_network_diagnostics}}.
#' Stores structural and topological summary statistics of the loaded isoform
#' interaction network.  Stored in the \code{@@network_diagnostics} slot of an
#' \linkS4class{IsoNetAnalyzeR} session object; retrieved with
#' \code{\link{get_network_diagnostics}}.
#'
#' @slot network_stats     Data frame with one row summarising the three
#'   network views and global topology metrics.  Columns:
#'   \describe{
#'     \item{full_isoforms}{Total isoform nodes in the full network.}
#'     \item{full_genes}{Unique genes in the full network.}
#'     \item{full_edges}{Total edges in the session network. For the
#'       bundled network (from 0.99.0) this equals the number of predicted
#'       interaction edges.}
#'     \item{full_predicted_edges}{Edges with \code{interact == 1}.}
#'     \item{pred_isoforms}{Isoform nodes in the predicted subnetwork.}
#'     \item{pred_genes}{Unique genes in the predicted subnetwork.}
#'     \item{pred_edges}{Edges in the predicted subnetwork.}
#'     \item{pred_density}{Edge density of the predicted isoform subnetwork.}
#'     \item{pred_avg_degree}{Mean isoform degree in the predicted subnetwork.}
#'     \item{pred_avg_clustering}{Mean local clustering coefficient.}
#'     \item{pred_n_components}{Number of connected components.}
#'     \item{pred_largest_component_n}{Node count of the largest component.}
#'     \item{pred_largest_component_pct}{Percentage of nodes in the largest
#'       component.}
#'     \item{pred_diameter}{Diameter of the largest component.}
#'     \item{pred_avg_path_length}{Mean shortest path length in the largest
#'       component.}
#'     \item{gene_nodes}{Nodes in the gene-level predicted network.}
#'     \item{gene_edges}{Edges in the gene-level predicted network.}
#'     \item{isolated_isoforms}{Isoforms with zero predicted edges.}
#'     \item{isolated_isoforms_pct}{Percentage of all isoforms that are
#'       isolated.}
#'   }
#'
#' @slot isoform_degree_table  Data frame with one row per isoform in the
#'   predicted subnetwork.  Columns: \code{isoform_id}, \code{gene_id},
#'   \code{degree}.  Isoforms with zero predicted edges are included with
#'   \code{degree = 0}.
#'
#' @slot gene_degree_table     Data frame with one row per gene in the
#'   gene-level predicted network.  Columns: \code{gene_id}, \code{degree}.
#'
#' @slot component_table       Data frame with one row per connected component
#'   in the predicted isoform subnetwork, sorted descending by size.  Columns:
#'   \code{component_id}, \code{n_nodes}, \code{n_edges}.
#'
#' @slot clustering_table      Data frame with one row per isoform in the
#'   predicted subnetwork.  Columns: \code{isoform_id}, \code{gene_id},
#'   \code{degree}, \code{clustering_coefficient}.
#'
#' @slot isoform_per_gene_table  Data frame with one row per gene present in
#'   the full network.  Columns: \code{gene_id}, \code{n_isoforms_total}
#'   (all isoforms for this gene in the full network),
#'   \code{n_isoforms_connected} (isoforms with at least one predicted edge),
#'   \code{n_isoforms_isolated} (isoforms with zero predicted edges),
#'   \code{connectivity_status} (character: \code{"all_connected"},
#'   \code{"partial"}, or \code{"all_isolated"}).
#'
#' @slot sticky_isoform_table  Data frame with one row per sticky isoform
#'   (degree above threshold).  Columns: \code{isoform_id}, \code{gene_id},
#'   \code{degree}.  Empty data frame if no sticky isoforms found.
#'
#' @slot sticky_gene_table     Data frame with one row per sticky gene.
#'   Columns: \code{gene_id}, \code{degree}.  Empty data frame if none found.
#'
#' @slot sticky_method         Character scalar.  Threshold method used:
#'   \code{"quantile"} or \code{"sd"}.
#'
#' @slot sticky_cutoff         Numeric scalar.  Cutoff value supplied to the
#'   threshold method.
#'
#' @slot sticky_isoform_threshold  Numeric scalar.  The computed degree
#'   threshold above which an isoform is classified as sticky.
#'
#' @slot sticky_gene_threshold     Numeric scalar.  The computed degree
#'   threshold above which a gene is classified as sticky.
#'
#' @slot created_at            \code{POSIXct} scalar.  Timestamp of the run.
#'
#' @return Objects of this class are created by \code{\link{run_network_diagnostics}} and returned by \code{\link{get_network_diagnostics}}.
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
#' get_network_diagnostics(toy)
#'
#' @exportClass NetworkDiagnostics
setClass(
  "NetworkDiagnostics",
  representation(
    network_stats            = "data.frame",
    isoform_degree_table     = "data.frame",
    gene_degree_table        = "data.frame",
    component_table          = "data.frame",
    clustering_table         = "data.frame",
    isoform_per_gene_table   = "data.frame",
    sticky_isoform_table     = "data.frame",
    sticky_gene_table        = "data.frame",
    sticky_method            = "character",
    sticky_cutoff            = "numeric",
    sticky_isoform_threshold = "numeric",
    sticky_gene_threshold    = "numeric",
    created_at               = "POSIXct"
  )
)


# ==============================================================================
#  Validator
# ==============================================================================

setValidity("NetworkDiagnostics", function(object) {
  msgs <- character()

  # network_stats must have exactly one row
  if (nrow(object@network_stats) != 1L)
    msgs <- c(msgs, "@network_stats must be a single-row data frame.")

  # sticky_method
  if (!object@sticky_method %in% c("quantile", "sd"))
    msgs <- c(msgs, "@sticky_method must be 'quantile' or 'sd'.")

  # sticky_cutoff
  if (length(object@sticky_cutoff) != 1L || is.na(object@sticky_cutoff))
    msgs <- c(msgs, "@sticky_cutoff must be a single non-NA numeric scalar.")

  # created_at
  if (length(object@created_at) != 1L)
    msgs <- c(msgs, "@created_at must be a length-1 POSIXct.")

  # Required columns for each table
  req <- list(
    isoform_degree_table   = c("isoform_id", "gene_id", "degree"),
    gene_degree_table      = c("gene_id", "degree"),
    component_table        = c("component_id", "n_nodes", "n_edges"),
    clustering_table       = c("isoform_id", "gene_id", "degree",
                               "clustering_coefficient"),
    isoform_per_gene_table = c("gene_id", "n_isoforms_total",
                               "n_isoforms_connected", "n_isoforms_isolated",
                               "connectivity_status"),
    sticky_isoform_table   = c("isoform_id", "gene_id", "degree"),
    sticky_gene_table      = c("gene_id", "degree")
  )

  for (slot_name in names(req)) {
    tbl     <- slot(object, slot_name)
    missing <- setdiff(req[[slot_name]], names(tbl))
    if (length(missing))
      msgs <- c(msgs, sprintf(
        "@%s is missing column(s): %s",
        slot_name, paste(missing, collapse = ", ")
      ))
  }

  if (length(msgs)) msgs else TRUE
})


# ==============================================================================
#  show()
# ==============================================================================

#' @describeIn NetworkDiagnostics Print a concise diagnostics summary.
#' @param object A \linkS4class{NetworkDiagnostics} object.
setMethod("show", "NetworkDiagnostics", function(object) {
  s <- object@network_stats

  cli::cli_rule("NetworkDiagnostics")
  .cli_aligned("Created    : {format(object@created_at, '%Y-%m-%d %H:%M:%S')}")
  .cli_aligned("")

  .cli_aligned("{.strong Full Network}")
  .cli_aligned("  Isoforms / genes     : {s$full_isoforms} / {s$full_genes}")
  .cli_aligned("  Total edges          : {s$full_edges}")
  .cli_aligned("  Predicted edges      : {s$full_predicted_edges}")
  .cli_aligned("  Isolated isoforms    : {s$isolated_isoforms} ({round(s$isolated_isoforms_pct, 1)}%)")
  .cli_aligned("")

  .cli_aligned("{.strong Predicted Isoform Network}")
  .cli_aligned("  Isoforms / genes     : {s$pred_isoforms} / {s$pred_genes}")
  .cli_aligned("  Edges                : {s$pred_edges}")
  .cli_aligned("  Density              : {round(s$pred_density, 4)}")
  .cli_aligned("  Avg degree           : {round(s$pred_avg_degree, 2)}")
  .cli_aligned("  Avg clustering       : {round(s$pred_avg_clustering, 4)}")
  .cli_aligned("  Components           : {s$pred_n_components}")
  .cli_aligned("  Largest component    : {s$pred_largest_component_n} nodes ({round(s$pred_largest_component_pct, 1)}%)")
  .cli_aligned("  Diameter             : {s$pred_diameter}")
  .cli_aligned("  Avg path length      : {round(s$pred_avg_path_length, 2)}")
  .cli_aligned("")

  .cli_aligned("{.strong Predicted Gene-Level Network}")
  .cli_aligned("  Genes / edges        : {s$gene_nodes} / {s$gene_edges}")
  .cli_aligned("")

  .cli_aligned("{.strong Sticky Nodes}  [{object@sticky_method} / cutoff = {object@sticky_cutoff}]")
  .cli_aligned("  Sticky isoforms      : {nrow(object@sticky_isoform_table)} (threshold degree > {round(object@sticky_isoform_threshold, 2)})")
  .cli_aligned("  Sticky genes         : {nrow(object@sticky_gene_table)} (threshold degree > {round(object@sticky_gene_threshold, 2)})")
  .cli_aligned("")

  .cli_aligned("{.strong Gene Connectivity Status}")
  status_counts <- table(object@isoform_per_gene_table$connectivity_status)
  # Same wording as panel F of plot_network_diagnostics()
  status_labels <- c(all_connected = "All connected",
                     partial       = "Partial",
                     all_isolated  = "All isolated")
  for (s_name in names(status_labels)) {
    n <- if (s_name %in% names(status_counts)) status_counts[[s_name]] else 0L
    label <- formatC(status_labels[[s_name]], width = -21L)
    .cli_aligned("  ", label, ": {n} gene{?s}")
  }

  cli::cli_rule()
  invisible(object)
})
