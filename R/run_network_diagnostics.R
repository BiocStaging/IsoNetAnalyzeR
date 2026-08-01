## run_network_diagnostics.R
## User-facing function to run a comprehensive network diagnostics analysis.
##
## Computes structural and topological summary statistics of the loaded isoform
## interaction network, degree distributions, connected component structure,
## sticky node identification, and gene-level isolation statistics.
##
## Results are stored in the @network_diagnostics slot of the IsoNetAnalyzeR
## session object.  Re-running overwrites the previous result.
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: run_network_diagnostics

#' @include isonetanalyzer_class.R
NULL
#' @include network_diagnostics_class.R
NULL


# ==============================================================================
#  Internal helpers
# ==============================================================================

#' Compute the degree threshold for sticky node detection
#'
#' @param deg    Named numeric vector of degrees.
#' @param method Character scalar. \code{"quantile"} or \code{"sd"}.
#' @param cutoff Numeric scalar. Quantile probability or SD multiplier.
#' @return Numeric scalar threshold.
#' @noRd
.sticky_threshold <- function(deg, method, cutoff) {
  if (length(deg) == 0L) return(NA_real_)
  switch(method,
         quantile = unname(stats::quantile(deg, cutoff)),
         sd       = mean(deg) + cutoff * stats::sd(deg),
         stop(sprintf("`sticky_method` must be 'quantile' or 'sd', got: '%s'.",
                      method), call. = FALSE)
  )
}


#' Compute average shortest path length of the largest component
#'
#' Falls back to NA if the component has only one node.
#'
#' @param g igraph object (should be the largest component subgraph).
#' @return Numeric scalar or NA.
#' @noRd
.avg_path_length <- function(g) {
  if (igraph::vcount(g) < 2L) return(NA_real_)
  igraph::mean_distance(g, directed = FALSE, unconnected = FALSE)
}


#' Compute diameter of the largest component
#'
#' Falls back to NA if the component has only one node.
#'
#' @param g igraph object (should be the largest component subgraph).
#' @return Integer scalar or NA.
#' @noRd
.diameter <- function(g) {
  if (igraph::vcount(g) < 2L) return(NA_integer_)
  igraph::diameter(g, directed = FALSE)
}


# ==============================================================================
#  Exported function
# ==============================================================================

#' Run a comprehensive network diagnostics analysis
#'
#' Characterises the isoform interaction network loaded in an
#' \linkS4class{IsoNetAnalyzeR} session object.  Computes:
#' \itemize{
#'   \item Structural summary statistics across three network views (full,
#'     predicted isoform, gene-level predicted).
#'   \item Global topology metrics: density, average degree, clustering
#'     coefficient, connected components, diameter, and average path length
#'     of the largest component.
#'   \item Per-isoform and per-gene degree tables.
#'   \item Connected component size table.
#'   \item Clustering coefficient vs degree table.
#'   \item Per-gene isoform connectivity status
#'     (\code{"all_connected"}, \code{"partial"}, \code{"all_isolated"}).
#'   \item Sticky isoform and gene tables based on a degree threshold.
#' }
#'
#' Results are stored in \code{obj@@network_diagnostics} as a
#' \linkS4class{NetworkDiagnostics} object and can be retrieved with
#' \code{\link{get_network_diagnostics}}.  Re-running overwrites any previous
#' result with a warning.
#'
#' @section Sticky node detection:
#' An isoform (or gene) is classified as sticky if its degree in the predicted
#' subnetwork exceeds a threshold derived from the overall degree distribution.
#' Two methods are supported:
#' \describe{
#'   \item{\code{"quantile"} (default)}{Threshold is the \code{sticky_cutoff}
#'     quantile of the degree distribution.  \code{sticky_cutoff = 0.95} means
#'     the top 5\% by degree are sticky.  Robust to the heavy-tailed
#'     distributions typical of protein interaction networks.}
#'   \item{\code{"sd"}}{Threshold is \code{mean + sticky_cutoff * SD}.
#'     More sensitive to extreme hubs inflating the mean.}
#' }
#'
#' @param obj           An \linkS4class{IsoNetAnalyzeR} session object.
#' @param sticky_method Character scalar. Threshold method for sticky node
#'   detection: \code{"quantile"} (default) or \code{"sd"}.
#' @param sticky_cutoff Numeric scalar. For \code{"quantile"}: a probability
#'   in \[0, 1\] (default \code{0.95}).  For \code{"sd"}: a positive multiplier
#'   (default \code{2}).
#'
#' @return The updated \linkS4class{IsoNetAnalyzeR} object, invisibly.
#'   Access the result with \code{get_network_diagnostics(obj)}.
#'
#' @examples
#' \dontrun{
#' obj <- IsoNetAnalyzeR("demo")
#'
#' # Default: quantile method, top 5%
#' obj <- run_network_diagnostics(obj)
#' get_network_diagnostics(obj)
#'
#' # SD method
#' obj <- run_network_diagnostics(obj, sticky_method = "sd", sticky_cutoff = 2)
#'
#' # Stricter quantile threshold
#' obj <- run_network_diagnostics(obj, sticky_cutoff = 0.99)
#' }
#'
#' @export
run_network_diagnostics <- function(obj,
                                    sticky_method = c("quantile", "sd"),
                                    sticky_cutoff = NULL) {

  # -- Input validation --------------------------------------------------------
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  sticky_method <- match.arg(sticky_method)

  if (is.null(sticky_cutoff)) {
    sticky_cutoff <- if (sticky_method == "quantile") 0.95 else 2
  }

  if (!is.numeric(sticky_cutoff) || length(sticky_cutoff) != 1L ||
      is.na(sticky_cutoff) || !is.finite(sticky_cutoff) || sticky_cutoff <= 0)
    stop("`sticky_cutoff` must be a single positive numeric scalar.",
         call. = FALSE)

  if (sticky_method == "quantile" && sticky_cutoff >= 1)
    stop("`sticky_cutoff` must be in (0, 1) when `sticky_method = 'quantile'`.",
         call. = FALSE)

  if (!is.null(obj@network_diagnostics))
    warning(
      "Overwriting existing network diagnostics result.",
      call. = FALSE
    )

  g <- obj@network

  .validate_network_contract(g)

  message("Running network diagnostics...")

  # -- Build subgraphs ---------------------------------------------------------
  message("  Building predicted subgraphs...")

  pred_edges <- which(igraph::edge_attr(g, "interact") == 1L)
  has_predicted_edges <- length(pred_edges) > 0L

  g_pred <- if (has_predicted_edges) {
    igraph::subgraph_from_edges(g, pred_edges, delete.vertices = TRUE)
  } else {
    igraph::make_empty_graph(n = 0L, directed = FALSE)
  }

  g_gene <- if (has_predicted_edges) {
    edges_df    <- igraph::as_data_frame(g_pred, what = "edges")
    iso_to_gene <- stats::setNames(igraph::V(g_pred)$gene,
                                   igraph::V(g_pred)$name)
    edges_df$from_gene <- iso_to_gene[edges_df$from]
    edges_df$to_gene   <- iso_to_gene[edges_df$to]
    igraph::simplify(
      igraph::graph_from_data_frame(
        edges_df[, c("from_gene", "to_gene")], directed = FALSE
      ),
      remove.multiple = TRUE, remove.loops = TRUE
    )
  } else {
    igraph::make_empty_graph(n = 0L, directed = FALSE)
  }

  # -- Connected components ----------------------------------------------------
  message("  Computing connected components...")

  if (has_predicted_edges) {
    comps      <- igraph::components(g_pred)
    comp_sizes <- sort(comps$csize, decreasing = TRUE)
    largest_id <- which.max(comps$csize)
    g_largest  <- igraph::induced_subgraph(
      g_pred, igraph::V(g_pred)[comps$membership == largest_id]
    )

    component_table <- data.frame(
      component_id = as.integer(seq_along(comp_sizes)),
      n_nodes      = as.integer(comp_sizes),
      n_edges      = as.integer(vapply(seq_along(comp_sizes), function(i) {
        vid <- igraph::V(g_pred)[comps$membership ==
                                   order(comps$csize, decreasing = TRUE)[i]]
        igraph::ecount(igraph::induced_subgraph(g_pred, vid))
      }, numeric(1L))),
      stringsAsFactors = FALSE
    )
  } else {
    comps <- list(no = 0L)
    g_largest <- igraph::make_empty_graph(n = 0L, directed = FALSE)
    component_table <- data.frame(
      component_id = integer(0), n_nodes = integer(0), n_edges = integer(0)
    )
  }

  # -- Degree tables -----------------------------------------------------------
  message("  Computing degree distributions...")

  # All isoforms in the full network, degree 0 for those not in g_pred
  all_isoforms <- igraph::V(g)$name
  all_genes    <- igraph::V(g)$gene
  pred_deg     <- igraph::degree(g_pred)

  isoform_degree_table <- data.frame(
    isoform_id = all_isoforms,
    gene_id    = all_genes,
    degree     = as.integer(pred_deg[all_isoforms]),
    stringsAsFactors = FALSE
  )
  isoform_degree_table$degree[is.na(isoform_degree_table$degree)] <- 0L

  gene_deg <- if (igraph::vcount(g_gene)) igraph::degree(g_gene) else numeric(0)
  gene_degree_table <- data.frame(
    gene_id = if (length(gene_deg)) names(gene_deg) else character(0),
    degree  = as.integer(gene_deg),
    stringsAsFactors = FALSE
  )

  # -- Clustering table --------------------------------------------------------
  message("  Computing clustering coefficients...")

  clust_vals  <- if (has_predicted_edges)
    igraph::transitivity(g_pred, type = "local", isolates = "zero") else numeric(0)
  pred_names  <- if (has_predicted_edges) igraph::V(g_pred)$name else character(0)
  pred_genes  <- if (has_predicted_edges) igraph::V(g_pred)$gene else character(0)
  pred_deg_v  <- if (has_predicted_edges) igraph::degree(g_pred) else numeric(0)

  clustering_table <- data.frame(
    isoform_id              = pred_names,
    gene_id                 = pred_genes,
    degree                  = as.integer(pred_deg_v),
    clustering_coefficient  = clust_vals,
    stringsAsFactors        = FALSE
  )

  # -- Isoform per gene table --------------------------------------------------
  message("  Computing per-gene connectivity status...")

  gene_ids_full <- unique(igraph::V(g)$gene)

  isoform_per_gene_table <- do.call(rbind, lapply(gene_ids_full, function(gid) {
    isos       <- igraph::V(g)$name[igraph::V(g)$gene == gid]
    n_total    <- length(isos)
    degrees    <- isoform_degree_table$degree[
      isoform_degree_table$isoform_id %in% isos
    ]
    n_conn     <- sum(degrees > 0L)
    n_iso_     <- sum(degrees == 0L)
    status     <- if (n_conn == 0L) "all_isolated" else
      if (n_iso_ == 0L) "all_connected" else "partial"
    data.frame(
      gene_id              = gid,
      n_isoforms_total     = n_total,
      n_isoforms_connected = n_conn,
      n_isoforms_isolated  = n_iso_,
      connectivity_status  = status,
      stringsAsFactors     = FALSE
    )
  }))
  rownames(isoform_per_gene_table) <- NULL

  # -- Network summary stats ---------------------------------------------------
  message("  Computing global topology metrics...")

  full_isoforms        <- igraph::vcount(g)
  full_genes           <- length(unique(igraph::V(g)$gene))
  full_edges           <- igraph::ecount(g)
  full_predicted_edges <- length(pred_edges)
  isolated_isoforms    <- sum(isoform_degree_table$degree == 0L)

  network_stats <- data.frame(
    # Full network
    full_isoforms            = full_isoforms,
    full_genes               = full_genes,
    full_edges               = full_edges,
    full_predicted_edges     = full_predicted_edges,
    # Predicted isoform network
    pred_isoforms            = igraph::vcount(g_pred),
    pred_genes               = length(unique(igraph::V(g_pred)$gene)),
    pred_edges               = igraph::ecount(g_pred),
    pred_density             = if (has_predicted_edges) igraph::edge_density(g_pred) else NA_real_,
    pred_avg_degree          = if (has_predicted_edges) mean(igraph::degree(g_pred)) else NA_real_,
    pred_avg_clustering      = if (has_predicted_edges) igraph::transitivity(g_pred, type = "average") else NA_real_,
    pred_n_components        = comps$no,
    pred_largest_component_n = igraph::vcount(g_largest),
    pred_largest_component_pct = if (has_predicted_edges)
      100 * igraph::vcount(g_largest) / igraph::vcount(g_pred) else NA_real_,
    pred_diameter            = .diameter(g_largest),
    pred_avg_path_length     = .avg_path_length(g_largest),
    # Gene-level network
    gene_nodes               = igraph::vcount(g_gene),
    gene_edges               = igraph::ecount(g_gene),
    # Isolation
    isolated_isoforms        = isolated_isoforms,
    isolated_isoforms_pct    = 100 * isolated_isoforms / full_isoforms,
    stringsAsFactors         = FALSE
  )

  # -- Sticky nodes ------------------------------------------------------------
  message("  Identifying sticky nodes...")

  iso_deg_nonzero <- isoform_degree_table[isoform_degree_table$degree > 0L, ]
  iso_threshold   <- .sticky_threshold(iso_deg_nonzero$degree,
                                       sticky_method, sticky_cutoff)

  sticky_isoform_table <- iso_deg_nonzero[
    iso_deg_nonzero$degree > iso_threshold,
    c("isoform_id", "gene_id", "degree"),
    drop = FALSE
  ]
  sticky_isoform_table <- sticky_isoform_table[
    order(sticky_isoform_table$degree, decreasing = TRUE), ]
  rownames(sticky_isoform_table) <- NULL

  gene_threshold <- .sticky_threshold(gene_degree_table$degree,
                                      sticky_method, sticky_cutoff)

  sticky_gene_table <- gene_degree_table[
    gene_degree_table$degree > gene_threshold,
    c("gene_id", "degree"),
    drop = FALSE
  ]
  sticky_gene_table <- sticky_gene_table[
    order(sticky_gene_table$degree, decreasing = TRUE), ]
  rownames(sticky_gene_table) <- NULL

  if (nrow(sticky_isoform_table) == 0L)
    message(sprintf(
      "  No sticky isoforms found above degree threshold (%.2f).",
      iso_threshold
    ))

  if (nrow(sticky_gene_table) == 0L)
    message(sprintf(
      "  No sticky genes found above degree threshold (%.2f).",
      gene_threshold
    ))

  # -- Assemble NetworkDiagnostics object -------------------------------------
  diag_obj <- new(
    "NetworkDiagnostics",
    network_stats            = network_stats,
    isoform_degree_table     = isoform_degree_table,
    gene_degree_table        = gene_degree_table,
    component_table          = component_table,
    clustering_table         = clustering_table,
    isoform_per_gene_table   = isoform_per_gene_table,
    sticky_isoform_table     = sticky_isoform_table,
    sticky_gene_table        = sticky_gene_table,
    sticky_method            = sticky_method,
    sticky_cutoff            = sticky_cutoff,
    sticky_isoform_threshold = iso_threshold,
    sticky_gene_threshold    = gene_threshold,
    created_at               = Sys.time()
  )

  validObject(diag_obj)

  obj@network_diagnostics <- diag_obj

  message(sprintf(
    "Network diagnostics complete. Access with get_network_diagnostics(obj)."
  ))

  invisible(obj)
}
