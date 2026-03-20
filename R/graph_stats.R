## graph_stats.R
## Network summary statistics and degree distribution utilities
##
## Internal helpers (starts with a ".") are not exported
## Export functions: network_summary, print_network_summary,
##  isoform_degree_distribution, gene_degree_distribution, find_sticky_isoforms,
##  find_sticky_genes


#-- Helper functions -----------------------------------------------------------

#' Extract the predicted-interaction subgraph
#'
#' Filters edges to those with \code{interact == 1} and returns the induced
#' subgraph (vertices with no remaining edges are dropped).
#'
#' @param g A validated igraph object.
#' @return An igraph object containing only predicted-interaction edges.
#' @noRd
.get_predicted_subgraph <- function(g) {
  interact_vec <- igraph::E(g, "interact")

  if (is.null(interact_vec)) {
    stop("Graph has no edge attribute 'interact'.", call. = FALSE)
  }

  pred_edges <- igraph::E(g)[which(interact_vec == 1)]

  if (length(pred_edges) == 0) {
    stop("No edges with `interact == 1` found in the graph.", call. = FALSE)
  }

  igraph::subgraph_from_edges(g, pred_edges, delete.vertices = TRUE)
}


#' Build a simplified gene-level graph from an isoform-level graph
#'
#' Maps each isoform (vertex) to its parent gene via the \code{gene} vertex
#' attribute, then constructs an undirected gene-level graph and removes
#' multi-edges and self-loops.
#'
#' @param g An igraph object with a \code{gene} vertex attribute.
#' @return A simplified, undirected igraph object at gene resolution.
#' @noRd
.build_gene_graph <- function(g) {
  if (is.null(igraph::V(g)$gene)) {
    stop("Vertex attribute 'gene' is missing from the graph.", call. = FALSE)
  }

  # Map transcript IDs -> gene IDs for both endpoints of each edge
  edges_df <- igraph::as_data_frame(g, what = "edges")
  iso_to_gene <- stats::setNames(igraph::V(g)$gene, igraph::V(g)$name)

  edges_df$from_gene <- iso_to_gene[edges_df$from]
  edges_df$to_gene   <- iso_to_gene[edges_df$to]

  g_gene <- igraph::graph_from_data_frame(
    edges_df[, c("from_gene", "to_gene")],
    directed = FALSE
  )

  # Collapse parallel edges and remove self-loops that arise from
  # isoforms of the same gene interacting with each other
  igraph::simplify(g_gene, remove.multiple = TRUE, remove.loops = TRUE)
}


#' Calculate degree-based sticky threshold and return qualifying nodes
#'
#' Shared logic for \code{find_sticky_genes} and \code{find_sticky_isoforms}.
#'
#' @param deg Named integer vector of degrees.
#' @param method \code{"sd"} or \code{"quantile"}.
#' @param cutoff Numeric scalar. Multiplier for \code{"sd"}, or probability
#'   in \[0, 1\] for \code{"quantile"}.
#' @param node_label Character scalar used in error messages ("gene" or
#'   "isoform").
#' @return A named list with \code{sticky} (named degree vector) and
#'   \code{threshold} (numeric scalar).
#' @noRd
.find_sticky <- function(deg, method, cutoff, node_label) {
  threshold <- switch(method,
                      sd       = mean(deg) + cutoff * stats::sd(deg),
                      quantile = stats::quantile(deg, cutoff),
                      stop(
                        "`method` must be 'sd' or 'quantile', got: '", method, "'",
                        call. = FALSE
                      )
  )

  sticky <- deg[deg > threshold]

  if (length(sticky) == 0) {
    message("No sticky ", node_label, "s found above threshold (", round(threshold, 2), ").")
  }

  list(sticky = sticky, threshold = threshold)
}


#-- Exported functions ---------------------------------------------------------

#' Comprehensive network summary statistics
#'
#' Reports summary statistics for three views of the network:
#' \itemize{
#'   \item The full isoform network (all edges).
#'   \item The predicted-interaction isoform subnetwork (\code{interact == 1}).
#'   \item The gene-level collapse of the predicted subnetwork.
#' }
#'
#' @param g An igraph object. If \code{NULL}, the internal graph is used.
#'
#' @return A named list with three sublists: \code{full_network},
#'   \code{predicted_isoform_network}, and \code{gene_level_predicted_network}.
#' @export
network_summary <- function(g = NULL) {
  g <- .resolve_graph(g)

  # Full network
  full_isoforms <- igraph::vcount(g)
  full_genes    <- length(unique(igraph::V(g)$gene))
  full_edges    <- igraph::ecount(g)
  full_predicted <- sum(igraph::E(g)$interact == 1, na.rm = TRUE)

  # Predicted isoform subnetwork
  g_pred      <- .get_predicted_subgraph(g)
  pred_isoforms <- igraph::vcount(g_pred)
  pred_genes    <- length(unique(igraph::V(g_pred)$gene))
  pred_edges    <- igraph::ecount(g_pred)

  # Gene-level collapse of predicted subnetwork
  g_gene     <- .build_gene_graph(g_pred)
  gene_nodes <- igraph::vcount(g_gene)
  gene_edges <- igraph::ecount(g_gene)

  list(
    full_network = list(
      isoforms        = full_isoforms,
      genes           = full_genes,
      total_edges     = full_edges,
      predicted_edges = full_predicted
    ),
    predicted_isoform_network = list(
      isoforms = pred_isoforms,
      genes    = pred_genes,
      edges    = pred_edges
    ),
    gene_level_predicted_network = list(
      genes = gene_nodes,
      edges = gene_edges
    )
  )
}


#' Print a formatted network summary
#'
#' A thin wrapper around \code{\link{network_summary}} that writes a
#' human-readable summary to the console.
#'
#' @param g An igraph object. If \code{NULL}, the internal package graph is
#'   used.
#' @return Invisibly returns the summary list produced by
#'   \code{\link{network_summary}}.
#' @export
print_network_summary <- function(g = NULL) {
  s <- network_summary(g)

  cat("Isoform Interaction Network Summary\n")
  cat("====================================\n\n")

  cat("Full Network:\n")
  cat(sprintf(
    "  Isoforms / genes : %d / %d\n",
    s$full_network$isoforms, s$full_network$genes
  ))
  cat(sprintf("  Total edges      : %d\n", s$full_network$total_edges))
  cat(sprintf("  Predicted edges  : %d\n\n", s$full_network$predicted_edges))

  cat("Predicted Isoform-Level Interaction Network:\n")
  cat(sprintf(
    "  Isoforms / genes : %d / %d\n",
    s$predicted_isoform_network$isoforms, s$predicted_isoform_network$genes
  ))
  cat(sprintf("  Predicted edges  : %d\n\n", s$predicted_isoform_network$edges))

  cat("Predicted Gene-Level Interaction Network:\n")
  cat(sprintf("  Unique genes     : %d\n", s$gene_level_predicted_network$genes))
  cat(sprintf("  Predicted edges  : %d\n", s$gene_level_predicted_network$edges))

  invisible(s)
}


#' Isoform degree distribution in the predicted subnetwork
#'
#' Returns the number of predicted interaction partners for each isoform that
#' appears in the predicted (\code{interact == 1}) subnetwork.
#'
#' @param g An igraph object. If \code{NULL}, the internal package graph is
#'   used.
#'
#' @return A named integer vector of degrees (names are transcript IDs).
#' @export
isoform_degree_distribution <- function(g = NULL) {
  g <- .resolve_graph(g)
  g_pred <- .get_predicted_subgraph(g)
  igraph::degree(g_pred)
}


#' Gene degree distribution in the predicted subnetwork
#'
#' Collapses the predicted isoform subnetwork to gene resolution, then returns
#' the number of unique gene-level interaction partners per gene.
#'
#' @param g An igraph object. If \code{NULL}, the internal package graph is
#'   used.
#'
#' @return A named integer vector of degrees (names are gene IDs).
#' @export
gene_degree_distribution <- function(g = NULL) {
  g <- .resolve_graph(g)
  g_pred <- .get_predicted_subgraph(g)
  g_gene <- .build_gene_graph(g_pred)
  igraph::degree(g_gene)
}


#' Identify sticky isoforms
#'
#' An isoform is classified as "sticky" if its degree in the predicted isoform
#' subnetwork exceeds a threshold derived from the overall degree distribution.
#'
#' @param g An igraph object. If \code{NULL}, the internal package graph is
#'   used.
#' @param method Threshold method: \code{"sd"} (mean + \code{cutoff} * SD) or
#'   \code{"quantile"} (the \code{cutoff} quantile).
#' @param cutoff Numeric scalar. SD multiplier for \code{"sd"}, or probability
#'   in \[0, 1\] for \code{"quantile"}.
#'
#' @return A named list with:
#' \describe{
#'   \item{sticky_isoforms}{Named integer vector of degrees for sticky isoforms.}
#'   \item{threshold}{The numeric threshold used.}
#' }
#' @export
find_sticky_isoforms <- function(g = NULL, method = "sd", cutoff = 2) {
  deg    <- isoform_degree_distribution(g)
  result <- .find_sticky(deg, method, cutoff, node_label = "isoform")
  list(sticky_isoforms = result$sticky, threshold = result$threshold)
}


#' Identify sticky genes
#'
#' A gene is classified as "sticky" if its degree in the predicted gene-level
#' network exceeds a threshold derived from the overall degree distribution.
#'
#' @param g An igraph object. If \code{NULL}, the internal package graph is
#'   used.
#' @param method Threshold method: \code{"sd"} (mean + \code{cutoff} * SD) or
#'   \code{"quantile"} (the \code{cutoff} quantile).
#' @param cutoff Numeric scalar. SD multiplier for \code{"sd"}, or probability
#'   in \[0, 1\] for \code{"quantile"}.
#'
#' @return A named list with:
#' \describe{
#'   \item{sticky_genes}{Named integer vector of degrees for sticky genes.}
#'   \item{threshold}{The numeric threshold used.}
#' }
#' @export
find_sticky_genes <- function(g = NULL, method = "sd", cutoff = 2) {
  deg    <- gene_degree_distribution(g)
  result <- .find_sticky(deg, method, cutoff, node_label = "gene")
  list(sticky_genes = result$sticky, threshold = result$threshold)
}
