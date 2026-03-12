#' Comprehensive network statistics
#'
#' Reports summary statistics for:
#'  - Full isoform network
#'  - Predicted interaction subnetwork
#'  - Gene-level predicted network
#'
#' @param g An igraph object. If NULL, the internal package network is used.
#' @return A named list of statistics.
#' @export
network_summary <- function(g = NULL) {

  # Use internal graph if not provided
  if (is.null(g)) {
    g <- isoform_network_graph
    if (!exists("g")) stop("Internal network not found.")
  }

  stopifnot(inherits(g, "igraph"))

  ## ---- FULL NETWORK ----
  full_isoforms <- igraph::vcount(g)
  full_genes <- length(unique(igraph::V(g)$gene))
  full_edges <- igraph::ecount(g)
  full_predicted <- sum(igraph::E(g)$interact == 1, na.rm = TRUE)

  ## ---- PREDICTED ISOFORM NETWORK ----
  #pred_edges <- igraph::E(g)[interact == 1]
  interact_vec <- igraph::get.edge.attribute(g, "interact")
  pred_edges <- igraph::E(g)[which(interact_vec == 1)]

  g_pred <- igraph::subgraph_from_edges(g,
                                        pred_edges,
                                        delete.vertices = TRUE)

  pred_isoforms <- igraph::vcount(g_pred)
  pred_genes <- length(unique(igraph::V(g_pred)$gene))
  pred_edges_n <- igraph::ecount(g_pred)

  ## ---- GENE-LEVEL NETWORK (predicted only) ----
  edges_df <- igraph::as_data_frame(g_pred, what = "edges")
  transcript_gene <- data.frame(
    enst = igraph::V(g_pred)$name,
    gene = igraph::V(g_pred)$gene,
    stringsAsFactors = FALSE
  )

  edges_df$from_gene <- transcript_gene$gene[
    match(edges_df$from, transcript_gene$enst)
  ]
  edges_df$to_gene <- transcript_gene$gene[
    match(edges_df$to, transcript_gene$enst)
  ]

  g_gene <- igraph::graph_from_data_frame(
    edges_df[, c("from_gene", "to_gene")],
    directed = FALSE
  )
  g_gene <- igraph::simplify(g_gene,
                             remove.multiple = TRUE,
                             remove.loops = TRUE)

  gene_nodes <- igraph::vcount(g_gene)
  gene_edges <- igraph::ecount(g_gene)

  list(
    full_network = list(
      isoforms = full_isoforms,
      genes = full_genes,
      total_edges = full_edges,
      predicted_edges = full_predicted
    ),
    predicted_isoform_network = list(
      isoforms = pred_isoforms,
      genes = pred_genes,
      edges = pred_edges_n
    ),
    gene_level_predicted_network = list(
      genes = gene_nodes,
      edges = gene_edges
    )
  )
}

#' Comprehensive network statistics
#'
#' Reports summary statistics for:
#'  - Full isoform network
#'  - Predicted interaction subnetwork
#'  - Gene-level predicted network
#'
#' @param g An igraph object. If NULL, the internal package network is used.
#' @return A named list of statistics.
#' @export
print_network_summary <- function(g = NULL) {
  s <- network_summary(g)

  cat("Isoform Interaction Network Summary\n")
  cat("===================================\n")

  cat("Full Network:\n")
  cat("Number of unique isoforms/genes: ",
      s$full_network$isoforms,
      "/",
      s$full_network$genes,
      "\n",
      sep = "")
  cat("Number of edges:",
      s$full_network$total_edges,
      "\n")
  cat("\n")

  cat("Predicted Isoform-Level Interaction Network:\n")
  cat("Number of unique isoforms/genes: ",
      s$predicted_isoform_network$isoforms,
      "/",
      s$predicted_isoform_network$genes,
      "\n",
      sep = "")
  cat("Number of predicted interactions:",
      s$predicted_isoform_network$edges,
      "\n")
  cat("\n")

  cat("Predicted Gene-Level Interaction Network:\n")
  cat("Number of unique genes:",
      s$gene_level_predicted_network$genes,
      "\n")
  cat("Number of predicted interactions:",
      s$gene_level_predicted_network$edges,
      "\n")
}

#' isoform_degree_distribution
#'
#' @param g graph object
#'
#' @returns number of interaction aprtners for each isoform
#' @export
isoform_degree_distribution <- function(g = NULL) {

  if (is.null(g)) {
    g <- isoform_network_graph
  }

  stopifnot(inherits(g, "igraph"))

  pred_edges <- igraph::E(g)[which(interact == 1)]

  g_pred <- igraph::subgraph_from_edges(g,
                                        pred_edges,
                                        delete.vertices = TRUE)

  deg <- igraph::degree(g_pred)

  return(deg)
}

#' gene_degree_distribution
#'
#' @param g graph object
#'
#' @returns number of interaction partner per gene
#' @export
gene_degree_distribution <- function(g = NULL) {

  if (is.null(g)) {
    g <- isoform_network_graph
  }

  stopifnot(inherits(g, "igraph"))

  pred_edges <- igraph::E(g)[which(interact == 1)]

  g_pred <- igraph::subgraph_from_edges(g,
                                        pred_edges,
                                        delete.vertices = TRUE)

  edges_df <- igraph::as_data_frame(g_pred, what = "edges")

  transcript_gene <- data.frame(
    enst = igraph::V(g)$name,
    gene = igraph::V(g)$gene,
    stringsAsFactors = FALSE
  )

  edges_df$from_gene <- transcript_gene$gene[
    match(edges_df$from, transcript_gene$enst)
  ]

  edges_df$to_gene <- transcript_gene$gene[
    match(edges_df$to, transcript_gene$enst)
  ]

  g_gene <- igraph::graph_from_data_frame(
    edges_df[, c("from_gene", "to_gene")],
    directed = FALSE
  )

  g_gene <- igraph::simplify(
    g_gene,
    remove.multiple = TRUE,
    remove.loops = TRUE
  )

  deg <- igraph::degree(g_gene)

  return(deg)
}

#' find_sticky_genes
#'
#' @param g graph object
#' @param method one of sd or quantile for cutoff calculation
#' @param cutoff cutoff to classify sticky genes
#'
#' @returns all sticky genes and the calculated cutoff number
#' @export
find_sticky_genes <- function(g = NULL, method = "sd", cutoff = 2) {

  deg <- gene_degree_distribution(g)

  if (method == "sd") {

    threshold <- mean(deg) + cutoff * stats::sd(deg)

  } else if (method == "quantile") {

    threshold <- stats::quantile(deg, cutoff)

  } else {

    stop("method must be 'sd' or 'quantile'")

  }

  sticky <- deg[deg > threshold]

  return(list(
    sticky_genes = sticky,
    threshold = threshold
  ))
}

#' find_sticky_isoforms
#'
#' @param g graph object
#' @param method one of sd or quantile
#' @param cutoff cutoff value
#'
#' @returns a list of isoforms that are classified as sticky and the used cutoff
#' @export
find_sticky_isoforms <- function(g = NULL, method = "sd", cutoff = 2) {

  deg <- isoform_degree_distribution(g)

  if (method == "sd") {

    threshold <- mean(deg) + cutoff * stats::sd(deg)

  } else if (method == "quantile") {

    threshold <- stats::quantile(deg, cutoff)

  } else {

    stop("method must be 'sd' or 'quantile'")

  }

  sticky <- deg[deg > threshold]

  return(list(
    sticky_isoforms = sticky,
    threshold = threshold
  ))
}

# Intern dev function
get_network <- function() {
  isoform_network_graph
}
