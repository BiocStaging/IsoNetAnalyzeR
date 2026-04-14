## virtual_pulldown_analysis.R
## Summary statistics and per-node analysis for virtual pulldown subgraphs.
##
## Internal helpers (prefixed with ".") are not exported.
## Export functions: summarize_subgraphs, analyze_subgraph


#-- Helper functions -----------------------------------------------------------

#' Extract seed vertices from a virtual pulldown subgraph
#'
#' Returns vertices whose \code{order} attribute is \code{0}, corresponding to
#' the seeds used when the pulldown was performed.
#'
#' @param g An igraph object with an \code{order} vertex attribute.
#' @return An igraph vertex sequence of seed vertices.
#' @noRd
.get_seeds <- function(g) {
  if (is.null(igraph::V(g)$order)) {
    stop(
      "Graph has no 'order' vertex attribute. ",
      "Was this subgraph produced by `virtual_pulldown()`?",
      call. = FALSE
    )
  }

  seeds <- igraph::V(g)[order == 0]

  if (length(seeds) == 0) {
    stop(
      "No seed vertices found (no vertex with order == 0). ",
      "Was this subgraph produced by `virtual_pulldown()`?",
      call. = FALSE
    )
  }

  seeds
}


#-- Exported functions ---------------------------------------------------------

#' Summarize a list of virtual pulldown subgraphs
#'
#' Computes a set of network-level statistics for each subgraph in
#' \code{vp_list}, returning one row per subgraph. When a subgraph has
#' multiple seed vertices, their names are collapsed into a single
#' comma-separated string.
#'
#' @param vp_list A named list of igraph objects, typically the output of
#'   repeated calls to \code{\link{virtual_pulldown}}.
#'
#' @return A data frame with one row per subgraph and the following columns:
#' \describe{
#'   \item{seed}{Comma-separated seed isoform name(s).}
#'   \item{nodes}{Total number of vertices.}
#'   \item{edges}{Total number of edges.}
#'   \item{density}{Edge density of the subgraph.}
#'   \item{avg_degree}{Mean vertex degree.}
#'   \item{clustering}{Average local clustering coefficient.}
#'   \item{starness}{Degree centralization score (0 = uniform, 1 = perfect star).}
#'   \item{genes}{Number of unique genes represented.}
#' }
#' @export
summarize_subgraphs <- function(vp_list) {
  if (!is.list(vp_list) || length(vp_list) == 0) {
    stop("`vp_list` must be a non-empty list of igraph objects.", call. = FALSE)
  }

  do.call(rbind, lapply(names(vp_list), function(nm) {
    g <- vp_list[[nm]]

    if (!inherits(g, "igraph")) {
      stop(
        "Element '", nm, "' in `vp_list` is not an igraph object.",
        call. = FALSE
      )
    }

    # Collapse multiple seeds into a single readable string
    seed_str <- paste(.get_seeds(g)$name, collapse = ", ")

    data.frame(
      seed        = seed_str,
      nodes       = igraph::vcount(g),
      edges       = igraph::ecount(g),
      density     = igraph::edge_density(g),
      avg_degree  = mean(igraph::degree(g)),
      clustering  = igraph::transitivity(g, type = "average"),
      starness    = igraph::centr_degree(g, normalized = TRUE)$centralization,
      genes       = length(unique(igraph::V(g)$gene)),
      stringsAsFactors = FALSE
    )
  }))
}


#' Per-node analysis of a virtual pulldown subgraph
#'
#' Returns node-level statistics for every vertex in a virtual pulldown
#' subgraph. When multiple seed vertices are present, \code{distance_from_seed}
#' is the average shortest-path distance across all seeds.
#'
#' @param graph An igraph object.
#'
#' @return A data frame with one row per vertex and the following columns:
#' \describe{
#'   \item{isoform}{Vertex name (transcript ID).}
#'   \item{degree}{Number of edges connected to the vertex.}
#'   \item{clustering}{Local clustering coefficient.}
#'   \item{distance_from_seed}{Average shortest-path distance from all seeds.
#'     \code{Inf} indicates the vertex is unreachable from the seed(s).}
#' }
#' @export
analyze_subgraph <- function(graph) {
  if (!inherits(graph, "igraph"))
    stop("`graph` must be an igraph object.", call. = FALSE)

  seeds <- .get_seeds(graph)

  # Compute shortest-path distances from every seed; result is a
  # (n_seeds x n_vertices) matrix — average across seeds column-wise
  dist_matrix      <- igraph::distances(graph, v = seeds)
  avg_dist_from_seed <- colMeans(dist_matrix)

  data.frame(
    isoform            = igraph::V(graph)$name,
    degree             = igraph::degree(graph),
    clustering         = igraph::transitivity(graph, type = "local"),
    distance_from_seed = avg_dist_from_seed,
    stringsAsFactors   = FALSE
  )
}
