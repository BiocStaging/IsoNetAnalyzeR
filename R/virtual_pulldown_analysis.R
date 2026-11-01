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


#' Compute mean pi and entropy of pi for a single vertex
#'
#' Extracts the \code{pi} edge attribute from all edges incident to
#' \code{vertex} in \code{graph} and returns their mean and Shannon entropy.
#' Only edges already present in the subgraph are considered — these are
#' exclusively \code{interact == 1} edges, so the statistics describe the
#' confidence profile of the vertex's \emph{actual} interactions rather than
#' all potential ones.
#'
#' Shannon entropy is computed after normalising the pi values to sum to 1
#' (treating them as a discrete probability mass):
#' \deqn{H = -\sum p_i \log_2 p_i}
#' where \eqn{p_i = \pi_i / \sum \pi_j}.
#'
#' Returns \code{NA} for both metrics when the vertex has no incident edges
#' with a \code{pi} attribute, and returns \code{NA} for entropy when only
#' one edge is present (entropy of a single-element distribution is
#' trivially 0 but uninformative).
#'
#' @param graph  An igraph object with a \code{pi} edge attribute.
#' @param vertex A single vertex name (character scalar) or igraph vertex ID.
#' @return A named numeric vector with elements \code{mean_pi} and
#'   \code{pi_entropy}.
#' @noRd
.pi_stats_for_vertex <- function(graph, vertex) {
  inc_edges <- igraph::incident(graph, vertex, mode = "all")

  pi_vals <- igraph::edge_attr(graph, "pi", index = inc_edges)

  if (is.null(pi_vals) || length(pi_vals) == 0L || all(is.na(pi_vals))) {
    return(c(mean_pi = NA_real_, pi_entropy = NA_real_))
  }

  pi_vals <- pi_vals[!is.na(pi_vals)]

  mean_pi <- mean(pi_vals)

  pi_entropy <- if (length(pi_vals) < 2L) {
    NA_real_
  } else {
    p <- pi_vals / sum(pi_vals)
    # Guard against log(0): any p == 0 contributes 0 to entropy by convention
    -sum(ifelse(p > 0, p * log2(p), 0))
  }

  c(mean_pi = round(mean_pi, 4L), pi_entropy = round(pi_entropy, 4L))
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
#' @examples
#' obj <- IsoNetAnalyzeR("example")
#' neighbourhoods <- list(
#'   cltcl1_up = virtual_pulldown(obj, "ENST00000621271", order = 1),
#'   cltcl1_down = virtual_pulldown(obj, "ENST00000427926", order = 1)
#' )
#' summarize_subgraphs(neighbourhoods)
#'
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

    n_nodes <- igraph::vcount(g)
    data.frame(
      seed        = seed_str,
      nodes       = n_nodes,
      edges       = igraph::ecount(g),
      density     = if (n_nodes < 2L) 0 else igraph::edge_density(g),
      avg_degree  = mean(igraph::degree(g)),
      clustering  = if (n_nodes < 2L) 0 else
        igraph::transitivity(g, type = "average", isolates = "zero"),
      starness    = if (n_nodes < 2L) NA_real_ else
        igraph::centr_degree(g, normalized = TRUE)$centralization,
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
#' Pi statistics (\code{mean_pi}, \code{pi_entropy}) are computed only for
#' seed vertices (order == 0), where they characterise the confidence profile
#' of the interactions that anchored the pulldown.  For non-seed vertices these
#' columns are \code{NA} — their edges were pulled in by BFS expansion and the
#' pi values reflect the network globally rather than any seed-specific signal.
#'
#' @param graph An igraph object with a \code{pi} edge attribute (preserved
#'   from the full network after \code{interact == 1} filtering).
#'
#' @return A data frame with one row per vertex and the following columns:
#' \describe{
#'   \item{isoform}{Vertex name (transcript ID).}
#'   \item{is_seed}{Logical flag copied directly from the graph's
#'     \code{order == 0} annotation.}
#'   \item{degree}{Number of edges connected to the vertex.}
#'   \item{clustering}{Local clustering coefficient.}
#'   \item{distance_from_seed}{Average shortest-path distance from all seeds.
#'     \code{Inf} indicates the vertex is unreachable from the seed(s).}
#'   \item{mean_pi}{Mean of the \code{pi} edge-attribute values across all
#'     incident edges.  \code{NA} for non-seed vertices.}
#'   \item{pi_entropy}{Shannon entropy of the normalised \code{pi} values
#'     across all incident edges (bits).  Measures how evenly interaction
#'     confidence is spread across a seed's partners: low entropy indicates a
#'     few high-confidence interactions; high entropy indicates many edges of
#'     similar moderate confidence.  \code{NA} for non-seed vertices and for
#'     seed vertices with only one incident edge.}
#' }
#' @examples
#' obj <- IsoNetAnalyzeR("example")
#' g <- virtual_pulldown(obj, "ENST00000621271", order = 1)
#' head(analyze_subgraph(g))
#'
#' @export
analyze_subgraph <- function(graph) {
  if (!inherits(graph, "igraph"))
    stop("`graph` must be an igraph object.", call. = FALSE)

  seeds <- .get_seeds(graph)

  # Compute shortest-path distances from every seed; result is a
  # (n_seeds x n_vertices) matrix — average across seeds column-wise
  dist_matrix        <- igraph::distances(graph, v = seeds)
  avg_dist_from_seed <- colMeans(dist_matrix)

  # Pi statistics for seed vertices only
  has_pi      <- !is.null(igraph::edge_attr(graph, "pi"))
  seed_names  <- seeds$name
  all_names   <- igraph::V(graph)$name

  pi_stats <- lapply(all_names, function(nm) {
    if (nm %in% seed_names && has_pi) {
      .pi_stats_for_vertex(graph, nm)
    } else {
      c(mean_pi = NA_real_, pi_entropy = NA_real_)
    }
  })

  data.frame(
    isoform            = all_names,
    is_seed            = all_names %in% seed_names,
    degree             = igraph::degree(graph),
    clustering         = igraph::transitivity(
      graph, type = "local", isolates = "zero"
    ),
    distance_from_seed = avg_dist_from_seed,
    mean_pi            = vapply(pi_stats, `[[`, numeric(1), "mean_pi"),
    pi_entropy         = vapply(pi_stats, `[[`, numeric(1), "pi_entropy"),
    stringsAsFactors   = FALSE
  )
}
