## virtual_pulldown.R
## Functions for performing a virtual pulldown on an isoform interaction network.
##
## Internal helpers (prefixed with ".") are not exported.
## Export functions: virtual_pulldown


#-- Helper functions -----------------------------------------------------------

#' Extract the predicted-interaction subgraph
#'
#' Filters the network to edges where \code{interact == 1} and returns the
#' induced subgraph. Isolated vertices are dropped. This is done once upfront
#' so that all BFS expansion and degree calculations operate only on
#' interaction edges, never on below-threshold edges.
#'
#' @param g   An igraph object with an \code{interact} edge attribute.
#' @return    An igraph object containing only predicted-interaction edges,
#'   with isolated vertices removed.
#' @noRd
.filter_interact_edges <- function(g) {
  interact_vec <- igraph::edge_attr(g, "interact")

  if (is.null(interact_vec))
    stop("Graph has no edge attribute 'interact'.", call. = FALSE)

  pred_edges <- which(interact_vec == 1L)

  if (length(pred_edges) == 0L)
    stop("No edges with `interact == 1` found in the graph.", call. = FALSE)

  igraph::subgraph_from_edges(g, pred_edges, delete.vertices = TRUE)
}


#' Annotate vertices with their BFS order relative to seed vertices
#'
#' Performs a single BFS from all seed vertices simultaneously using
#' \code{igraph::distances()}, which runs entirely at the C level and is
#' substantially faster than looping over \code{igraph::neighbors()} in R.
#' Each vertex is assigned the shortest-path distance from the nearest seed.
#' Vertices not reachable within \code{order} hops are left as \code{NA}.
#'
#' @param g             An igraph object (interaction-edge filtered).
#' @param seed_vertices Character vector of seed vertex names.
#' @param order         Integer. Maximum number of hops from any seed.
#'
#' @return The input graph with an \code{order} vertex attribute added.
#' @noRd
.get_neighbor <- function(g, seed_vertices, order) {

  seed_idx <- which(igraph::V(g)$name %in% seed_vertices)

  if (length(seed_idx) == 0L)
    stop("None of the seed vertices were found in the (interaction-filtered) network.",
         call. = FALSE)

  # igraph::distances() runs BFS entirely in C — one call covers all seeds
  # simultaneously, returning a (n_seeds x n_vertices) matrix.
  # Column-wise minimum gives each vertex's distance from its nearest seed.
  dist_mat <- igraph::distances(
    g,
    v         = seed_idx,
    mode      = "all",
    algorithm = "unweighted"   # pure BFS — fastest option for unweighted graphs
  )

  min_dist <- if (is.matrix(dist_mat)) apply(dist_mat, 2L, min) else as.vector(dist_mat)

  # Vertices unreachable within `order` hops stay NA
  igraph::V(g)$order <- ifelse(min_dist <= order, min_dist, NA_real_)
  g
}


#-- Exported functions ---------------------------------------------------------

#' Virtual pulldown
#'
#' Performs a virtual pulldown by first restricting the network to predicted
#' interaction edges (\code{interact == 1}), expanding from a set of seed
#' isoforms up to \code{order} hops, then filtering out non-specific nodes
#' whose edges are predominantly external to the subnetwork.
#'
#' The network is taken directly from the \linkS4class{IsoNetAnalyzeR} session
#' object. The cutoff defaults to the value stored in
#' \code{obj@@metadata[["network_cutoff"]]} but can be overridden to explore
#' alternative thresholds.
#'
#' @param obj           An \linkS4class{IsoNetAnalyzeR} session object.
#' @param seed_vertices Character vector of seed isoform vertex names.
#'   All supplied vertices are treated as order-0 anchors.
#' @param order         Integer. Number of interaction hops to expand from
#'   the seeds.
#' @param cutoff        Numeric scalar in \[0, 1\]. Minimum fraction of edges
#'   that must be internal to the subnetwork for a node to be retained. When
#'   \code{NULL} (default) the value stored in
#'   \code{obj@@metadata[["network_cutoff"]]} is used.
#'
#' @return An igraph object representing the filtered virtual pulldown
#'   subnetwork. Vertices carry an \code{order} attribute (hop distance from
#'   the nearest seed) and a \code{frac_internal} attribute (fraction of
#'   their interaction edges that remain internal to the subnetwork).
#' @export
virtual_pulldown <- function(obj, seed_vertices, order, cutoff = NULL) {

  # ── Input validation ────────────────────────────────────────────────────────
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.character(seed_vertices) || length(seed_vertices) == 0L ||
      any(nchar(trimws(seed_vertices)) == 0L))
    stop("`seed_vertices` must be a non-empty character vector with no blank strings.",
         call. = FALSE)

  if (!is.numeric(order) || length(order) != 1L || order < 1L)
    stop("`order` must be a positive integer scalar.", call. = FALSE)

  order <- as.integer(order)

  # ── Cutoff resolution ───────────────────────────────────────────────────────
  if (is.null(cutoff)) {
    cutoff <- obj@metadata[["network_cutoff"]]
    if (is.null(cutoff) || is.na(cutoff))
      stop(
        "`cutoff` not supplied and `obj@metadata$network_cutoff` is not set. ",
        "Please provide a numeric cutoff explicitly.",
        call. = FALSE
      )
  } else {
    if (!is.numeric(cutoff) || length(cutoff) != 1L ||
        cutoff < 0 || cutoff > 1)
      stop("`cutoff` must be a single numeric value in [0, 1].", call. = FALSE)
  }

  # ── Step 1: restrict to interaction edges upfront ──────────────────────────
  # Working on this smaller graph means BFS and degree calculations never
  # touch below-threshold edges — the single biggest performance win.
  g_interact <- .filter_interact_edges(obj@network)

  # ── Step 2: BFS expansion — single C-level distances() call ────────────────
  g_interact <- .get_neighbor(g_interact, seed_vertices, order)

  # ── Step 3: extract subnetwork vertices via integer indices (fast) ──────────
  sub_idx <- which(!is.na(igraph::V(g_interact)$order))

  if (length(sub_idx) == 0L)
    stop("No vertices found within the requested hop distance from the seeds.",
         call. = FALSE)

  sub_names <- igraph::V(g_interact)$name[sub_idx]

  # ── Step 4: compute degrees — one induced_subgraph call, integer indices ────
  # Degree in the full interaction graph (external + internal edges combined)
  deg_full <- igraph::degree(g_interact, v = sub_idx, mode = "all")

  # Degree within the candidate subnetwork (internal edges only)
  g_sub        <- igraph::induced_subgraph(g_interact, vids = sub_idx)
  deg_internal <- igraph::degree(g_sub, mode = "all")

  frac_internal <- ifelse(deg_full > 0L, deg_internal / deg_full, NA_real_)

  # ── Step 5: apply cutoff — retain nodes with sufficient internal edges ──────
  keep_mask  <- !is.na(frac_internal) & frac_internal >= cutoff
  keep_names <- sub_names[keep_mask]

  if (length(keep_names) == 0L)
    stop(
      "No nodes passed the cutoff (", round(cutoff, 4L), "). ",
      "Consider lowering the cutoff or increasing the expansion order.",
      call. = FALSE
    )

  keep_idx_sub   <- which(igraph::V(g_sub)$name %in% keep_names)
  filtered_graph <- igraph::induced_subgraph(g_sub, vids = keep_idx_sub)

  # ── Step 6: attach frac_internal to the filtered graph ─────────────────────
  igraph::V(filtered_graph)$frac_internal <- frac_internal[
    match(igraph::V(filtered_graph)$name, sub_names)
  ]

  filtered_graph
}
