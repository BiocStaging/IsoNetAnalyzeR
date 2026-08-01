## virtual_pulldown.R
## Functions for performing a virtual pulldown on an isoform interaction network.
##
## Internal helpers (prefixed with ".") are not exported.
## Export functions: virtual_pulldown


#-- Helper functions -----------------------------------------------------------

#' Extract the predicted-interaction subgraph
#'
#' Filters the network to edges where \code{interact == 1} and returns the
#' induced subgraph.  Isolated vertices — those present in the full network
#' but carrying no \code{interact == 1} edges — are \strong{retained}.  This
#' is intentional: an isoform that exists in the network but has no predicted
#' edge is a valid single-node result and must remain reachable by the
#' seed-matching and cutoff logic downstream.
#'
#' @param g   An igraph object with an \code{interact} edge attribute.
#' @return    An igraph object containing only predicted-interaction edges,
#'   with all original vertices preserved (including those now isolated).
#' @noRd
.filter_interact_edges <- function(g) {
  interact_vec <- igraph::edge_attr(g, "interact")

  if (is.null(interact_vec))
    stop("Graph has no edge attribute 'interact'.", call. = FALSE)

  pred_edges <- which(interact_vec == 1L)

  if (length(pred_edges) == 0L)
    return(igraph::delete_edges(g, igraph::E(g)))

  # delete.vertices = FALSE: vertices with no interact==1 edges are kept so
  # that seed isoforms with zero interactions are still found and produce a
  # single-node subgraph rather than silently disappearing.
  igraph::subgraph_from_edges(g, pred_edges, delete.vertices = FALSE)
}


#' Annotate vertices with their BFS order relative to seed vertices
#'
#' Performs a single BFS from all seed vertices simultaneously using
#' \code{igraph::distances()}, which runs entirely at the C level and is
#' substantially faster than looping over \code{igraph::neighbors()} in R.
#' Each vertex is assigned the shortest-path distance from the nearest seed.
#' Vertices not reachable within \code{order} hops are left as \code{NA}.
#'
#' Isolated seed vertices (degree zero in the interaction-filtered graph)
#' receive \code{order = 0} automatically because their distance from
#' themselves is 0, which satisfies \code{min_dist <= order} for any
#' \code{order >= 0}.  They will therefore always be included in the
#' candidate subnetwork and are subsequently protected from the cutoff filter
#' in \code{virtual_pulldown()}.
#'
#' @param g             An igraph object (interaction-edge filtered, with
#'   isolated vertices retained).
#' @param seed_vertices Character vector of seed vertex names.
#' @param order         Integer. Maximum number of hops from any seed.
#'
#' @return The input graph with an \code{order} vertex attribute added.
#' @noRd
.get_neighbor <- function(g, seed_vertices, order) {

  seed_idx <- which(igraph::V(g)$name %in% seed_vertices)

  if (length(seed_idx) == 0L)
    stop("None of the seed vertices were found in the network.",
         call. = FALSE)

  # igraph::distances() runs BFS entirely in C — one call covers all seeds
  # simultaneously, returning a (n_seeds x n_vertices) matrix.
  # Column-wise minimum gives each vertex's distance from its nearest seed.
  # For isolated seeds, distances() returns 0 for the seed itself and Inf
  # for all other vertices — so isolated seeds correctly get order = 0.
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
#'   Seed vertices with zero \code{interact == 1} edges are retained as
#'   isolated nodes (\code{frac_internal = NA}, \code{degree = 0}),
#'   representing isoforms present in the pinned network but incident to no
#'   predicted edge in that network.
#' @export
virtual_pulldown <- function(obj, seed_vertices, order, cutoff = NULL) {

  # ── Input validation ────────────────────────────────────────────────────────
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.character(seed_vertices) || length(seed_vertices) == 0L ||
      anyNA(seed_vertices) || any(!nzchar(trimws(seed_vertices))))
    stop("`seed_vertices` must be a non-empty character vector with no blank strings.",
         call. = FALSE)

  order <- .validate_positive_integer(order, "order")

  # ── Cutoff resolution ───────────────────────────────────────────────────────
  if (is.null(cutoff)) {
    cutoff <- obj@metadata[["network_cutoff"]]
    if (is.null(cutoff))
      stop(
        "`cutoff` not supplied and `obj@metadata$network_cutoff` is not set. ",
        "Please provide a numeric cutoff explicitly.",
        call. = FALSE
      )
    cutoff <- .validate_cutoff(cutoff, "obj@metadata$network_cutoff")
  } else {
    cutoff <- .validate_cutoff(cutoff)
  }

  # ── Step 1: restrict to interaction edges upfront ──────────────────────────
  # Working on this smaller graph means BFS and degree calculations never
  # touch below-threshold edges — the single biggest performance win.
  # Isolated vertices (isoforms in the network with no interact==1 edges) are
  # deliberately retained so that zero-interaction seeds produce a valid
  # single-node subgraph rather than disappearing before BFS.
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
  # Degree in the full interaction graph (external + internal edges combined).
  # Isolated seeds have deg_full = 0; frac_internal is set to NA for them,
  # but they are explicitly protected in Step 5 below.
  deg_full <- igraph::degree(g_interact, v = sub_idx, mode = "all")

  # Degree within the candidate subnetwork (internal edges only)
  g_sub        <- igraph::induced_subgraph(g_interact, vids = sub_idx)
  deg_internal <- igraph::degree(g_sub, mode = "all")

  frac_internal <- ifelse(deg_full > 0L, deg_internal / deg_full, NA_real_)

  # ── Step 5: apply cutoff — retain nodes with sufficient internal edges ──────
  # Seeds (order == 0) are always kept regardless of frac_internal.  The
  # cutoff is designed to remove non-specific *neighbours*, not the seeds
  # themselves.  An isolated seed (frac_internal = NA, degree = 0) represents
  # a valid result — the isoform is in the network but is incident to no
  # predicted edge in this network — and must not be silently discarded.
  is_seed   <- igraph::V(g_interact)$order[sub_idx] == 0L
  keep_mask <- is_seed | (!is.na(frac_internal) & frac_internal >= cutoff)
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
