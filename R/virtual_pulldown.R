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
  # edge_attr(g)$interact: edge_attr(g, name) first builds E(g), ~30 ms here
  interact_vec <- igraph::edge_attr(g)$interact

  if (is.null(interact_vec))
    stop("Graph has no edge attribute 'interact'.", call. = FALSE)

  # The bundled network has predicted edges only; skip copying it
  if (all(interact_vec == 1L)) return(g)

  pred_edges <- which(interact_vec == 1L)

  if (length(pred_edges) == 0L)
    return(igraph::delete_edges(g, igraph::E(g)))

  # delete.vertices = FALSE: vertices with no interact==1 edges are kept so
  # that seed isoforms with zero interactions are still found and produce a
  # single-node subgraph rather than silently disappearing.
  igraph::subgraph_from_edges(g, pred_edges, delete.vertices = FALSE)
}


#' The candidate subnetwork: vertices within \code{order} hops of a seed
#'
#' Finds the vertices within \code{order} hops of any seed with
#' \code{igraph::ego()}, takes the subgraph they induce and gives each vertex
#' its shortest-path distance from the nearest seed, computed with
#' \code{igraph::distances()} inside that subgraph. Both run at the C level.
#' An isolated seed (no predicted edge) is its own only candidate and gets
#' \code{order = 0}; \code{virtual_pulldown()} always keeps seeds.
#'
#' @param g             An igraph object (interaction-edge filtered, with
#'   isolated vertices retained).
#' @param seed_vertices Character vector of seed vertex names.
#' @param order         Integer. Maximum number of hops from any seed.
#'
#' @return A list: \code{graph}, the candidate subgraph with an \code{order}
#'   vertex attribute, and \code{vertices}, the candidates' vertex indices in
#'   \code{g} (in the subgraph's vertex order).
#' @noRd
.get_neighbor <- function(g, seed_vertices, order) {

  seed_idx <- which(igraph::V(g)$name %in% seed_vertices)

  if (length(seed_idx) == 0L)
    stop("None of the seed vertices were found in the network.",
         call. = FALSE)

  candidates <- sort(unique(unlist(
    igraph::ego(g, order = order, nodes = seed_idx, mode = "all"),
    use.names = FALSE
  )))
  # induced_subgraph() keeps the vertex order, so vertex i is candidates[i]
  sub <- igraph::induced_subgraph(g, candidates)
  # Every vertex on a shortest path of length <= order is itself within
  # `order` hops, so nearest-seed distances inside the subgraph are exact
  dist_mat <- igraph::distances(
    sub,
    v         = match(seed_idx, candidates),
    mode      = "all",
    algorithm = "unweighted"   # pure BFS — fastest option for unweighted graphs
  )
  igraph::V(sub)$order <- apply(dist_mat, 2L, min)
  list(graph = sub, vertices = candidates)
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
#' @param resolve_duplicates Logical. Resolve versioned and protein-identical
#'   redundant transcript IDs against the active network. Default \code{TRUE}.
#'
#' @return An igraph object representing the filtered virtual pulldown
#'   subnetwork. Vertices carry an \code{order} attribute (hop distance from
#'   the nearest seed) and a \code{frac_internal} attribute (fraction of
#'   their interaction edges that remain internal to the subnetwork).
#'   Seed vertices with zero \code{interact == 1} edges are retained as
#'   isolated nodes (\code{frac_internal = NA}, \code{degree = 0}),
#'   representing isoforms present in the pinned network but incident to no
#'   predicted edge in that network.
#' @examples
#' obj <- IsoNetAnalyzeR("example")
#' g <- virtual_pulldown(obj, "ENST00000621271", order = 1)
#' c(nodes = igraph::vcount(g), edges = igraph::ecount(g))
#'
#' @export
virtual_pulldown <- function(obj, seed_vertices, order, cutoff = NULL,
                             resolve_duplicates = TRUE) {

  # ── Input validation ────────────────────────────────────────────────────────
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.character(seed_vertices) || length(seed_vertices) == 0L ||
      anyNA(seed_vertices) || any(!nzchar(trimws(seed_vertices))))
    stop("`seed_vertices` must be a non-empty character vector with no blank strings.",
         call. = FALSE)

  if (!is.logical(resolve_duplicates) || length(resolve_duplicates) != 1L ||
      is.na(resolve_duplicates))
    stop("`resolve_duplicates` must be a single non-NA logical value.",
         call. = FALSE)

  resolution_audit <- NULL
  if (resolve_duplicates) {
    seed_vertices <- resolve_duplicate_isoforms(
      seed_vertices, network = obj@network
    )
    resolution_audit <- attr(seed_vertices, "dedup_audit")
    .warn_duplicate_resolution(resolution_audit)
    seed_vertices <- unique(as.character(seed_vertices))
  }

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

  # ── Steps 2-3: candidate subnetwork within `order` hops of the seeds ───────
  candidates <- .get_neighbor(g_interact, seed_vertices, order)
  g_sub      <- candidates$graph
  sub_names  <- igraph::V(g_sub)$name

  # ── Step 4: compute degrees ─────────────────────────────────────────────────
  # Degree in the full interaction graph (external + internal edges combined).
  # Isolated seeds have deg_full = 0; frac_internal is set to NA for them,
  # but they are explicitly protected in Step 5 below.
  deg_full <- igraph::degree(g_interact, v = candidates$vertices, mode = "all")

  # Degree within the candidate subnetwork (internal edges only)
  deg_internal <- igraph::degree(g_sub, mode = "all")

  frac_internal <- ifelse(deg_full > 0L, deg_internal / deg_full, NA_real_)

  # ── Step 5: apply cutoff — retain nodes with sufficient internal edges ──────
  # Seeds (order == 0) are always kept regardless of frac_internal.  The
  # cutoff is designed to remove non-specific *neighbours*, not the seeds
  # themselves.  An isolated seed (frac_internal = NA, degree = 0) represents
  # a valid result — the isoform is in the network but is incident to no
  # predicted edge in this network — and must not be silently discarded.
  is_seed   <- igraph::V(g_sub)$order == 0L
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

  if (!is.null(resolution_audit)) {
    filtered_graph <- igraph::set_graph_attr(
      filtered_graph, "duplicate_resolution_audit", resolution_audit
    )
  }

  filtered_graph
}
