## virtual_pulldown.R
## Functions for performing a virtual pulldown on an isoform interaction network.
##
## Internal helpers (prefixed with ".") are not exported.
## Public API: virtual_pulldown


#-- Helper functions -----------------------------------------------------------

#' Annotate vertices with their neighbor order relative to seed vertices
#'
#' Performs a breadth-first expansion from \code{seed_vertices} up to
#' \code{order} hops and assigns each reached vertex its distance from the
#' seed as the \code{order} vertex attribute. Vertices not reached are left
#' as \code{NA}.
#'
#' @param network_graph An igraph object.
#' @param seed_vertices Character vector of seed vertex names.
#' @param order Integer. Maximum number of hops from the seed.
#'
#' @return The input graph with an \code{order} vertex attribute added.
#' @noRd
.get_neighbor <- function(network_graph, seed_vertices, order) {

  # Initialise layer map: layer 1 holds the seed vertices (order = 0)
  layer_map <- vector("list", order + 1)
  layer_map[[1]] <- igraph::V(network_graph)[name %in% seed_vertices]

  # Iteratively expand one hop at a time, tracking already-visited vertices
  layer_map <- purrr::reduce(
    .x = seq_len(order),
    .f = function(acc, i) {
      prev <- acc[[i]]
      new  <- unique(unlist(igraph::neighbors(network_graph, prev, mode = "all")))
      # Drop vertices already assigned to an earlier layer
      new  <- setdiff(new, unlist(acc))

      if (length(new) == 0) {
        message(sprintf("No new neighbors found at layer %d. Stopping expansion.", i + 1))
        return(acc)
      }

      acc[[i + 1]] <- new
      acc
    },
    .init = layer_map
  )

  # Flatten layer map to a named vector of order values (0 = seed)
  vertex_order <- unlist(
    lapply(seq_along(layer_map), function(i) {
      v <- layer_map[[i]]
      stats::setNames(
        rep(i - 1, length(v)),
        igraph::V(network_graph)[v]$name
      )
    })
  )

  # Assign order as a vertex attribute; unreached vertices remain NA
  igraph::V(network_graph)$order <- NA
  igraph::V(network_graph)[names(vertex_order)]$order <- vertex_order

  network_graph
}


#-- Exported functions ---------------------------------------------------------

#' Virtual pulldown
#'
#' Performs a virtual pulldown by expanding a set of seed isoforms up to
#' \code{order} hops in the network, then filtering out highly connected
#' (non-specific) nodes based on the fraction of their edges that remain
#' internal to the subnetwork.
#'
#' @param network_graph An igraph object. If \code{NULL}, the internal package
#'   graph is used.
#' @param seed_vertices Character vector of seed isoform vertex names.
#' @param order Integer. Number of interaction hops to expand from the seed.
#' @param cutoff Numeric scalar in \[0, 1\]. Minimum fraction of edges that must
#'   be internal to the subnetwork for a node to be retained. Nodes below this
#'   threshold are considered non-specific and removed.
#'
#' @return An igraph object representing the filtered virtual pulldown
#'   subnetwork. Vertices carry a \code{frac_internal} attribute indicating
#'   their fraction of internal edges.
#' @export
virtual_pulldown <- function(network_graph = NULL, seed_vertices, order, cutoff) {

  network_graph <- .resolve_graph(network_graph)

  # Annotate all vertices with their hop distance from the seed
  network_graph <- .get_neighbor(network_graph, seed_vertices, order)

  # Extract names of all vertices within the requested hop distance
  sub_network_names <- purrr::keep(
    igraph::V(network_graph),
    ~ !is.na(.x$order) && .x$order <= order
  )$name

  # Build the initial subgraph from those vertices
  sub_network_graph <- igraph::induced_subgraph(
    network_graph,
    vids = igraph::V(network_graph)[name %in% sub_network_names]
  )

  # Compute degree in the full network and internally to identify sticky nodes
  deg_full     <- igraph::degree(network_graph,     v    = sub_network_names, mode = "all")
  deg_internal <- igraph::degree(sub_network_graph, mode = "all")

  frac_internal <- ifelse(deg_full > 0, deg_internal / deg_full, NA_real_)

  # Retain only nodes whose edges are sufficiently internal (>= cutoff)
  nodes_to_keep <- sub_network_names[!is.na(frac_internal) & frac_internal >= cutoff]

  filtered_graph <- igraph::induced_subgraph(
    sub_network_graph,
    vids = igraph::V(sub_network_graph)[name %in% nodes_to_keep]
  )

  # Attach frac_internal as a vertex attribute on the filtered graph
  igraph::V(filtered_graph)$frac_internal <- frac_internal[
    match(igraph::V(filtered_graph)$name, sub_network_names)
  ]

  filtered_graph
}
