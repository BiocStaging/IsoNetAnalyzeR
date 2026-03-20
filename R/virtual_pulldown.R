#' get_neighbor
#'
#' @param network_graph graph object
#' @param seed_vertices seed vertex name
#' @param order order of neighbors that should be looked at
#'
#' @returns graph object with the degree of neighbor assigned
#' @export
get_neighbor <- function(network_graph, seed_vertices, order) {
  layer_map <- vector("list", order)
  layer_map[[1]] <- igraph::V(network_graph)[name %in% seed_vertices]  # store vertex objects, not names

  # Reduce over the number of layers
  layer_map <- purrr::reduce(
    .x = seq_len(order),
    .f = function(acc, i) {
      prev <- acc[[i]]  # vertex objects in previous layer
      new <- unique(unlist(igraph::neighbors(network_graph, prev, mode = "all")))
      new <- setdiff(new, unlist(acc))  # remove already seen

      # DEBUG: check if no new neighbors found
      if (length(new) == 0) {
        message(sprintf("DEBUG: No new neighbors found at layer %d. Aborting.", i + 1))
        return(acc)
      }

      acc[[i + 1]] <- new
      acc
    },
    .init = layer_map
  )

  # Flatten to a vertex -> order mapping
  vertex_order <- unlist(
    lapply(seq_along(layer_map), function(i) {
      v <- layer_map[[i]]
      stats::setNames(rep(i - 1, length(v)), igraph::V(network_graph)[v]$name)  # use names
    })
  )

  # Assign vertex order as an attribute
  igraph::V(network_graph)$order <- NA
  igraph::V(network_graph)[names(vertex_order)]$order <- vertex_order

  network_graph
}

#' virtual_pulldown
#'
#' @param network_graph graph object
#' @param seed_vertices seed isoform
#' @param order order of neighbors that should be returned
#' @param cutoff value to remove highly connected isoforms
#'
#' @returns a subgraph (virtual pulldown)
#' @export
virtual_pulldown <- function(network_graph, seed_vertices, order, cutoff) {
  network_graph <- get_neighbor(network_graph, seed_vertices, order)

  sub_network_names <- purrr::keep(igraph::V(network_graph),
                                   ~ !is.na(.x$order) && .x$order <= order)$name

  sub_network_graph <- igraph::induced_subgraph(network_graph,
                                        vids = igraph::V(network_graph)[name %in% sub_network_names])

  mode <- "all"

  deg_full     <- igraph::degree(network_graph, v = sub_network_names, mode = mode)
  deg_internal <- igraph::degree(sub_network_graph,            mode = mode)

  res <- data.frame(
    node          = sub_network_names,
    deg_internal  = deg_internal,
    deg_full      = deg_full,
    frac_internal = ifelse(deg_full > 0, deg_internal / deg_full, NA_real_),
    stringsAsFactors = FALSE
  )

  nodes_to_keep <- res$node[!is.na(res$frac_internal) & res$frac_internal >= cutoff]

  filteredSubGraph <- igraph::induced_subgraph(sub_network_graph, vids = igraph::V(sub_network_graph)[name %in% nodes_to_keep])

  igraph::V(filteredSubGraph)$frac_internal <- res$frac_internal[match(igraph::V(filteredSubGraph)$name, res$node)]

  filteredSubGraph
}

#' plot_vp
#'
#' internal dev only
#'
#' @param network_graph vp graph object
#'
#' @returns plots a graph
plot_vp <- function(network_graph) {
  network_graph %>%
    ggraph(., layout = "igraph", algorithm = "kk") +
    geom_edge_link() +
    geom_node_point(aes(color = factor(order),
                        size = 12)) +
    geom_node_text(aes(label = name), repel = TRUE) +
    theme_void()
}
