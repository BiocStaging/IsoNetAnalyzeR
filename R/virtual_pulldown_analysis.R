
#' summarize_subgraphs
#'
#' @param vp_list
#'
#' @returns
#' @export
summarize_subgraphs <- function(vp_list) {

  do.call(rbind, lapply(names(vp_list), function(name) {

    g <- vp_list[[name]]

    data.frame(
      seed = igraph::V(g)[order == 0]$name,
      nodes = igraph::vcount(g),
      edges = igraph::ecount(g),
      density = igraph::edge_density(g),
      avg_degree = mean(igraph::degree(g)),
      clustering = igraph::transitivity(g, type = "average"),
      starness = igraph::centr_degree(g, normalized = TRUE)$centralization,
      genes = length(unique(igraph::V(g)$gene))
    )

  }))
}


#' analyze_subgraph
#'
#' @param graph
#'
#' @returns
#' @export
analyze_subgraph <- function(graph) {
  seed = igraph::V(graph)[order == 0]$name

  data.frame(
    isoform = igraph::V(graph)$name,
    degree = igraph::degree(graph),
    clustering = igraph::transitivity(graph, type = "local"),
    distance_from_seed = igraph::distances(graph, v = seed)[1,]
  )

}
