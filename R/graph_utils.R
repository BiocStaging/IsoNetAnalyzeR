## graph_utils.R
## Shared graph-building utilities used by both component_analysis.R
## and plot_network.R.
##
## Internal helpers (prefixed with ".") are not exported.

#' @include isonetanalyzer_class.R
#' @include diff_pulldown_analysis_class.R
NULL

# ==============================================================================
#  Internal: build the merged cross-gene graph
# ==============================================================================

#' Merge subgraphs across genes into a single annotated graph
#'
#' For each gene, extracts the requested direction subgraph(s), collects
#' vertices and edges, and annotates every vertex with the gene(s) whose
#' pulldown it belongs to.  Vertices appearing in more than one gene's
#' subgraph are labelled \code{"multi-gene"}.
#'
#' @param gene_results Named list of GeneDiffPulldown objects.
#' @param direction    \code{"up"}, \code{"down"}, or \code{"both"}.
#' @param level        \code{"isoform"} or \code{"gene"}.
#'
#' @return A list with:
#'   \item{graph}{igraph with vertex attributes \code{source_gene}
#'     (gene ID or \code{"multi-gene"}) and \code{is_seed} (logical).}
#'   \item{seed_names}{Character vector of all seed vertex names.}
#' @noRd
.build_cross_gene_graph <- function(gene_results, direction, level) {

  # Collect per-gene subgraphs -------------------------------------------------
  .get_graphs <- function(gr, direction) {
    graphs <- list()
    if (direction %in% c("up", "both") && !is.null(gr@up_graph))
      graphs <- c(graphs, list(gr@up_graph))
    if (direction %in% c("down", "both") && !is.null(gr@down_graph))
      graphs <- c(graphs, list(gr@down_graph))
    graphs
  }

  # For each gene, gather: vertex names, edges, and seed isoform names
  per_gene_vnames <- list()   # gene_id -> character vector of vertex names
  all_edge_dfs    <- list()
  all_seed_names  <- character(0)

  for (gid in names(gene_results)) {
    gr      <- gene_results[[gid]]
    graphs  <- .get_graphs(gr, direction)

    gene_vnames <- character(0)
    for (g in graphs) {
      gene_vnames <- union(gene_vnames, igraph::V(g)$name)

      edf <- igraph::as_data_frame(g, what = "edges")[, c("from", "to"),
                                                      drop = FALSE]
      if (nrow(edf) > 0)
        all_edge_dfs <- c(all_edge_dfs, list(edf))
    }

    per_gene_vnames[[gid]] <- gene_vnames

    # Collect seed isoforms for this gene (depending on direction)
    if (direction %in% c("up", "both"))
      all_seed_names <- union(all_seed_names, gr@up_isoforms)
    if (direction %in% c("down", "both"))
      all_seed_names <- union(all_seed_names, gr@down_isoforms)
  }

  # Determine source_gene per vertex ------------------------------------------
  all_vnames <- unique(unlist(per_gene_vnames, use.names = FALSE))

  if (length(all_vnames) == 0L)
    return(list(
      graph      = igraph::make_empty_graph(directed = FALSE),
      seed_names = character(0)
    ))

  # One row per (pulldown gene, vertex) pair; vertices are unique per gene
  member_gene   <- rep(names(per_gene_vnames), lengths(per_gene_vnames))
  member_vertex <- unlist(per_gene_vnames, use.names = FALSE)

  # A vertex reached from exactly one gene's pulldown keeps that gene
  .source_of <- function(key, member, levels) {
    n <- tabulate(match(key, levels), length(levels))
    stats::setNames(ifelse(n == 1L, member[match(levels, key)], "multi-gene"),
                    levels)
  }
  source_gene <- .source_of(member_vertex, member_gene, all_vnames)

  # Carry the gene (parent gene) vertex attribute through if available --------
  gene_attr <- rep(NA_character_, length(all_vnames))
  names(gene_attr) <- all_vnames

  for (gid in names(gene_results)) {
    gr     <- gene_results[[gid]]
    graphs <- .get_graphs(gr, direction)
    for (g in graphs) {
      g_gene <- igraph::V(g)$gene
      if (!is.null(g_gene)) {
        g_names <- igraph::V(g)$name
        still_na <- is.na(gene_attr[g_names])
        gene_attr[g_names[still_na]] <- g_gene[still_na]
      }
    }
  }

  # Build union edge data frame (deduplicated) --------------------------------
  if (length(all_edge_dfs) > 0) {
    combined_edf <- do.call(rbind, all_edge_dfs)
    # Canonical key for deduplication (undirected)
    ekey <- paste(pmin(combined_edf$from, combined_edf$to),
                  pmax(combined_edf$from, combined_edf$to),
                  sep = ":::")
    combined_edf <- combined_edf[!duplicated(ekey), , drop = FALSE]
  } else {
    combined_edf <- data.frame(from = character(0), to = character(0),
                               stringsAsFactors = FALSE)
  }

  # Assemble igraph -----------------------------------------------------------
  g <- igraph::graph_from_data_frame(
    combined_edf,
    directed = FALSE,
    vertices = data.frame(name = all_vnames, stringsAsFactors = FALSE)
  )

  igraph::V(g)$source_gene <- source_gene[igraph::V(g)$name]
  igraph::V(g)$gene        <- gene_attr[igraph::V(g)$name]
  igraph::V(g)$is_seed     <- igraph::V(g)$name %in% all_seed_names

  # Gene-level contraction if requested ----------------------------------------
  if (level == "gene") {
    g <- .contract_to_gene(g)

    # Re-derive source_gene at gene level: a gene node is "multi-gene" if
    # isoforms from different pulldown genes contributed to it
    gv_names <- igraph::V(g)$name
    parent <- unique(data.frame(gene = unname(gene_attr[member_vertex]),
                                member = member_gene))
    parent <- parent[!is.na(parent$gene), , drop = FALSE]

    igraph::V(g)$source_gene <- unname(
      .source_of(parent$gene, parent$member, gv_names)
    )

    # Re-derive is_seed at gene level
    seed_genes <- unique(stats::na.omit(
      gene_attr[all_seed_names[all_seed_names %in% names(gene_attr)]]
    ))
    igraph::V(g)$is_seed <- gv_names %in% seed_genes
  }

  list(
    graph      = g,
    seed_names = if (level == "gene") {
      igraph::V(g)$name[igraph::V(g)$is_seed]
    } else {
      all_seed_names
    }
  )
}
