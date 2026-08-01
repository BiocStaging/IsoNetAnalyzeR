## component_analysis.R
## S4 class and user-facing function for connected-component analysis of a
## DifferentialPulldownAnalysis result.
##
## Results are stored inside the parent DifferentialPulldownAnalysis object
## (in its @component_analyses slot), which is then written back into the
## IsoNetAnalyzeR session object.
##
## Tier 1 (structural):
##   component_table  - size, seed enrichment, direction mixing, shared nodes
##   node_table       - per-node membership, role, degree
##
## Tier 2 (topological):
##   topology_table   - diameter, radius, degree centralization per component
##   hub_table        - betweenness centrality for shared (multi-gene) nodes,
##                      ranked within their component
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: ComponentAnalysis class, run_component_analysis,
##           show(ComponentAnalysis), get_component_overview,
##           get_component

#' @include isonetanalyzer_class.R
#' @include diff_pulldown_analysis_class.R
#' @include graph_utils.R
NULL


# ==============================================================================
#  S4 class: ComponentAnalysis
# ==============================================================================

#' ComponentAnalysis - connected-component analysis of a differential pulldown
#'
#' Produced by \code{\link{run_component_analysis}}.  Stores Tier-1 structural
#' and Tier-2 topological summary tables derived from decomposing the merged
#' cross-gene network of a \linkS4class{DifferentialPulldownAnalysis} run.
#'
#' @slot label            Character scalar.  User-supplied label, used as the
#'   \code{component_analysis} argument in
#'   \code{get_component_overview(obj, analysis, component_analysis)}.
#' @slot parent_label     Character scalar.  Label of the parent
#'   \linkS4class{DifferentialPulldownAnalysis} run.
#'
#' @slot component_table  Data frame (Tier 1).  One row per retained connected
#'   component. Columns:
#'   \code{component_id}, \code{n_nodes}, \code{n_edges}, \code{density},
#'   \code{avg_degree}, \code{clustering},
#'   \code{n_seeds}, \code{n_neighbours}, \code{seed_fraction},
#'   \code{n_seeds_up}, \code{n_seeds_down},
#'   \code{n_genes_up}, \code{n_genes_down},
#'   \code{mixed} (logical), \code{n_shared_nodes}, \code{shared_fraction}.
#'
#' @slot node_table       Data frame (Tier 1).  One row per node across all
#'   retained components. Columns:
#'   \code{node}, \code{component_id}, \code{source_gene}, \code{is_seed},
#'   \code{seed_direction} (\code{"up"}/\code{"down"}/\code{"both"}/\code{NA}),
#'   \code{degree}.
#'
#' @slot topology_table   Data frame (Tier 2).  One row per retained component.
#'   Augments \code{component_table} with:
#'   \code{component_id}, \code{diameter}, \code{radius},
#'   \code{degree_centralization},
#'   \code{topology_class} (\code{"clique-like"}, \code{"hub-spoke"}, or
#'   \code{"chain-like"} - see Details).
#'
#' @slot hub_table        Data frame (Tier 2).  One row per shared
#'   (\code{source_gene == "multi-gene"}) node across all retained components,
#'   ranked by betweenness centrality within its component. Columns:
#'   \code{node}, \code{component_id}, \code{degree},
#'   \code{betweenness} (raw), \code{betweenness_norm} (normalised to
#'   \[0, 1\] within each component), \code{betweenness_rank} (1 = highest).
#'
#' @slot direction        Character scalar: \code{"up"}, \code{"down"}, or
#'   \code{"both"}.
#' @slot level            Character scalar: \code{"isoform"} or \code{"gene"}.
#' @slot ora_results   Named list of \linkS4class{ComponentOraResult}
#'   objects, keyed by user-supplied label.  Populated by
#'   \code{\link{run_component_ora}}.  Empty list by default.
#' @slot genes            Character vector.  Gene IDs included (length 0 = all).
#' @slot created_at       \code{POSIXct} scalar.  Timestamp of the run.
#'
#' @details
#' \strong{Topology classification} (\code{topology_class}):
#' Each component is heuristically classified based on its degree
#' centralization score:
#' \describe{
#'   \item{\code{"hub-spoke"}}{Centralization > 0.6: degree is concentrated
#'     in one or a few nodes.}
#'   \item{\code{"clique-like"}}{Centralization < 0.2: connectivity is
#'     relatively uniformly distributed.}
#'   \item{\code{"chain-like"}}{Centralization between 0.2 and 0.6: a
#'     intermediate degree concentration.}
#' }
#' These thresholds are descriptive heuristics and do not establish a
#' regulator, complex, pathway, or mechanism.
#'
#' @exportClass ComponentAnalysis
setClass(
  "ComponentAnalysis",
  representation(
    label            = "character",
    parent_label     = "character",
    # Tier 1
    component_table  = "data.frame",
    node_table       = "data.frame",
    # Tier 2
    topology_table   = "data.frame",
    hub_table        = "data.frame",
    # Params
    direction        = "character",
    level            = "character",
    genes            = "character",
    ora_results      = "list",
    created_at       = "POSIXct"
  )
)


# ==============================================================================
#  Validator
# ==============================================================================

setValidity("ComponentAnalysis", function(object) {
  msgs <- character()

  if (length(object@label) != 1L || nchar(trimws(object@label)) == 0L)
    msgs <- c(msgs, "@label must be a non-empty character scalar.")

  if (length(object@parent_label) != 1L ||
      nchar(trimws(object@parent_label)) == 0L)
    msgs <- c(msgs, "@parent_label must be a non-empty character scalar.")

  if (!object@direction %in% c("up", "down", "both"))
    msgs <- c(msgs, "@direction must be 'up', 'down', or 'both'.")

  if (!object@level %in% c("isoform", "gene"))
    msgs <- c(msgs, "@level must be 'isoform' or 'gene'.")

  if (length(object@created_at) != 1L)
    msgs <- c(msgs, "@created_at must be a length-1 POSIXct.")

  if (length(msgs)) msgs else TRUE
})


# ==============================================================================
#  show()
# ==============================================================================

#' @describeIn ComponentAnalysis Print a concise summary.
#' @param object A \linkS4class{ComponentAnalysis} object.
setMethod("show", "ComponentAnalysis", function(object) {
  ct <- object@component_table
  tt <- object@topology_table
  ht <- object@hub_table

  n_comp  <- nrow(ct)
  n_mixed <- if (n_comp > 0L) sum(ct$mixed, na.rm = TRUE) else 0L

  cli::cli_rule("ComponentAnalysis")
  cli::cli_text("Label        : {object@label}")
  cli::cli_text("Parent run   : {object@parent_label}")
  cli::cli_text("Created      : {format(object@created_at, '%Y-%m-%d %H:%M:%S')}")
  cli::cli_text("Direction    : {object@direction}    Level: {object@level}")
  cli::cli_text("Genes        : {if (length(object@genes) == 0L) 'all' else length(object@genes)}")
  cli::cli_text("Components   : {n_comp} retained  ({n_mixed} mixed up/down)")

  if (n_comp > 0L) {
    cli::cli_text(
      "Size range   : {min(ct$n_nodes)} - {max(ct$n_nodes)} nodes per component"
    )

    if (nrow(tt) > 0L) {
      tclass_counts <- table(tt$topology_class)
      tclass_str    <- paste(
        names(tclass_counts), tclass_counts, sep = ": ", collapse = "  |  "
      )
      cli::cli_text("Topology     : {tclass_str}")
    }

    if (nrow(ht) > 0L)
      cli::cli_text(
        "Shared hubs  : {nrow(ht)} scored  (top: {ht$node[1]}, betweenness_norm = {round(ht$betweenness_norm[1], 3)})"
      )
  }
  cli::cli_rule()
  invisible(object)
})


# ==============================================================================
#  Exported: run_component_analysis
# ==============================================================================

#' Run a connected-component analysis on a differential pulldown result
#'
#' Decomposes the merged cross-gene network of an existing
#' \linkS4class{DifferentialPulldownAnalysis} run into its connected components
#' and computes both Tier-1 structural statistics and Tier-2 topological
#' statistics for each component.  The result is stored inside the parent
#' analysis and can be accessed with
#' \code{get_component_overview(obj, analysis, label)}.
#'
#' @section Tier-1 statistics (structural):
#' \itemize{
#'   \item \strong{Size}: nodes, edges, density, average degree, clustering
#'     coefficient.
#'   \item \strong{Seed enrichment}: count and fraction of seed vs. neighbour
#'     nodes.
#'   \item \strong{Direction mixing}: up/down seed counts and contributing
#'     genes; the \code{mixed} flag records components containing both labels.
#'   \item \strong{Multi-gene hubs}: count and fraction of shared nodes.
#' }
#'
#' @section Tier-2 statistics (topological):
#' \itemize{
#'   \item \strong{Diameter}: longest shortest path in the component.  Large
#'     diameter suggests a signalling cascade rather than a tight complex.
#'   \item \strong{Radius}: minimum eccentricity across all nodes.
#'   \item \strong{Degree centralization}: how unequally degree is distributed
#'     (0 = perfectly uniform, 1 = pure star graph).
#'   \item \strong{Topology class}: heuristic classification into
#'     \code{"hub-spoke"}, \code{"chain-like"}, or \code{"clique-like"} based
#'     on degree centralization (see \linkS4class{ComponentAnalysis} Details).
#'   \item \strong{Hub betweenness}: for every shared (multi-gene) node,
#'     betweenness centrality computed within its component, normalised to
#'     \[0, 1\], and ranked.  Nodes with high normalised betweenness are
#'     structural bridges across the sub-communities pulled in by different
#'     genes. Biological relevance requires external evidence.
#' }
#'
#' @param obj          An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis     Character scalar.  Label of an existing
#'   \linkS4class{DifferentialPulldownAnalysis}; use \code{list_analyses(obj)} to see available labels.
#' @param label        Character scalar.  Key under which to store the result.
#' @param direction    \code{"both"} (default), \code{"up"}, or \code{"down"}.
#' @param level        \code{"isoform"} (default) or \code{"gene"}.
#' @param genes        Character vector of gene IDs to include, or \code{NULL}
#'   for all genes.
#' @param min_component_nodes Integer or \code{NULL}.  Drop components with
#'   fewer nodes than this.
#' @param min_shared_nodes    Integer or \code{NULL}.  Drop components with
#'   fewer shared (multi-gene) nodes than this.
#' @param min_edges           Integer or \code{NULL}.  Drop components with
#'   fewer edges than this.
#'
#' @return An updated \linkS4class{IsoNetAnalyzeR} object (invisibly).
#'   Access results with \code{\link{get_component_overview}} (combined
#'   Tier-1 and Tier-2 table) and \code{\link{get_component}} (per-node
#'   table for a single component).
#'
#' @examples
#' \dontrun{
#' obj <- IsoNetAnalyzeR("demo")
#' obj <- run_differential_pulldown(obj, "cond_A", "de_isoforms.csv", order = 1)
#'
#' # Default: isoform level, both directions, all genes
#' obj <- run_component_analysis(obj, "cond_A", label = "comp_both")
#'
#' # show() summary
#' obj   # the updated session object; the component analysis is listed
#'
#' # Combined Tier-1 + Tier-2 overview (one row per component)
#' overview <- get_component_overview(obj, "cond_A", "comp_both")
#'
#' # Components containing both up and down seed labels
#' overview[overview$mixed, ]
#'
#' # Components with a descriptive hub-spoke topology
#' overview[!is.na(overview$topology_class) &
#'            overview$topology_class == "hub-spoke", ]
#'
#' # Top shared-node hubs by betweenness across all components
#' comp1 <- get_component(obj, "cond_A", "comp_both", component_id = 1)
#' comp1[!is.na(comp1$betweenness_norm), ]
#'
#' # Full per-component overview is already a joined table
#' head(overview[order(-overview$n_nodes), ])
#' }
#'
#' @export
run_component_analysis <- function(obj,
                                   analysis,
                                   label,
                                   direction           = c("both", "up", "down"),
                                   level               = c("isoform", "gene"),
                                   genes               = NULL,
                                   min_component_nodes = NULL,
                                   min_shared_nodes    = NULL,
                                   min_edges           = NULL) {

  # -- Input validation --------------------------------------------------------
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.character(analysis) || length(analysis) != 1L || is.na(analysis) ||
      !nzchar(trimws(analysis)))
    stop("`analysis` must be a single non-empty character string.", call. = FALSE)

  if (!analysis %in% names(obj@analyses))
    stop(sprintf(
      "No analysis labelled '%s'. Available: %s",
      analysis, paste(names(obj@analyses), collapse = ", ")
    ), call. = FALSE)

  diff_result <- obj@analyses[[analysis]]
  if (!is(diff_result, "DifferentialPulldownAnalysis"))
    stop(sprintf(
      "Analysis '%s' is not a DifferentialPulldownAnalysis object.", analysis
    ), call. = FALSE)

  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(trimws(label)))
    stop("`label` must be a single non-empty character string.", call. = FALSE)

  if (label %in% names(diff_result@component_analyses))
    warning(sprintf(
      "A component analysis labelled '%s' already exists in '%s' and will be overwritten.",
      label, analysis
    ), call. = FALSE)

  direction <- match.arg(direction)
  level     <- match.arg(level)

  .check_int_filter <- function(x, name) {
    if (!is.null(x)) .validate_positive_integer(x, name)
  }
  .check_int_filter(min_component_nodes, "min_component_nodes")
  .check_int_filter(min_shared_nodes,    "min_shared_nodes")
  .check_int_filter(min_edges,           "min_edges")

  # -- Gene subset -------------------------------------------------------------
  all_gene_ids <- names(diff_result@gene_results)

  if (is.null(genes)) {
    genes_used <- all_gene_ids
  } else {
    bad <- setdiff(genes, all_gene_ids)
    if (length(bad))
      stop(sprintf(
        "Gene(s) not found in analysis '%s': %s\nAvailable: %s",
        analysis, paste(bad, collapse = ", "),
        paste(all_gene_ids, collapse = ", ")
      ), call. = FALSE)
    genes_used <- genes
  }

  gene_results_subset <- diff_result@gene_results[genes_used]

  # -- Build and annotate merged graph -----------------------------------------
  message(sprintf(
    "Building merged cross-gene graph for '%s' (direction = '%s', level = '%s') ...",
    analysis, direction, level
  ))

  built <- .build_cross_gene_graph(gene_results_subset, direction, level)
  g     <- built$graph

  # -- Empty-graph early exit --------------------------------------------------
  if (igraph::vcount(g) == 0L) {
    message("Merged graph is empty - no components to analyse.")
    return(.store_component_result(
      obj, analysis, label, direction, level, genes,
      .empty_component_table(), .empty_node_table(),
      .empty_topology_table(),  .empty_hub_table()
    ))
  }

  g <- .annotate_seed_direction(g, gene_results_subset, direction, level)

  # -- Decompose and filter ----------------------------------------------------
  message("Decomposing into connected components ...")
  comps <- igraph::decompose(g)
  comps <- .apply_component_filters(comps, min_component_nodes,
                                    min_shared_nodes, min_edges)

  if (length(comps) == 0L) {
    message("All components were removed by the requested filters.")
    return(.store_component_result(
      obj, analysis, label, direction, level, genes,
      .empty_component_table(), .empty_node_table(),
      .empty_topology_table(),  .empty_hub_table()
    ))
  }

  message(sprintf("%d component(s) retained.", length(comps)))

  # -- Sort components by size then edges (both descending) --------------------
  # Ensures component_id = 1 is always the largest component, making IDs
  # stable and interpretable across direction/level/gene subsets.
  # Edge count breaks ties between components with equal node counts.
  comp_sizes <- vapply(comps, igraph::vcount, double(1))
  comp_edges <- vapply(comps, igraph::ecount, double(1))
  comps      <- comps[order(comp_sizes, comp_edges, decreasing = TRUE)]

  # -- Build all four tables ---------------------------------------------------
  message("Computing Tier-1 structural statistics ...")
  component_table <- .build_component_table(comps, gene_results_subset, direction, level)
  node_table      <- .build_node_table(comps)

  message("Computing Tier-2 topological statistics ...")
  topology_table  <- .build_topology_table(comps)
  hub_table       <- .build_hub_table(comps)

  n_mixed <- sum(component_table$mixed, na.rm = TRUE)
  message(sprintf(
    "  %d node(s) total  |  %d mixed component(s)  |  %d shared hub(s) scored.",
    nrow(node_table), n_mixed, nrow(hub_table)
  ))

  # -- Store and return --------------------------------------------------------
  .store_component_result(
    obj, analysis, label, direction, level, genes,
    component_table, node_table, topology_table, hub_table
  )
}


# ==============================================================================
#  Internal: store result and return updated obj
# ==============================================================================

#' Build the ComponentAnalysis object, slot it into obj, and return obj
#' @noRd
.store_component_result <- function(obj, analysis, label, direction, level,
                                    genes, component_table, node_table,
                                    topology_table, hub_table) {
  comp_result <- new(
    "ComponentAnalysis",
    label           = as.character(label),
    parent_label    = as.character(analysis),
    component_table = component_table,
    node_table      = node_table,
    topology_table  = topology_table,
    hub_table       = hub_table,
    direction       = direction,
    level           = level,
    genes           = if (is.null(genes)) character(0) else as.character(genes),
    ora_results     = list(),
    created_at      = Sys.time()
  )
  validObject(comp_result)

  obj@analyses[[analysis]]@component_analyses[[label]] <- comp_result

  message(sprintf(
    "Done. Access the results with:\n  get_component_overview(obj, \"%s\", \"%s\")",
    analysis, label
  ))

  invisible(obj)
}


# ==============================================================================
#  Internal: annotate seed_direction vertex attribute
# ==============================================================================

#' Tag each seed vertex with its direction ("up", "down", or "both")
#' @noRd
.annotate_seed_direction <- function(g, gene_results, direction, level) {

  vnames <- igraph::V(g)$name

  if (level == "gene") {
    # At gene level, vertex names are gene IDs. A gene node is a seed if it
    # appears as gr@gene_id for a result that has isoforms in the given direction.
    up_seed_genes   <- character(0)
    down_seed_genes <- character(0)
    for (gr in gene_results) {
      if (direction %in% c("up", "both") && !is.null(gr@up_graph) &&
          length(gr@up_isoforms) > 0L)
        up_seed_genes   <- union(up_seed_genes,   gr@gene_id)
      if (direction %in% c("down", "both") && !is.null(gr@down_graph) &&
          length(gr@down_isoforms) > 0L)
        down_seed_genes <- union(down_seed_genes, gr@gene_id)
    }
    check_up   <- vnames %in% up_seed_genes
    check_down <- vnames %in% down_seed_genes
  } else {
    # At isoform level, vertex names are isoform IDs matched directly.
    up_seeds   <- character(0)
    down_seeds <- character(0)
    for (gr in gene_results) {
      if (direction %in% c("up", "both") && !is.null(gr@up_graph))
        up_seeds   <- union(up_seeds,   gr@up_isoforms)
      if (direction %in% c("down", "both") && !is.null(gr@down_graph))
        down_seeds <- union(down_seeds, gr@down_isoforms)
    }
    check_up   <- vnames %in% up_seeds
    check_down <- vnames %in% down_seeds
  }

  seed_dir <- vapply(seq_along(vnames), function(i) {
    if      (check_up[i] && check_down[i]) "both"
    else if (check_up[i])                  "up"
    else if (check_down[i])                "down"
    else                                   NA_character_
  }, character(1))

  igraph::V(g)$seed_direction <- seed_dir
  g
}


# ==============================================================================
#  Internal: component filter
# ==============================================================================

#' Drop components that do not meet size / edge / shared-node thresholds
#' @noRd
.apply_component_filters <- function(comps, min_component_nodes,
                                     min_shared_nodes, min_edges) {
  n_before   <- length(comps)
  keep_flags <- logical(n_before)

  for (i in seq_along(comps)) {
    cg <- comps[[i]]

    if (!is.null(min_component_nodes) &&
        igraph::vcount(cg) < min_component_nodes)               next
    if (!is.null(min_edges) &&
        igraph::ecount(cg) < min_edges)                         next
    if (!is.null(min_shared_nodes) &&
        sum(igraph::V(cg)$source_gene == "multi-gene",
            na.rm = TRUE) < min_shared_nodes)                   next

    keep_flags[i] <- TRUE
  }

  n_dropped <- n_before - sum(keep_flags)
  if (n_dropped > 0L)
    message(sprintf(
      "  Filters removed %d of %d component(s); %d retained.",
      n_dropped, n_before, sum(keep_flags)
    ))

  comps[keep_flags]
}


# ==============================================================================
#  Internal: Tier-1 tables
# ==============================================================================

#' Build the per-component structural summary (Tier 1)
#' @noRd
.build_component_table <- function(comps, gene_results, direction, level = c("isoform", "gene")) {

  level <- match.arg(level)

  rows <- lapply(seq_along(comps), function(i) {
    cg     <- comps[[i]]
    vattrs <- igraph::as_data_frame(cg, what = "vertices")

    n_nodes <- igraph::vcount(cg)
    n_edges <- igraph::ecount(cg)
    density <- if (n_nodes < 2L) 0 else igraph::edge_density(cg)
    avg_deg <- mean(igraph::degree(cg))
    clust   <- if (n_nodes < 2L) 0 else
      igraph::transitivity(cg, type = "average", isolates = "zero")

    is_seed   <- !is.na(vattrs$seed_direction)
    n_seeds   <- sum(is_seed)
    seed_frac <- if (n_nodes > 0L) n_seeds / n_nodes else NA_real_

    seed_dir     <- vattrs$seed_direction[is_seed]
    n_seeds_up   <- sum(seed_dir %in% c("up",   "both"))
    n_seeds_down <- sum(seed_dir %in% c("down", "both"))

    seed_vnames  <- vattrs$name[is_seed]
    seed_dirs_by_name <- stats::setNames(seed_dir, seed_vnames)
    n_genes_up   <- 0L
    n_genes_down <- 0L
    for (gr in gene_results) {
      if (level == "gene") {
        gene_seed_dir <- seed_dirs_by_name[gr@gene_id]
        if (any(gene_seed_dir %in% c("up", "both")))
          n_genes_up <- n_genes_up + 1L
        if (any(gene_seed_dir %in% c("down", "both")))
          n_genes_down <- n_genes_down + 1L
      } else {
        if (any(gr@up_isoforms   %in% seed_vnames)) n_genes_up   <- n_genes_up   + 1L
        if (any(gr@down_isoforms %in% seed_vnames)) n_genes_down <- n_genes_down + 1L
      }
    }

    n_shared    <- sum(vattrs$source_gene == "multi-gene", na.rm = TRUE)
    shared_frac <- if (n_nodes > 0L) n_shared / n_nodes else NA_real_

    data.frame(
      component_id    = i,
      n_nodes         = n_nodes,
      n_edges         = n_edges,
      density         = round(density,   4L),
      avg_degree      = round(avg_deg,   3L),
      clustering      = round(clust,     3L),
      n_seeds         = n_seeds,
      n_neighbours    = n_nodes - n_seeds,
      seed_fraction   = round(seed_frac, 3L),
      n_seeds_up      = n_seeds_up,
      n_seeds_down    = n_seeds_down,
      n_genes_up      = n_genes_up,
      n_genes_down    = n_genes_down,
      mixed           = n_seeds_up > 0L && n_seeds_down > 0L,
      n_shared_nodes  = n_shared,
      shared_fraction = round(shared_frac, 3L),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


#' Build the per-node membership table (Tier 1)
#' @noRd
.build_node_table <- function(comps) {

  rows <- lapply(seq_along(comps), function(i) {
    cg <- comps[[i]]
    data.frame(
      node           = igraph::V(cg)$name,
      component_id   = i,
      source_gene    = igraph::V(cg)$source_gene,
      is_seed        = !is.na(igraph::V(cg)$seed_direction),
      seed_direction = igraph::V(cg)$seed_direction,
      degree         = igraph::degree(cg),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


# ==============================================================================
#  Internal: Tier-2 tables
# ==============================================================================

#' Build the per-component topology table (Tier 2)
#'
#' Computes diameter, radius, and degree centralization for each component,
#' then assigns a heuristic topology class.
#'
#' Single-node components receive diameter = 0, radius = 0, and
#' degree_centralization = NA (classification is meaningless for a singleton).
#'
#' @noRd
.build_topology_table <- function(comps) {

  rows <- lapply(seq_along(comps), function(i) {
    cg <- comps[[i]]

    if (igraph::vcount(cg) < 2L) {
      diam  <- 0L
      rad   <- 0L
      centr <- NA_real_
    } else {
      diam  <- igraph::diameter(cg, directed = FALSE, unconnected = FALSE)
      rad   <- igraph::radius(cg, mode = "all")
      centr <- igraph::centr_degree(cg, mode = "all",
                                    normalized = TRUE)$centralization
    }

    # Heuristic topology classification
    tclass <- if (is.na(centr)) {
      NA_character_
    } else if (centr > 0.6) {
      "hub-spoke"
    } else if (centr < 0.2) {
      "clique-like"
    } else {
      "chain-like"
    }

    data.frame(
      component_id          = i,
      diameter              = diam,
      radius                = rad,
      degree_centralization = round(centr, 4L),
      topology_class        = tclass,
      stringsAsFactors      = FALSE
    )
  })

  do.call(rbind, rows)
}


#' Build the shared-node betweenness hub table (Tier 2)
#'
#' For every shared (multi-gene) node, computes betweenness centrality within
#' its component, normalises to \[0, 1\] within that component, and ranks.
#' The table is sorted globally by betweenness_norm descending so the most
#' critical cross-gene bridges appear at the top.
#'
#' Betweenness is computed on the full component (including non-shared nodes)
#' so that scores correctly reflect how many inter-community shortest paths
#' each shared node lies on.  Only shared nodes are then retained in the
#' output.
#'
#' @noRd
.build_hub_table <- function(comps) {

  rows <- lapply(seq_along(comps), function(i) {
    cg <- comps[[i]]

    source_gene <- igraph::V(cg)$source_gene
    is_shared   <- !is.na(source_gene) & source_gene == "multi-gene"

    if (!any(is_shared)) return(NULL)

    # Full-component betweenness (unweighted, undirected)
    btw     <- igraph::betweenness(cg, directed = FALSE, normalized = FALSE)
    max_btw <- max(btw, na.rm = TRUE)

    # Normalise within this component; guard against all-zero case
    btw_norm <- if (max_btw > 0) btw / max_btw else ifelse(is_shared, 1, 0)

    shared_idx  <- which(is_shared)
    shared_btw  <- btw[shared_idx]
    shared_norm <- btw_norm[shared_idx]
    shared_rank <- rank(-shared_btw, ties.method = "min")

    data.frame(
      node             = igraph::V(cg)$name[shared_idx],
      component_id     = i,
      degree           = igraph::degree(cg)[shared_idx],
      betweenness      = round(shared_btw,  3L),
      betweenness_norm = round(shared_norm, 4L),
      betweenness_rank = as.integer(shared_rank),
      stringsAsFactors = FALSE
    )
  })

  combined <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])

  if (is.null(combined) || nrow(combined) == 0L)
    return(.empty_hub_table())

  # Global sort: highest normalised betweenness first, ties broken by component
  combined[order(-combined$betweenness_norm, combined$component_id), ]
}


# ==============================================================================
#  Internal: empty table constructors
# ==============================================================================

#' @noRd
.empty_component_table <- function() {
  data.frame(
    component_id = integer(0), n_nodes = integer(0), n_edges = integer(0),
    density = numeric(0), avg_degree = numeric(0), clustering = numeric(0),
    n_seeds = integer(0), n_neighbours = integer(0),
    seed_fraction = numeric(0), n_seeds_up = integer(0),
    n_seeds_down = integer(0), n_genes_up = integer(0),
    n_genes_down = integer(0), mixed = logical(0),
    n_shared_nodes = integer(0), shared_fraction = numeric(0),
    stringsAsFactors = FALSE
  )
}

#' @noRd
.empty_node_table <- function() {
  data.frame(
    node = character(0), component_id = integer(0),
    source_gene = character(0), is_seed = logical(0),
    seed_direction = character(0), degree = integer(0),
    stringsAsFactors = FALSE
  )
}

#' @noRd
.empty_topology_table <- function() {
  data.frame(
    component_id = integer(0), diameter = integer(0),
    radius = numeric(0), degree_centralization = numeric(0),
    topology_class = character(0),
    stringsAsFactors = FALSE
  )
}

#' @noRd
.empty_hub_table <- function() {
  data.frame(
    node = character(0), component_id = integer(0),
    degree = integer(0), betweenness = numeric(0),
    betweenness_norm = numeric(0), betweenness_rank = integer(0),
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
#  Internal helper: resolve a ComponentAnalysis from an IsoNetAnalyzeR object
# ==============================================================================

#' Resolve and validate a ComponentAnalysis from an IsoNetAnalyzeR object
#'
#' Shared guard used by \code{get_component_overview} and
#' \code{get_component}.  Stops with an informative message if either the
#' parent \linkS4class{DifferentialPulldownAnalysis} or the
#' \linkS4class{ComponentAnalysis} label is missing.
#'
#' @param obj                An \linkS4class{IsoNetAnalyzeR} object.
#' @param analysis           Character scalar. Label of the parent
#'   \linkS4class{DifferentialPulldownAnalysis}.
#' @param component_analysis Character scalar. Label of the
#'   \linkS4class{ComponentAnalysis} stored inside the parent analysis.
#'
#' @return A \linkS4class{ComponentAnalysis} object.
#' @noRd
.resolve_component_analysis <- function(obj, analysis, component_analysis) {

  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.character(analysis) || length(analysis) != 1L ||
      nchar(trimws(analysis)) == 0L)
    stop("`analysis` must be a single non-empty character string.", call. = FALSE)

  if (!analysis %in% names(obj@analyses))
    stop(sprintf(
      "No analysis labelled '%s'. Available: %s",
      analysis, paste(names(obj@analyses), collapse = ", ")
    ), call. = FALSE)

  diff_result <- obj@analyses[[analysis]]

  if (!is(diff_result, "DifferentialPulldownAnalysis"))
    stop(sprintf(
      "Analysis '%s' is not a DifferentialPulldownAnalysis object.", analysis
    ), call. = FALSE)

  if (!is.character(component_analysis) || length(component_analysis) != 1L ||
      nchar(trimws(component_analysis)) == 0L)
    stop("`component_analysis` must be a single non-empty character string.",
         call. = FALSE)

  if (!component_analysis %in% names(diff_result@component_analyses))
    stop(sprintf(
      "No component analysis labelled '%s' in analysis '%s'. Available: %s",
      component_analysis, analysis,
      if (length(diff_result@component_analyses) == 0L) "none"
      else paste(names(diff_result@component_analyses), collapse = ", ")
    ), call. = FALSE)

  ca <- diff_result@component_analyses[[component_analysis]]

  if (!is(ca, "ComponentAnalysis"))
    stop(sprintf(
      "Component analysis '%s' is not a ComponentAnalysis object.",
      component_analysis
    ), call. = FALSE)

  ca
}


# ==============================================================================
#  Exported: get_component_overview
# ==============================================================================

#' Get a merged per-component overview table from a ComponentAnalysis result
#'
#' Returns a single wide data frame with one row per retained connected
#' component, combining all Tier-1 structural columns
#' (\code{@@component_table}) and all Tier-2 topological columns
#' (\code{@@topology_table}) side by side via a left-join on
#' \code{component_id}.
#'
#' This is the recommended starting point for exploring a
#' \linkS4class{ComponentAnalysis} result: scan the table for mixed
#' components (\code{mixed == TRUE}), hub-spoke topology
#' (\code{topology_class == "hub-spoke"}), or high shared-node fractions
#' (\code{shared_fraction}), then drill into individual components with
#' \code{\link{get_component}}.
#'
#' @param obj                An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis           Character scalar.  Label of the parent
#'   \linkS4class{DifferentialPulldownAnalysis}; use \code{list_analyses(obj)}
#'   to see available labels.
#' @param component_analysis Character scalar.  Label of the
#'   \linkS4class{ComponentAnalysis}; use \code{list_analyses(obj)} then
#'   inspect the result to see available component analysis labels.
#'
#' @return A data frame with one row per retained component and the following
#'   columns (Tier-1 first, then Tier-2):
#'   \code{component_id}, \code{n_nodes}, \code{n_edges}, \code{density},
#'   \code{avg_degree}, \code{clustering}, \code{n_seeds},
#'   \code{n_neighbours}, \code{seed_fraction}, \code{n_seeds_up},
#'   \code{n_seeds_down}, \code{n_genes_up}, \code{n_genes_down},
#'   \code{mixed}, \code{n_shared_nodes}, \code{shared_fraction},
#'   \code{diameter}, \code{radius}, \code{degree_centralization},
#'   \code{topology_class}.
#'   Topology columns are \code{NA} for any component absent from
#'   \code{@@topology_table} (e.g. single-node components where topology
#'   metrics are undefined).
#'
#' @examples
#' \dontrun{
#' obj <- IsoNetAnalyzeR("demo")
#' obj <- run_differential_pulldown(obj, "cond_A", input, order = 1)
#' obj <- run_component_analysis(obj, "cond_A", label = "comp_iso")
#'
#' overview <- get_component_overview(obj, "cond_A", "comp_iso")
#'
#' # Components containing both seed labels
#' overview[overview$mixed, ]
#'
#' # Components with a descriptive hub-spoke topology
#' overview[!is.na(overview$topology_class) &
#'            overview$topology_class == "hub-spoke", ]
#'
#' # Sort by number of nodes descending
#' overview[order(-overview$n_nodes), ]
#' }
#'
#' @export
get_component_overview <- function(obj, analysis, component_analysis) {

  ca <- .resolve_component_analysis(obj, analysis, component_analysis)

  ct <- ca@component_table
  tt <- ca@topology_table

  if (nrow(ct) == 0L) {
    message("Component table is empty - no components were retained.")
    return(ct)
  }

  # Left-join: every component row is kept; topology columns are NA when absent
  # (e.g. single-node components for which topology is undefined).
  # Base-R merge with all.x = TRUE is a left join.
  result <- merge(ct, tt, by = "component_id", all.x = TRUE,
                  sort = FALSE)

  # Restore the original component_id order (merge may reorder)
  result <- result[order(result$component_id), ]
  rownames(result) <- NULL

  result
}


# ==============================================================================
#  Exported: get_component
# ==============================================================================

#' Get a merged per-node table for a single connected component
#'
#' Returns a single data frame with one row per node in the requested
#' component, combining the Tier-1 node membership columns
#' (\code{@@node_table}) with the Tier-2 hub betweenness columns
#' (\code{@@hub_table}) via a left-join on \code{node} and
#' \code{component_id}.
#'
#' Nodes that are not shared across genes (\code{source_gene != "multi-gene"})
#' receive \code{NA} in the hub columns (\code{betweenness},
#' \code{betweenness_norm}, \code{betweenness_rank}), as betweenness is only
#' computed for shared nodes.  This makes it easy to filter for hub nodes with
#' \code{!is.na(betweenness_norm)} or to sort all nodes by their hub score.
#'
#' @param obj                An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis           Character scalar.  Label of the parent
#'   \linkS4class{DifferentialPulldownAnalysis}; use \code{list_analyses(obj)}
#'   to see available labels.
#' @param component_analysis Character scalar.  Label of the
#'   \linkS4class{ComponentAnalysis}; use \code{list_analyses(obj)} then
#'   inspect the result to see available component analysis labels.
#' @param component_id       Integer scalar.  The component to retrieve, as
#'   listed in \code{component_id} column of \code{\link{get_component_overview}}.
#'   Component IDs are assigned in descending order of node count (edge count
#'   as tiebreaker), so \code{component_id = 1} is always the largest
#'   component.
#'
#' @return A data frame with one row per node in the requested component and
#'   the following columns:
#'   \code{node}, \code{component_id}, \code{source_gene}, \code{is_seed},
#'   \code{seed_direction}, \code{degree},
#'   \code{betweenness} (\code{NA} for non-shared nodes),
#'   \code{betweenness_norm} (\code{NA} for non-shared nodes),
#'   \code{betweenness_rank} (\code{NA} for non-shared nodes).
#'   Rows are sorted by \code{betweenness_norm} descending (hub nodes first),
#'   with \code{NA}s at the end.
#'
#' @examples
#' \dontrun{
#' obj <- IsoNetAnalyzeR("demo")
#' obj <- run_differential_pulldown(obj, "cond_A", input, order = 1)
#' obj <- run_component_analysis(obj, "cond_A", label = "comp_iso")
#'
#' # Inspect the largest component
#' comp1 <- get_component(obj, "cond_A", "comp_iso", component_id = 1)
#' comp1
#'
#' # Seed nodes only
#' comp1[comp1$is_seed, ]
#'
#' # Shared hub nodes only (have a betweenness score)
#' comp1[!is.na(comp1$betweenness_norm), ]
#'
#' # Top hub by betweenness within this component
#' comp1[which.min(comp1$betweenness_rank), ]
#' }
#'
#' @export
get_component <- function(obj, analysis, component_analysis, component_id) {

  ca <- .resolve_component_analysis(obj, analysis, component_analysis)

  # -- Validate component_id --------------------------------------------------
  component_id <- .validate_positive_integer(component_id, "component_id")

  valid_ids <- ca@component_table$component_id
  if (!component_id %in% valid_ids)
    stop(sprintf(
      "No component with id %d. Available component_id values: %s",
      component_id,
      paste(sort(valid_ids), collapse = ", ")
    ), call. = FALSE)

  # -- Subset node_table to this component ------------------------------------
  nt <- ca@node_table[ca@node_table$component_id == component_id, ,
                      drop = FALSE]

  if (nrow(nt) == 0L) {
    message(sprintf("No nodes found for component_id %d.", component_id))
    return(nt)
  }

  # -- Subset hub_table to this component ------------------------------------
  ht <- ca@hub_table[ca@hub_table$component_id == component_id, ,
                     drop = FALSE]

  # Keep only the hub-specific columns for the join (drop degree - already in nt)
  hub_cols <- c("node", "component_id", "betweenness", "betweenness_norm",
                "betweenness_rank")
  ht <- ht[, hub_cols[hub_cols %in% names(ht)], drop = FALSE]

  # -- Left-join: all nodes kept, hub columns NA for non-shared nodes ---------
  result <- merge(nt, ht, by = c("node", "component_id"), all.x = TRUE,
                  sort = FALSE)

  # Sort: hub nodes (non-NA betweenness_norm) first, ranked by betweenness_norm
  # descending; non-hub nodes follow sorted by degree descending.
  is_hub   <- !is.na(result$betweenness_norm)
  hub_rows <- result[is_hub,  , drop = FALSE]
  other    <- result[!is_hub, , drop = FALSE]

  if (nrow(hub_rows) > 0L)
    hub_rows <- hub_rows[order(-hub_rows$betweenness_norm), , drop = FALSE]
  if (nrow(other) > 0L)
    other <- other[order(-other$degree), , drop = FALSE]

  result <- rbind(hub_rows, other)
  rownames(result) <- NULL

  result
}
