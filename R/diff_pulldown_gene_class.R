## diff_pulldown_gene_class.R
## S4 class for per-gene differential virtual pulldown results.
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: GeneDiffPulldown class, plot_network methods


# ==============================================================================
#  Class definition
# ==============================================================================

#' GeneDiffPulldown - per-gene differential pulldown result
#'
#' Stores the virtual pulldown subgraphs and statistics for a single gene's
#' up- and down-regulated isoforms.  One \linkS4class{GeneDiffPulldown} object
#' is created per gene in a \linkS4class{DifferentialPulldownAnalysis} run.
#'
#' @slot gene_id       Character scalar. The gene identifier from the input
#'   table.
#' @slot up_isoforms   Character vector. Isoform IDs from the input table that
#'   were marked as up-regulated for this gene.
#' @slot down_isoforms Character vector. Isoform IDs from the input table that
#'   were marked as down-regulated for this gene.
#' @slot up_graph      An \code{igraph} object. Virtual pulldown subgraph seeded
#'   from \code{up_isoforms}, or \code{NULL} when no up-regulated isoforms were
#'   present in the network.
#' @slot down_graph    An \code{igraph} object. Virtual pulldown subgraph seeded
#'   from \code{down_isoforms}, or \code{NULL} when no down-regulated isoforms
#'   were present in the network.
#' @slot up_stats      Data frame. Per-node statistics for \code{up_graph}
#'   (output of \code{\link{analyze_subgraph}}), or a zero-row data frame when
#'   \code{up_graph} is \code{NULL}.
#' @slot down_stats    Data frame. Per-node statistics for \code{down_graph},
#'   or a zero-row data frame when \code{down_graph} is \code{NULL}.
#' @slot graph_summary Data frame with one row per direction (\code{"up"},
#'   \code{"down"}) summarising network-level properties: nodes, edges,
#'   density, avg_degree, clustering, starness, genes.  Rows are \code{NA}
#'   when the corresponding graph is absent.
#' @slot gene_summary  Named list of gene-level overview scalars:
#'   \code{n_up_input}, \code{n_down_input} (isoforms supplied in the table),
#'   \code{n_up_network}, \code{n_down_network} (isoforms found in the network),
#'   \code{has_up}, \code{has_down} (logical flags for subgraph availability).
#'
#' @exportClass GeneDiffPulldown
setClass(
  "GeneDiffPulldown",
  representation(
    gene_id       = "character",
    up_isoforms   = "character",
    down_isoforms = "character",
    up_graph      = "ANY",          # igraph or NULL
    down_graph    = "ANY",          # igraph or NULL
    up_stats      = "data.frame",
    down_stats    = "data.frame",
    graph_summary = "data.frame",
    gene_summary  = "list"
  )
)


# ==============================================================================
#  Validator
# ==============================================================================

setValidity("GeneDiffPulldown", function(object) {
  msgs <- character()

  if (length(object@gene_id) != 1L || nchar(trimws(object@gene_id)) == 0L)
    msgs <- c(msgs, "@gene_id must be a non-empty character scalar.")

  if (!is.character(object@up_isoforms))
    msgs <- c(msgs, "@up_isoforms must be a character vector.")

  if (!is.character(object@down_isoforms))
    msgs <- c(msgs, "@down_isoforms must be a character vector.")

  # up_graph and down_graph may be NULL (direction absent) or igraph
  if (!is.null(object@up_graph) && !inherits(object@up_graph, "igraph"))
    msgs <- c(msgs, "@up_graph must be an igraph object or NULL.")

  if (!is.null(object@down_graph) && !inherits(object@down_graph, "igraph"))
    msgs <- c(msgs, "@down_graph must be an igraph object or NULL.")

  if (length(msgs)) msgs else TRUE
})


# ==============================================================================
#  Internal constructor
# ==============================================================================

#' Build a GeneDiffPulldown object
#'
#' Internal constructor called by \code{\link{run_differential_pulldown}}.
#' Assembles all slots from pre-computed components.
#'
#' @param gene_id          Character scalar.
#' @param up_isoforms      Character vector of input up-regulated isoform IDs.
#' @param down_isoforms    Character vector of input down-regulated isoform IDs.
#' @param up_graph         igraph or NULL.
#' @param down_graph       igraph or NULL.
#' @param up_stats         Data frame from \code{analyze_subgraph}, or zero-row.
#' @param down_stats       Data frame from \code{analyze_subgraph}, or zero-row.
#' @param network_vertices Character vector of all vertex names in the session
#'   network.  Forwarded to \code{.gene_overview()} to compute coverage counts.
#'
#' @return A \linkS4class{GeneDiffPulldown} object.
#' @noRd
.new_gene_diff_pulldown <- function(gene_id, up_isoforms, down_isoforms,
                                    up_graph, down_graph,
                                    up_stats, down_stats,
                                    network_vertices = character(0)) {

  graph_summary <- .summarise_direction_pair(up_graph, down_graph)
  gene_summary  <- .gene_overview(up_isoforms, down_isoforms,
                                  up_graph, down_graph,
                                  network_vertices = network_vertices)

  obj <- new(
    "GeneDiffPulldown",
    gene_id       = as.character(gene_id),
    up_isoforms   = as.character(up_isoforms),
    down_isoforms = as.character(down_isoforms),
    up_graph      = up_graph,
    down_graph    = down_graph,
    up_stats      = up_stats,
    down_stats    = down_stats,
    graph_summary = graph_summary,
    gene_summary  = gene_summary
  )

  validObject(obj)
  obj
}


# ==============================================================================
#  Internal helpers for slot construction
# ==============================================================================

#' Build the graph_summary data frame for one gene
#'
#' Returns a two-row data frame (up / down).  Rows for absent graphs are NA.
#'
#' @param up_graph   igraph or NULL.
#' @param down_graph igraph or NULL.
#' @return data.frame with columns: direction, nodes, edges, density,
#'   avg_degree, clustering, starness, genes.
#' @noRd
.summarise_direction_pair <- function(up_graph, down_graph) {

  .stats_row <- function(g, direction) {
    if (is.null(g)) {
      return(data.frame(
        direction  = direction,
        nodes      = NA_integer_,
        edges      = NA_integer_,
        density    = NA_real_,
        avg_degree = NA_real_,
        clustering = NA_real_,
        starness   = NA_real_,
        genes      = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      direction  = direction,
      nodes      = igraph::vcount(g),
      edges      = igraph::ecount(g),
      density    = igraph::edge_density(g),
      avg_degree = mean(igraph::degree(g)),
      clustering = igraph::transitivity(g, type = "average"),
      starness   = igraph::centr_degree(g, normalized = TRUE)$centralization,
      genes      = length(unique(igraph::V(g)$gene)),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    .stats_row(up_graph,   "up"),
    .stats_row(down_graph, "down")
  )
}


#' Build the gene_summary list for one gene
#'
#' @param up_isoforms     Character vector from input table.
#' @param down_isoforms   Character vector from input table.
#' @param up_graph        igraph or NULL.
#' @param down_graph      igraph or NULL.
#' @param network_vertices Character vector of all vertex names in the session
#'   network.  Used to count how many input isoforms are actually present in
#'   the network before BFS expansion.
#' @return Named list.
#' @noRd
.gene_overview <- function(up_isoforms, down_isoforms, up_graph, down_graph,
                           network_vertices = character(0)) {
  list(
    n_up_input      = length(up_isoforms),
    n_down_input    = length(down_isoforms),
    n_up_found      = sum(up_isoforms   %in% network_vertices),
    n_down_found    = sum(down_isoforms %in% network_vertices),
    n_up_network    = if (!is.null(up_graph))   igraph::vcount(up_graph)   else 0L,
    n_down_network  = if (!is.null(down_graph)) igraph::vcount(down_graph) else 0L,
    has_up          = !is.null(up_graph),
    has_down        = !is.null(down_graph)
  )
}


# ==============================================================================
#  show()
# ==============================================================================

#' @describeIn GeneDiffPulldown Print a concise per-gene summary.
#' @param object A \linkS4class{GeneDiffPulldown} object.
setMethod("show", "GeneDiffPulldown", function(object) {
  gs <- object@gene_summary

  cli::cli_rule("GeneDiffPulldown")
  cli::cli_text("Gene       : {object@gene_id}")

  up_desc <- if (gs$has_up) {
    paste0("subgraph with ",
           object@graph_summary$nodes[object@graph_summary$direction == "up"],
           " nodes, ",
           object@graph_summary$edges[object@graph_summary$direction == "up"],
           " edges")
  } else {
    "no network hits"
  }

  down_desc <- if (gs$has_down) {
    paste0("subgraph with ",
           object@graph_summary$nodes[object@graph_summary$direction == "down"],
           " nodes, ",
           object@graph_summary$edges[object@graph_summary$direction == "down"],
           " edges")
  } else {
    "no network hits"
  }

  cli::cli_text("Up         : {gs$n_up_input} isoform(s) input
               {cli::symbol$arrow_right} {up_desc}")
  cli::cli_text("Down       : {gs$n_down_input} isoform(s) input
               {cli::symbol$arrow_right} {down_desc}")
  cli::cli_rule()
  invisible(object)
})
