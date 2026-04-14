## diff_pulldown_analysis_class.R
## S4 class for the full differential pulldown analysis result.
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: DifferentialPulldownAnalysis class, show method


# ==============================================================================
#  Class definition
# ==============================================================================

#' DifferentialPulldownAnalysis — full differential pulldown analysis result
#'
#' Top-level container produced by \code{\link{run_differential_pulldown}}.
#' Stores one \linkS4class{GeneDiffPulldown} object per gene, the validated
#' input table, the run parameters, and cross-gene summary statistics.
#'
#' @slot label        Character scalar. A user-supplied label for this
#'   analysis run, used as the key in \code{obj@@analyses}.
#' @slot input_table  Data frame. The validated long-format input table with
#'   columns \code{gene_id}, \code{isoform_id}, \code{direction}
#'   (\code{"up"} / \code{"down"}).
#' @slot gene_results Named list of \linkS4class{GeneDiffPulldown} objects,
#'   one per gene in \code{input_table}.  Names are gene IDs.
#' @slot cross_gene_summary Data frame with one row per gene summarising
#'   up/down subgraph properties side by side.  Columns:
#'   \code{gene_id}, \code{n_up_input}, \code{n_down_input},
#'   \code{n_up_network}, \code{n_down_network},
#'   \code{up_nodes}, \code{up_edges}, \code{up_density},
#'   \code{up_avg_degree}, \code{up_clustering},
#'   \code{down_nodes}, \code{down_edges}, \code{down_density},
#'   \code{down_avg_degree}, \code{down_clustering}.
#' @slot order        Integer scalar. Number of BFS hops used.
#' @slot cutoff       Numeric scalar. Internal-edge fraction cutoff used.
#' @slot created_at   \code{POSIXct} scalar. Timestamp of the run.
#'
#' @exportClass DifferentialPulldownAnalysis
setClass(
  "DifferentialPulldownAnalysis",
  representation(
    label               = "character",
    input_table         = "data.frame",
    gene_results        = "list",
    cross_gene_summary  = "data.frame",
    order               = "integer",
    cutoff              = "numeric",
    created_at          = "POSIXct"
  )
)


# ==============================================================================
#  Validator
# ==============================================================================

setValidity("DifferentialPulldownAnalysis", function(object) {
  msgs <- character()

  if (length(object@label) != 1L || nchar(trimws(object@label)) == 0L)
    msgs <- c(msgs, "@label must be a non-empty character scalar.")

  required_cols <- c("gene_id", "isoform_id", "direction")
  missing_cols  <- setdiff(required_cols, names(object@input_table))
  if (length(missing_cols))
    msgs <- c(msgs, paste0(
      "@input_table is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    ))

  bad_genes <- vapply(
    object@gene_results,
    function(x) !is(x, "GeneDiffPulldown"),
    logical(1)
  )
  if (any(bad_genes))
    msgs <- c(msgs, paste0(
      "@gene_results contains non-GeneDiffPulldown element(s): ",
      paste(names(object@gene_results)[bad_genes], collapse = ", ")
    ))

  if (length(object@order) != 1L || is.na(object@order) || object@order < 1L)
    msgs <- c(msgs, "@order must be a positive integer scalar.")

  if (length(object@cutoff) != 1L || is.na(object@cutoff) ||
      object@cutoff < 0 || object@cutoff > 1)
    msgs <- c(msgs, "@cutoff must be a single numeric value in [0, 1].")

  if (length(object@created_at) != 1L)
    msgs <- c(msgs, "@created_at must be a length-1 POSIXct.")

  if (length(msgs)) msgs else TRUE
})


# ==============================================================================
#  Internal constructor
# ==============================================================================

#' Build a DifferentialPulldownAnalysis object
#'
#' Internal constructor called by \code{\link{run_differential_pulldown}}.
#'
#' @param label        Character scalar.
#' @param input_table  Validated data frame.
#' @param gene_results Named list of GeneDiffPulldown objects.
#' @param order        Integer scalar.
#' @param cutoff       Numeric scalar.
#'
#' @return A \linkS4class{DifferentialPulldownAnalysis} object.
#' @noRd
.new_diff_pulldown_analysis <- function(label, input_table,
                                        gene_results, order, cutoff) {

  cross_gene_summary <- .build_cross_gene_summary(gene_results)

  obj <- new(
    "DifferentialPulldownAnalysis",
    label              = as.character(label),
    input_table        = input_table,
    gene_results       = gene_results,
    cross_gene_summary = cross_gene_summary,
    order              = as.integer(order),
    cutoff             = as.numeric(cutoff),
    created_at         = Sys.time()
  )

  validObject(obj)
  obj
}


# ==============================================================================
#  Internal: cross-gene summary builder
# ==============================================================================

#' Build the cross-gene summary data frame
#'
#' Iterates over all per-gene results and assembles a single wide data frame
#' with one row per gene.
#'
#' @param gene_results Named list of GeneDiffPulldown objects.
#' @return data.frame.
#' @noRd
.build_cross_gene_summary <- function(gene_results) {

  do.call(rbind, lapply(gene_results, function(gr) {
    gs  <- gr@gene_summary
    gdf <- gr@graph_summary

    up_row   <- gdf[gdf$direction == "up",   , drop = FALSE]
    down_row <- gdf[gdf$direction == "down", , drop = FALSE]

    data.frame(
      gene_id          = gr@gene_id,
      n_up_input       = gs$n_up_input,
      n_down_input     = gs$n_down_input,
      n_up_network     = gs$n_up_network,
      n_down_network   = gs$n_down_network,
      up_nodes         = up_row$nodes,
      up_edges         = up_row$edges,
      up_density       = up_row$density,
      up_avg_degree    = up_row$avg_degree,
      up_clustering    = up_row$clustering,
      down_nodes       = down_row$nodes,
      down_edges       = down_row$edges,
      down_density     = down_row$density,
      down_avg_degree  = down_row$avg_degree,
      down_clustering  = down_row$clustering,
      stringsAsFactors = FALSE
    )
  }))
}


# ==============================================================================
#  show()
# ==============================================================================

#' @describeIn DifferentialPulldownAnalysis Print a concise run summary.
setMethod("show", "DifferentialPulldownAnalysis", function(object) {

  n_genes    <- length(object@gene_results)
  n_up_ok    <- sum(vapply(object@gene_results,
                           function(x) x@gene_summary$has_up,  logical(1)))
  n_down_ok  <- sum(vapply(object@gene_results,
                           function(x) x@gene_summary$has_down, logical(1)))

  cat("── DifferentialPulldownAnalysis ─────────────────────────────────────────\n")
  cat(sprintf("  Label      : %s\n", object@label))
  cat(sprintf("  Created    : %s\n", format(object@created_at, "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("  Order      : %d    Cutoff : %.4f\n", object@order, object@cutoff))
  cat(sprintf("  Genes      : %d total\n", n_genes))
  cat(sprintf("  Up graphs  : %d / %d genes produced a subgraph\n", n_up_ok,   n_genes))
  cat(sprintf("  Down graphs: %d / %d genes produced a subgraph\n", n_down_ok, n_genes))
  cat(sprintf("  Input rows : %d isoform entries\n", nrow(object@input_table)))
  cat("─────────────────────────────────────────────────────────────────────────\n")
  invisible(object)
})
