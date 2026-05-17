## diff_pulldown_analysis_class.R
## S4 class for the full differential pulldown analysis result.
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: DifferentialPulldownAnalysis class, show method

#' @include diff_pulldown_gene_class.R
#' @include isonetanalyzer_class.R
NULL

# ==============================================================================
#  Class definition
# ==============================================================================

#' DifferentialPulldownAnalysis - full differential pulldown analysis result
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
#'   \code{gene_id},
#'   \code{n_up_input}, \code{n_down_input} (isoforms in the input table),
#'   \code{n_up_found}, \code{n_down_found} (input isoforms present in the
#'     network before BFS expansion),
#'   \code{pct_up_found}, \code{pct_down_found} (percentage of input isoforms
#'     found in the network; \code{NA} when \code{n_input == 0}),
#'   \code{up_nodes}, \code{up_edges}, \code{up_density},
#'   \code{up_avg_degree}, \code{up_clustering},
#'   \code{up_starness} (degree centralization of the up subgraph; 0 = uniform,
#'     1 = perfect hub-and-spoke topology),
#'   \code{up_genes} (number of unique genes represented in the up subgraph),
#'   \code{down_nodes}, \code{down_edges}, \code{down_density},
#'   \code{down_avg_degree}, \code{down_clustering},
#'   \code{down_starness} (degree centralization of the down subgraph),
#'   \code{down_genes} (number of unique genes represented in the down subgraph),
#'   \code{n_shared_nodes} (nodes present in both up and down subgraphs),
#'   \code{n_up_unique} (nodes exclusive to the up subgraph),
#'   \code{n_down_unique} (nodes exclusive to the down subgraph),
#'   \code{jaccard} (Jaccard index of node-set overlap between up and down
#'     subgraphs: \code{|up ∩ down| / |up ∪ down|}; \code{NA} when either
#'     subgraph is absent),
#'   \code{n_edges_gained} (edges present in the up subgraph but absent from
#'     the down subgraph — interactions gained when switching to the up
#'     isoform),
#'   \code{n_edges_lost} (edges present in the down subgraph but absent from
#'     the up subgraph — interactions lost when switching to the up isoform),
#'   \code{n_edges_shared} (edges present in both subgraphs; together with
#'     \code{n_edges_gained} and \code{n_edges_lost} this gives a complete
#'     edge-level accounting of rewiring).  All three edge columns are
#'     \code{NA} when either subgraph is absent,
#'   \code{up_seed_mean_pi} (mean of per-seed mean-pi values across all up
#'     seeds; summarises average interaction confidence of the up seed set),
#'   \code{up_seed_pi_entropy} (mean of per-seed pi-entropy values across all
#'     up seeds; summarises how evenly interaction confidence is distributed
#'     across partners for the up seed set),
#'   \code{down_seed_mean_pi}, \code{down_seed_pi_entropy} (analogous metrics
#'     for the down seed set).  All four columns are \code{NA} when the
#'     corresponding subgraph is absent or carries no \code{pi} edge attribute.
#' @slot gsea_results  Named list of \linkS4class{FgseaResult} objects,
#'   keyed by user-supplied label.  Populated by \code{\link{run_fgsea}}.
#'   Empty list by default.
#' @slot component_analyses Named list of \linkS4class{ComponentAnalysis}
#'   objects, keyed by user-supplied label.  Populated by
#'   \code{\link{run_component_analysis}}.  Empty list by default.
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
    gsea_results        = "list",
    component_analyses  = "list",
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
    gsea_results       = list(),
    component_analyses = list(),
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
#' Note: \code{n_up_network} / \code{n_down_network} (formerly in
#' \code{gene_summary}) are identical to \code{up_nodes} / \code{down_nodes}
#' from \code{graph_summary} — both call \code{igraph::vcount()} on the same
#' graph — and have been removed to avoid duplication.
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

    # Percentage of input isoforms that were found in the network (pre-BFS).
    # NA when a direction had no input isoforms to avoid a 0/0 division.
    pct_up_found   <- if (gs$n_up_input   > 0L)
      round(100 * gs$n_up_found   / gs$n_up_input,   1L) else NA_real_
    pct_down_found <- if (gs$n_down_input > 0L)
      round(100 * gs$n_down_found / gs$n_down_input, 1L) else NA_real_

    # ---- Jaccard node-set overlap between up and down subgraphs -------------
    # Computed from the igraph objects directly so the three counts are always
    # internally consistent (shared + up_unique + down_unique == union size).
    has_up   <- !is.null(gr@up_graph)
    has_down <- !is.null(gr@down_graph)

    if (has_up && has_down) {
      up_names   <- igraph::V(gr@up_graph)$name
      down_names <- igraph::V(gr@down_graph)$name
      n_shared      <- length(intersect(up_names, down_names))
      n_up_unique   <- length(setdiff(up_names,   down_names))
      n_down_unique <- length(setdiff(down_names, up_names))
      union_size    <- n_shared + n_up_unique + n_down_unique
      jaccard       <- if (union_size > 0L)
        round(n_shared / union_size, 4L) else NA_real_
    } else {
      # Cannot compute overlap when one or both subgraphs are absent.
      n_shared      <- NA_integer_
      n_up_unique   <- NA_integer_
      n_down_unique <- NA_integer_
      jaccard       <- NA_real_
    }

    # ---- Edge-level gain / loss between up and down subgraphs ---------------
    # An edge is "gained" in the up condition if it is present in up_graph but
    # absent from down_graph, and "lost" if present in down_graph but absent
    # from up_graph.  Because the graphs are undirected, endpoints are sorted
    # before hashing so that A--B and B--A collapse to the same key.
    .edge_keys <- function(g) {
      if (is.null(g) || igraph::ecount(g) == 0L) return(character(0))
      edf <- igraph::as_data_frame(g, what = "edges")
      paste(pmin(edf$from, edf$to), pmax(edf$from, edf$to), sep = ":::")
    }

    if (has_up && has_down) {
      up_keys   <- .edge_keys(gr@up_graph)
      down_keys <- .edge_keys(gr@down_graph)
      n_edges_gained  <- length(setdiff(up_keys,   down_keys))
      n_edges_lost    <- length(setdiff(down_keys, up_keys))
      n_edges_shared  <- length(intersect(up_keys, down_keys))
    } else {
      n_edges_gained  <- NA_integer_
      n_edges_lost    <- NA_integer_
      n_edges_shared  <- NA_integer_
    }

    # ---- Aggregated seed pi stats -------------------------------------------
    # up_stats / down_stats carry mean_pi and pi_entropy for seed rows only
    # (NA for neighbour nodes). Here we average across seeds per direction so
    # the cross-gene table captures the overall confidence quality of each
    # gene's seed set. Averaging per-seed entropies (rather than pooling all
    # pi values) keeps each seed's distribution intact.
    .agg_seed_pi <- function(stats_df) {
      if (nrow(stats_df) == 0L)
        return(c(mean_pi = NA_real_, pi_entropy = NA_real_))
      seed_rows <- stats_df[!is.na(stats_df$mean_pi), , drop = FALSE]
      if (nrow(seed_rows) == 0L)
        return(c(mean_pi = NA_real_, pi_entropy = NA_real_))
      c(
        mean_pi    = round(mean(seed_rows$mean_pi,    na.rm = TRUE), 4L),
        pi_entropy = round(mean(seed_rows$pi_entropy, na.rm = TRUE), 4L)
      )
    }

    up_pi   <- .agg_seed_pi(gr@up_stats)
    down_pi <- .agg_seed_pi(gr@down_stats)

    data.frame(
      gene_id          = gr@gene_id,
      # ---- input coverage --------------------------------------------------
      n_up_input       = gs$n_up_input,
      n_down_input     = gs$n_down_input,
      n_up_found       = gs$n_up_found,
      n_down_found     = gs$n_down_found,
      pct_up_found     = pct_up_found,
      pct_down_found   = pct_down_found,
      # ---- subgraph topology (up_nodes / down_nodes are the authoritative
      #      node counts; n_up/down_network was an identical duplicate) -------
      up_nodes         = up_row$nodes,
      up_edges         = up_row$edges,
      up_density       = up_row$density,
      up_avg_degree    = up_row$avg_degree,
      up_clustering    = up_row$clustering,
      up_starness      = up_row$starness,
      up_genes         = up_row$genes,
      down_nodes       = down_row$nodes,
      down_edges       = down_row$edges,
      down_density     = down_row$density,
      down_avg_degree  = down_row$avg_degree,
      down_clustering  = down_row$clustering,
      down_starness    = down_row$starness,
      down_genes       = down_row$genes,
      # ---- node-set overlap (Jaccard) --------------------------------------
      n_shared_nodes   = n_shared,
      n_up_unique      = n_up_unique,
      n_down_unique    = n_down_unique,
      jaccard          = jaccard,
      # ---- edge-level gain / loss ------------------------------------------
      n_edges_gained   = n_edges_gained,
      n_edges_lost     = n_edges_lost,
      n_edges_shared   = n_edges_shared,
      # ---- aggregated seed pi stats ----------------------------------------
      up_seed_mean_pi      = up_pi[["mean_pi"]],
      up_seed_pi_entropy   = up_pi[["pi_entropy"]],
      down_seed_mean_pi    = down_pi[["mean_pi"]],
      down_seed_pi_entropy = down_pi[["pi_entropy"]],
      stringsAsFactors = FALSE
    )
  }))
}


# ==============================================================================
#  show()
# ==============================================================================

#' @describeIn DifferentialPulldownAnalysis Print a concise run summary.
#' @param object A \linkS4class{DifferentialPulldownAnalysis} object.
setMethod("show", "DifferentialPulldownAnalysis", function(object) {

  n_genes    <- length(object@gene_results)
  n_up_ok    <- sum(vapply(object@gene_results,
                           function(x) x@gene_summary$has_up,  logical(1)))
  n_down_ok  <- sum(vapply(object@gene_results,
                           function(x) x@gene_summary$has_down, logical(1)))

  cli::cli_rule("DifferentialPulldownAnalysis")
  cli::cli_text("Label      : {object@label}")
  cli::cli_text("Created    : {format(object@created_at, '%Y-%m-%d %H:%M:%S')}")
  cli::cli_text("Order      : {object@order}    Cutoff : {round(object@cutoff, 4)}")
  cli::cli_text("Genes      : {n_genes} total")
  cli::cli_text("Up graphs  : {n_up_ok} / {n_genes} genes produced a subgraph")
  cli::cli_text("Down graphs: {n_down_ok} / {n_genes} genes produced a subgraph")
  cli::cli_text("Input rows : {nrow(object@input_table)} isoform entries")
  cli::cli_rule()
  invisible(object)
})


# ==============================================================================
#  Internal helper: resolve a DifferentialPulldownAnalysis from IsoNetAnalyzeR
# ==============================================================================

#' Resolve and validate a DifferentialPulldownAnalysis from an IsoNetAnalyzeR
#'
#' Shared guard used by all IsoNetAnalyzeR-level differential pulldown
#' accessors.  Stops with an informative message if the label is missing or
#' points to the wrong class.
#'
#' @param obj      An \linkS4class{IsoNetAnalyzeR} object.
#' @param analysis Character scalar. Label of the analysis to retrieve.
#' @return A \linkS4class{DifferentialPulldownAnalysis} object.
#' @noRd
.resolve_diff_analysis <- function(obj, analysis) {
  if (!analysis %in% names(obj@analyses))
    stop(sprintf(
      "No analysis labelled '%s'. Available: %s",
      analysis, paste(names(obj@analyses), collapse = ", ")
    ), call. = FALSE)

  result <- obj@analyses[[analysis]]

  if (!is(result, "DifferentialPulldownAnalysis"))
    stop(sprintf(
      "Analysis '%s' is not a DifferentialPulldownAnalysis object.", analysis
    ), call. = FALSE)

  result
}


# ==============================================================================
#  get_cross_gene_summary
# ==============================================================================

#' Get the cross-gene summary table from a differential pulldown analysis
#'
#' Returns the wide-format summary data frame with one row per gene,
#' containing subgraph topology metrics for both the up- and down-regulated
#' conditions side by side.
#'
#' @param obj      An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis Character scalar. Label of the
#'   \linkS4class{DifferentialPulldownAnalysis} stored in \code{obj@@analyses}.
#'
#' @return A data frame with one row per gene. See
#'   \linkS4class{DifferentialPulldownAnalysis} for a full description of
#'   columns.
#'
#' @examples
#' \dontrun{
#' cgs <- get_cross_gene_summary(obj, "my_run")
#' head(cgs)
#' }
#'
#' @export
setGeneric("get_cross_gene_summary",
           function(obj, analysis) standardGeneric("get_cross_gene_summary"))

#' @describeIn get_cross_gene_summary Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_cross_gene_summary", "IsoNetAnalyzeR",
          function(obj, analysis) {
            .resolve_diff_analysis(obj, analysis)@cross_gene_summary
          })


# ==============================================================================
#  get_gene_result
# ==============================================================================

#' Get the per-gene result object from a differential pulldown analysis
#'
#' Returns the \linkS4class{GeneDiffPulldown} object for a single gene,
#' which contains the up and down subgraphs, per-node statistics, and
#' graph-level summaries for that gene.
#'
#' @param obj      An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis Character scalar. Label of the
#'   \linkS4class{DifferentialPulldownAnalysis} stored in \code{obj@@analyses}.
#' @param gene_id  Character scalar. Ensembl gene ID (or whichever identifier
#'   was used in the input table).
#'
#' @return A \linkS4class{GeneDiffPulldown} object.
#'
#' @examples
#' \dontrun{
#' gr <- get_gene_result(obj, "my_run", "ENSG00000120333")
#' gr@@up_stats
#' gr@@down_stats
#' }
#'
#' @export
setGeneric("get_gene_result",
           function(obj, analysis, gene_id) standardGeneric("get_gene_result"))

#' @describeIn get_gene_result Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_gene_result", "IsoNetAnalyzeR",
          function(obj, analysis, gene_id) {
            diff_result <- .resolve_diff_analysis(obj, analysis)

            if (!gene_id %in% names(diff_result@gene_results))
              stop(sprintf(
                "No result for gene '%s' in analysis '%s'. Available: %s",
                gene_id, analysis,
                paste(names(diff_result@gene_results), collapse = ", ")
              ), call. = FALSE)

            diff_result@gene_results[[gene_id]]
          })


# ==============================================================================
#  get_isoform_stats
# ==============================================================================

#' Get combined isoform-level statistics for a gene from a differential
#' pulldown analysis
#'
#' Merges the per-node statistics from the up- and down-regulated subgraphs of
#' a single gene into one long-format data frame, with a \code{direction}
#' column (\code{"up"} / \code{"down"}) identifying the condition.  A
#' \code{is_seed} column flags isoforms that were supplied as seeds (i.e.
#' distance from seed == 0).  Isoforms that appear in only one condition are
#' represented in that condition's rows only; no artificial NAs are introduced
#' for the missing direction.
#'
#' @param obj      An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis Character scalar. Label of the
#'   \linkS4class{DifferentialPulldownAnalysis} stored in \code{obj@@analyses}.
#' @param gene_id  Character scalar. Ensembl gene ID (or whichever identifier
#'   was used in the input table).
#'
#' @return A data frame with one row per isoform per direction and columns:
#'   \code{direction}, \code{isoform}, \code{is_seed}, \code{degree},
#'   \code{clustering}, \code{distance_from_seed}, \code{mean_pi},
#'   \code{pi_entropy}.
#'
#' @examples
#' \dontrun{
#' stats <- get_isoform_stats(obj, "my_run", "ENSG00000120333")
#' # All isoforms in the up condition
#' stats[stats$direction == "up", ]
#' # Seed isoforms only
#' stats[stats$is_seed, ]
#' }
#'
#' @export
setGeneric("get_isoform_stats",
           function(obj, analysis, gene_id) standardGeneric("get_isoform_stats"))

#' @describeIn get_isoform_stats Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_isoform_stats", "IsoNetAnalyzeR",
          function(obj, analysis, gene_id) {
            gr <- get_gene_result(obj, analysis, gene_id)

            .prep <- function(stats_df, direction) {
              if (nrow(stats_df) == 0L) return(NULL)
              df            <- stats_df
              df$direction  <- direction
              df$is_seed    <- df$distance_from_seed == 0
              # Return with direction and is_seed as the first columns for readability
              df[, c("direction", "isoform", "is_seed",
                     "degree", "clustering", "distance_from_seed",
                     "mean_pi", "pi_entropy")]
            }

            up_df   <- .prep(gr@up_stats,   "up")
            down_df <- .prep(gr@down_stats, "down")

            result <- do.call(rbind, list(up_df, down_df))

            if (is.null(result) || nrow(result) == 0L)
              stop(sprintf(
                "No isoform statistics available for gene '%s' in analysis '%s'.",
                gene_id, analysis
              ), call. = FALSE)

            rownames(result) <- NULL
            result
          })
