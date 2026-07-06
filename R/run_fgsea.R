## run_fgsea.R
## Gene Set Enrichment Analysis (GSEA) on a DifferentialPulldownAnalysis result.
##
## Genes are ranked by a per-gene network topology delta (up minus down)
## derived from the cross_gene_summary table, then tested against GO gene sets
## (or a user-supplied pathway list) with fgsea::fgsea().
##
## Two ranking tiers are supported:
##   Tier 1 - built-in presets via `rank_by`  (delta of any paired up/down
##             column in cross_gene_summary)
##   Tier 3 - full override via `stats`  (user supplies a named numeric vector)
##
## Results are stored inside the parent DifferentialPulldownAnalysis object
## under @gsea_results[[label]], which is written back into the IsoNetAnalyzeR
## session object.
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: FgseaResult class, run_fgsea, get_gsea_result, show(FgseaResult)

#' @include isonetanalyzer_class.R
NULL
#' @include diff_pulldown_analysis_class.R
NULL


# ==============================================================================
#  Valid rank_by presets
# ==============================================================================

#' Named character vector of all valid rank_by presets
#'
#' Each name is the user-facing preset string; the value is a short description
#' used in messages and stored in the result object.
#' @noRd
.FGSEA_RANK_PRESETS <- c(
  delta_nodes          = "up_nodes - down_nodes",
  delta_edges          = "up_edges - down_edges",
  delta_avg_degree     = "up_avg_degree - down_avg_degree",
  delta_density        = "up_density - down_density",
  delta_clustering     = "up_clustering - down_clustering",
  delta_starness       = "up_starness - down_starness",
  delta_genes          = "up_genes - down_genes",
  delta_seed_mean_pi   = "up_seed_mean_pi - down_seed_mean_pi",
  delta_seed_pi_entropy = "up_seed_pi_entropy - down_seed_pi_entropy",
  edges_gained_vs_lost = "n_edges_gained - n_edges_lost"
)


# ==============================================================================
#  S4 class: FgseaResult
# ==============================================================================

#' FgseaResult - GSEA result on a differential pulldown analysis
#'
#' Produced by \code{\link{run_fgsea}}.  Stores the ranked gene vector, the
#' fgsea result table, and all run parameters.  One \linkS4class{FgseaResult}
#' object can be stored per user-supplied label inside a
#' \linkS4class{DifferentialPulldownAnalysis}; retrieve with
#' \code{get_gsea_result(obj, analysis, label)}.
#'
#' @slot label        Character scalar.  User-supplied label, used as the
#'   \code{gsea_label} argument in
#'   \code{get_gsea_result(obj, analysis, gsea_label)}.
#' @slot parent_label Character scalar.  Label of the parent
#'   \linkS4class{DifferentialPulldownAnalysis} run.
#' @slot result_table Data frame.  Output of \code{fgsea::fgsea()} with
#'   standard columns: \code{pathway}, \code{pval}, \code{padj}, \code{ES},
#'   \code{NES}, \code{nMoreExtreme}, \code{size}, \code{leadingEdge}.
#'   Sorted by \code{padj} ascending.  Access via
#'   \code{get_gsea_table(obj, analysis, gsea_label)}.
#' @slot ranked_genes Named numeric vector.  The gene-level rank statistic
#'   passed to fgsea, sorted descending.  Names are gene IDs matching the
#'   pathway gene sets.  Access as a data frame via
#'   \code{get_gsea_stats(obj, analysis, gsea_label)}.
#' @slot rank_by      Character scalar.  The preset used (\code{rank_by}
#'   argument), or \code{"custom"} when the user supplied \code{stats}
#'   directly.
#' @slot go_category  Character scalar.  GO subcollection(s) fetched from
#'   msigdbr, or \code{NA} when a user-supplied pathway list was used.
#' @slot params       Named list.  Run parameters: \code{minSize},
#'   \code{maxSize}, \code{nperm}, \code{n_genes_ranked},
#'   \code{n_genes_dropped_na}.
#' @slot created_at   \code{POSIXct} scalar.  Timestamp of the run.
#'
#' @exportClass FgseaResult
setClass(
  "FgseaResult",
  representation(
    label        = "character",
    parent_label = "character",
    result_table = "data.frame",
    ranked_genes = "numeric",
    rank_by      = "character",
    go_category  = "character",
    params       = "list",
    created_at   = "POSIXct"
  )
)


# ==============================================================================
#  show()
# ==============================================================================

#' @describeIn FgseaResult Print a concise summary.
#' @param object A \linkS4class{FgseaResult} object.
setMethod("show", "FgseaResult", function(object) {
  rt      <- object@result_table
  n_sig   <- if (nrow(rt) > 0L) sum(rt$padj < 0.05, na.rm = TRUE) else 0L
  top_row <- if (nrow(rt) > 0L) rt[1L, ] else NULL

  cli::cli_rule("FgseaResult")
  cli::cli_text("Label        : {object@label}")
  cli::cli_text("Parent run   : {object@parent_label}")
  cli::cli_text("Created      : {format(object@created_at, '%Y-%m-%d %H:%M:%S')}")
  cli::cli_text("Rank by      : {object@rank_by}")
  cli::cli_text("GO category  : {object@go_category}")
  cli::cli_text("Genes ranked : {object@params$n_genes_ranked}  ({object@params$n_genes_dropped_na} dropped - NA rank)")
  cli::cli_text("Pathways     : {nrow(rt)} tested  |  {n_sig} significant (padj < 0.05)")

  if (!is.null(top_row))
    cli::cli_text(
      "Top pathway  : {top_row$pathway}  (NES = {round(top_row$NES, 3)}, padj = {signif(top_row$padj, 3)})"
    )

  cli::cli_text(
    "Params       : minSize = {object@params$minSize}  maxSize = {object@params$maxSize}  nperm = {object@params$nperm}"
  )
  cli::cli_rule()
  invisible(object)
})


# ==============================================================================
#  Exported: run_fgsea
# ==============================================================================

#' Run GSEA on a differential pulldown analysis
#'
#' Ranks genes by a network-topology delta derived from the
#' \code{cross_gene_summary} table of an existing
#' \linkS4class{DifferentialPulldownAnalysis} run, then tests the ranked list
#' against GO gene sets (or a user-supplied pathway list) using
#' \code{fgsea::fgsea()}.
#'
#' @section Ranking tiers:
#' \describe{
#'   \item{Tier 1 - built-in preset (\code{rank_by})}{
#'     A named delta computed as the up-condition metric minus the
#'     down-condition metric from \code{cross_gene_summary}.  Positive values
#'     indicate a gene whose up-regulated isoform produced a larger network
#'     effect; negative values indicate the opposite.
#'
#'     Available presets:
#'     \itemize{
#'       \item \code{"delta_nodes"} - up_nodes \eqn{-} down_nodes \strong{(default)}
#'       \item \code{"delta_edges"} - up_edges \eqn{-} down_edges
#'       \item \code{"delta_avg_degree"} - up_avg_degree \eqn{-} down_avg_degree
#'       \item \code{"delta_density"} - up_density \eqn{-} down_density
#'       \item \code{"delta_clustering"} - up_clustering \eqn{-} down_clustering
#'       \item \code{"delta_starness"} - up_starness \eqn{-} down_starness
#'       \item \code{"delta_genes"} - up_genes \eqn{-} down_genes
#'       \item \code{"delta_seed_mean_pi"} - up_seed_mean_pi \eqn{-} down_seed_mean_pi
#'       \item \code{"delta_seed_pi_entropy"} - up_seed_pi_entropy \eqn{-} down_seed_pi_entropy
#'       \item \code{"edges_gained_vs_lost"} - n_edges_gained \eqn{-} n_edges_lost
#'     }
#'   }
#'   \item{Tier 3 - full override (\code{stats})}{
#'     Pass a named numeric vector directly (names must be gene IDs matching
#'     the pathway gene sets).  This completely bypasses the built-in ranking
#'     logic and \code{rank_by} is ignored.  Useful for complex custom rankings:
#'
#'     \preformatted{
#'     cgs      <- get_cross_gene_summary(obj, "my_run")
#'     my_stats <- setNames(
#'       (cgs$up_nodes - cgs$down_nodes) * 0.7 +
#'       (cgs$up_edges - cgs$down_edges) * 0.3,
#'       cgs$gene_id
#'     )
#'     obj <- run_fgsea(obj, "my_run", label = "weighted", stats = my_stats)
#'     }
#'   }
#' }
#'
#' @section GO pathways:
#' When \code{pathways = NULL} (default), gene sets are fetched from MSigDB
#' via \code{msigdbr} using category \code{"C5"} (GO gene sets).  The
#' \code{go_category} argument selects the subcollection(s):
#' \code{"GO:BP"} (biological process), \code{"GO:MF"} (molecular function),
#' \code{"GO:CC"} (cellular component), or \code{"ALL"} for all three.
#' Gene sets are keyed by Ensembl gene ID (\code{ensembl_gene}).
#'
#' @param obj        An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis   Character scalar.  Label of the
#'   \linkS4class{DifferentialPulldownAnalysis}; use \code{list_analyses(obj)}
#'   to see available labels.
#' @param label      Character scalar.  Key under which the result will be
#'   stored; retrieve it afterwards with
#'   \code{get_gsea_result(obj, analysis, label)}.
#' @param rank_by    Character scalar.  One of the built-in ranking presets
#'   (see Ranking tiers above).  Ignored when \code{stats} is supplied.
#'   Default \code{"delta_nodes"}.
#' @param stats      Named numeric vector or \code{NULL}.  When supplied,
#'   this vector is used directly as the fgsea rank statistic, bypassing
#'   \code{rank_by}.  Names must be gene IDs.  Default \code{NULL}.
#' @param pathways   Named list of character vectors (gene IDs per pathway),
#'   or \code{NULL} (default) to fetch GO gene sets automatically via
#'   \code{msigdbr}.
#' @param go_category Character scalar.  GO subcollection to fetch when
#'   \code{pathways = NULL}.  One of \code{"GO:BP"}, \code{"GO:MF"},
#'   \code{"GO:CC"}, or \code{"ALL"}.  Default \code{"GO:BP"}.
#' @param minSize    Integer.  Minimum gene set size after filtering to genes
#'   present in the ranked list.  Default \code{15}.
#' @param maxSize    Integer.  Maximum gene set size.  Default \code{500}.
#' @param nperm      Integer.  Number of permutations for p-value estimation.
#'   Default \code{1000}.
#' @param ...        Additional arguments passed to \code{fgsea::fgsea()}.
#'
#' @return The \linkS4class{IsoNetAnalyzeR} session object, invisibly, with the
#'   \linkS4class{FgseaResult} stored under \code{label} inside the analysis.
#'   Access via \code{get_gsea_result(obj, analysis, label)}.
#'
#' @examples
#' \dontrun{
#' # Tier 1 - built-in preset
#' obj <- run_fgsea(obj, "my_run", label = "nodes_gsea")
#' obj <- run_fgsea(obj, "my_run", label = "edges_gsea", rank_by = "delta_edges")
#'
#' # Tier 3 - custom ranking
#' cgs      <- get_cross_gene_summary(obj, "my_run")
#' my_stats <- setNames(
#'   (cgs$up_nodes - cgs$down_nodes) * 0.7 +
#'   (cgs$up_edges - cgs$down_edges) * 0.3,
#'   cgs$gene_id
#' )
#' obj <- run_fgsea(obj, "my_run", label = "weighted", stats = my_stats)
#'
#' # Retrieve filtered results (use get_gsea_table / get_gsea_stats)
#' get_gsea_table(obj, "my_run", "nodes_gsea", padj_cutoff = 0.05)
#' get_gsea_stats(obj, "my_run", "nodes_gsea")
#' }
#'
#' @export
run_fgsea <- function(obj,
                      analysis,
                      label,
                      rank_by      = "delta_nodes",
                      stats        = NULL,
                      pathways     = NULL,
                      go_category  = c("GO:BP", "GO:MF", "GO:CC", "ALL"),
                      minSize      = 15L,
                      maxSize      = 500L,
                      nperm        = 1000L,
                      ...) {

  # -- Package checks ----------------------------------------------------------
  if (!requireNamespace("fgsea",    quietly = TRUE))
    stop("Package 'fgsea' is required. Install with:\n",
         "  BiocManager::install('fgsea')", call. = FALSE)

  # -- Input validation --------------------------------------------------------
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  diff_res <- .resolve_diff_analysis(obj, analysis)

  if (!is.character(label) || length(label) != 1L || nchar(trimws(label)) == 0L)
    stop("`label` must be a non-empty character scalar.", call. = FALSE)

  go_category <- match.arg(go_category)
  minSize     <- as.integer(minSize)
  maxSize     <- as.integer(maxSize)
  nperm       <- as.integer(nperm)

  # -- Build ranked stats vector -----------------------------------------------
  if (!is.null(stats)) {
    # Tier 3: user-supplied override
    if (!is.numeric(stats) || is.null(names(stats)))
      stop("`stats` must be a named numeric vector.", call. = FALSE)
    if (any(is.null(names(stats))) || any(nchar(names(stats)) == 0L))
      stop("All elements of `stats` must have non-empty names (gene IDs).",
           call. = FALSE)

    rank_label <- "custom"
    stats_vec  <- stats

    message(sprintf(
      "Using user-supplied stats vector (%d genes).", length(stats_vec)
    ))

  } else {
    # Tier 1: built-in preset
    rank_by <- match.arg(rank_by, names(.FGSEA_RANK_PRESETS))
    rank_label <- rank_by

    cgs        <- diff_res@cross_gene_summary
    stats_vec  <- .compute_rank_preset(cgs, rank_by)

    message(sprintf(
      "Ranking genes by preset '%s'  [%s].",
      rank_by, .FGSEA_RANK_PRESETS[rank_by]
    ))
  }

  # Drop NAs and report
  n_total   <- length(stats_vec)
  na_mask   <- is.na(stats_vec)
  n_dropped <- sum(na_mask)

  if (n_dropped > 0L)
    message(sprintf(
      "  Dropping %d gene(s) with NA rank statistic (%d remain).",
      n_dropped, n_total - n_dropped
    ))

  stats_vec <- stats_vec[!na_mask]

  if (length(stats_vec) == 0L)
    stop("No genes with a non-NA rank statistic - cannot run GSEA.", call. = FALSE)

  # fgsea requires the stats vector sorted descending
  stats_vec <- sort(stats_vec, decreasing = TRUE)

  # -- Build / fetch pathways --------------------------------------------------
  if (is.null(pathways)) {
    pathways   <- .fetch_go_pathways(go_category)
    go_cat_str <- go_category
  } else {
    if (!is.list(pathways) || is.null(names(pathways)))
      stop("`pathways` must be a named list of character vectors.", call. = FALSE)
    go_cat_str <- NA_character_
    message(sprintf("Using user-supplied pathway list (%d gene sets).",
                    length(pathways)))
  }

  # -- Run fgsea ---------------------------------------------------------------
  message(sprintf(
    "Running fgsea: %d genes ranked, %d pathway(s), minSize = %d, maxSize = %d, nperm = %d ...",
    length(stats_vec), length(pathways), minSize, maxSize, nperm
  ))

  result_tbl <- fgsea::fgsea(
    pathways = pathways,
    stats    = stats_vec,
    minSize  = minSize,
    maxSize  = maxSize,
    nperm    = nperm,
    ...
  )

  # Sort by adjusted p-value ascending
  result_tbl <- result_tbl[order(result_tbl$padj, result_tbl$pval), ]

  n_sig <- sum(result_tbl$padj < 0.05, na.rm = TRUE)
  message(sprintf(
    "Done. %d pathway(s) tested; %d significant at padj < 0.05.",
    nrow(result_tbl), n_sig
  ))

  # -- Assemble result object --------------------------------------------------
  gsea_result <- new(
    "FgseaResult",
    label        = label,
    parent_label = analysis,
    result_table = as.data.frame(result_tbl),
    ranked_genes = stats_vec,
    rank_by      = rank_label,
    go_category  = go_cat_str,
    params       = list(
      minSize            = minSize,
      maxSize            = maxSize,
      nperm              = nperm,
      n_genes_ranked     = length(stats_vec),
      n_genes_dropped_na = n_dropped
    ),
    created_at   = Sys.time()
  )

  # -- Store and return --------------------------------------------------------
  diff_res@gsea_results[[label]] <- gsea_result
  obj@analyses[[analysis]]       <- diff_res

  message(sprintf(paste0(
    "Done. Access results with:\n",
    "  get_gsea_result(obj, \"%s\", \"%s\")"),
    analysis, label
  ))

  invisible(obj)
}


# ==============================================================================
#  Internal: compute built-in rank presets
# ==============================================================================

#' Compute a named rank statistic vector from a cross_gene_summary data frame
#'
#' @param cgs     Data frame. cross_gene_summary from a
#'   DifferentialPulldownAnalysis.
#' @param rank_by Character scalar. One of the keys in .FGSEA_RANK_PRESETS.
#' @return Named numeric vector (names = gene_id).
#' @noRd
.compute_rank_preset <- function(cgs, rank_by) {

  # Column pairs for each preset: list(up_col, down_col) - or a special case
  preset_cols <- list(
    delta_nodes           = c("up_nodes",          "down_nodes"),
    delta_edges           = c("up_edges",          "down_edges"),
    delta_avg_degree      = c("up_avg_degree",     "down_avg_degree"),
    delta_density         = c("up_density",        "down_density"),
    delta_clustering      = c("up_clustering",     "down_clustering"),
    delta_starness        = c("up_starness",        "down_starness"),
    delta_genes           = c("up_genes",           "down_genes"),
    delta_seed_mean_pi    = c("up_seed_mean_pi",    "down_seed_mean_pi"),
    delta_seed_pi_entropy = c("up_seed_pi_entropy", "down_seed_pi_entropy"),
    edges_gained_vs_lost  = c("n_edges_gained",     "n_edges_lost")
  )

  cols <- preset_cols[[rank_by]]
  up_col   <- cols[1L]
  down_col <- cols[2L]

  # Guard: check columns are present
  missing <- setdiff(c(up_col, down_col), names(cgs))
  if (length(missing))
    stop(sprintf(
      "cross_gene_summary is missing column(s) required for preset '%s': %s",
      rank_by, paste(missing, collapse = ", ")
    ), call. = FALSE)

  stat <- cgs[[up_col]] - cgs[[down_col]]
  setNames(stat, cgs$gene_id)
}


# ==============================================================================
#  Internal: fetch GO pathways from msigdbr
# ==============================================================================

#' Fetch GO gene sets from MSigDB via msigdbr
#'
#' Downloads the C5 collection (GO gene sets) for Homo sapiens and returns a
#' named list of Ensembl gene ID vectors, one entry per pathway.
#'
#' @param go_category Character scalar. One of "GO:BP", "GO:MF", "GO:CC",
#'   or "ALL".
#' @return Named list of character vectors (Ensembl gene IDs per pathway).
#' @noRd
.fetch_go_pathways <- function(go_category) {

  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Package 'msigdbr' is required to fetch GO pathways automatically.\n",
         "Install with:  install.packages('msigdbr')\n",
         "Or supply your own pathway list via the `pathways` argument.",
         call. = FALSE)

  message(sprintf("Fetching GO gene sets from MSigDB (C5, subcollection: %s) ...",
                  go_category))

  if (go_category == "ALL") {
    go_sets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5")
    # Restrict to GO subcollections only (exclude HPO)
    go_sets <- go_sets[go_sets$gs_subcollection %in% c("GO:BP", "GO:MF", "GO:CC"), ]
  } else {
    go_sets <- msigdbr::msigdbr(
      species       = "Homo sapiens",
      collection    = "C5",
      subcollection = go_category
    )
  }

  # Drop rows with missing Ensembl IDs
  go_sets <- go_sets[!is.na(go_sets$ensembl_gene) &
                       nchar(go_sets$ensembl_gene) > 0L, ]

  if (nrow(go_sets) == 0L)
    stop(sprintf(
      "msigdbr returned no gene sets for subcollection '%s'. ",
      go_category
    ), call. = FALSE)

  pathways <- split(go_sets$ensembl_gene, go_sets$gs_name)

  message(sprintf("  %d gene set(s) retrieved.", length(pathways)))
  pathways
}


# ==============================================================================
#  Exported: get_gsea_result
# ==============================================================================

#' Retrieve a stored GSEA result from a differential pulldown analysis
#'
#' @param obj        An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis   Character scalar.  Label of the
#'   \linkS4class{DifferentialPulldownAnalysis}.
#' @param gsea_label Character scalar.  Label of the
#'   \linkS4class{FgseaResult} stored inside the analysis.
#'
#' @return A \linkS4class{FgseaResult} object.
#'
#' @examples
#' \dontrun{
#' # show() summary
#' get_gsea_result(obj, "my_run", "nodes_gsea")
#'
#' # For data frames, use the dedicated table accessors:
#' get_gsea_table(obj, "my_run", "nodes_gsea", padj_cutoff = 0.05)
#' get_gsea_stats(obj, "my_run", "nodes_gsea")
#' }
#'
#' @export
setGeneric("get_gsea_result",
           function(obj, analysis, gsea_label)
             standardGeneric("get_gsea_result"))

#' @describeIn get_gsea_result Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_gsea_result", "IsoNetAnalyzeR",
          function(obj, analysis, gsea_label) {

            diff_res <- .resolve_diff_analysis(obj, analysis)

            if (!gsea_label %in% names(diff_res@gsea_results))
              stop(sprintf(
                "No GSEA result labelled '%s' in analysis '%s'. Available: %s",
                gsea_label, analysis,
                paste(names(diff_res@gsea_results), collapse = ", ")
              ), call. = FALSE)

            diff_res@gsea_results[[gsea_label]]
          })
