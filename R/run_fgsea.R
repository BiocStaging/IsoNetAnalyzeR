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

#' The rank_by presets and the cross_gene_summary columns they subtract
#'
#' Each name is the user-facing preset string; the value is the pair of
#' columns whose difference (first minus second) ranks the genes.
#' @noRd
.FGSEA_RANK_PRESETS <- list(
  delta_nodes           = c("up_nodes",           "down_nodes"),
  delta_edges           = c("up_edges",           "down_edges"),
  delta_avg_degree      = c("up_avg_degree",      "down_avg_degree"),
  delta_density         = c("up_density",         "down_density"),
  delta_clustering      = c("up_clustering",      "down_clustering"),
  delta_starness        = c("up_starness",        "down_starness"),
  delta_genes           = c("up_genes",           "down_genes"),
  delta_seed_mean_pi    = c("up_seed_mean_pi",    "down_seed_mean_pi"),
  delta_seed_pi_entropy = c("up_seed_pi_entropy", "down_seed_pi_entropy"),
  edges_gained_vs_lost  = c("n_edges_gained",     "n_edges_lost")
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
#' @slot result_table \code{DataFrame}.  Output of \code{fgsea::fgsea()} with
#'   standard columns: \code{pathway}, \code{pval}, \code{padj}, \code{ES},
#'   \code{NES}, \code{size}, \code{leadingEdge}, and \code{log2err}
#'   (\code{fgseaMultilevel}) or \code{nMoreExtreme} (\code{fgseaSimple}).
#'   Sorted by \code{padj} ascending.  Access via
#'   \code{get_gsea_table(obj, analysis, gsea_label)}.
#' @slot ranked_genes Named numeric vector.  The gene-level rank statistic
#'   passed to fgsea, sorted descending.  Names are gene IDs matching the
#'   pathway gene sets.  Access as a table via
#'   \code{get_gsea_stats(obj, analysis, gsea_label)}.
#' @slot rank_by      Character scalar.  The preset used (\code{rank_by}
#'   argument), or \code{"custom"} when the user supplied \code{stats}
#'   directly.
#' @slot go_category  Character scalar.  GO subcollection(s) and MSigDB
#'   release fetched with the msigdb package, or \code{NA} when gene sets
#'   were supplied.
#' @slot params       Named list.  Run parameters: \code{minSize},
#'   \code{maxSize}, \code{nperm} (\code{NULL} unless set),
#'   \code{method} (\code{"fgseaMultilevel"} or \code{"fgseaSimple"}),
#'   \code{n_genes_ranked}, \code{n_genes_dropped_na}.
#' @slot created_at   \code{POSIXct} scalar.  Timestamp of the run.
#'
#' @return Objects of this class are created by \code{\link{run_fgsea}} and returned by \code{\link{get_gsea_result}}.
#'
#' @examples
#' seeds <- read.csv(system.file("extdata", "gse141432_t2d_vs_lean_seeds.csv",
#'                               package = "IsoNetAnalyzeR"))
#' # Six genes for speed; the vignette analyses all 65. With few ranked genes,
#' # minSize must be small.
#' genes <- c("CLTCL1", "WAS", "VAV1", "HLA-DMA", "ICOS", "PRKACB")
#' obj <- IsoNetAnalyzeR("example")
#' obj <- run_differential_pulldown(
#'   obj, "t2d", seeds[seeds$gene_name %in% genes, ],
#'   order = 1, confirm = FALSE, verbose = FALSE
#' )
#' # Two small gene sets of gene symbols, for illustration. MSigDB collections
#' # are available as GeneSetCollections from the msigdb package; see the
#' # vignette.
#' gene_sets <- GSEABase::GeneSetCollection(list(
#'   GSEABase::GeneSet(c("WAS", "VAV1", "ICOS"), setName = "set_A",
#'                     geneIdType = GSEABase::SymbolIdentifier()),
#'   GSEABase::GeneSet(c("CLTCL1", "PRKACB", "HLA-DMA"), setName = "set_B",
#'                     geneIdType = GSEABase::SymbolIdentifier())
#' ))
#' set.seed(42)
#' obj <- run_fgsea(obj, "t2d", label = "example_gsea", pathways = gene_sets,
#'                  minSize = 3)
#' get_gsea_result(obj, "t2d", "example_gsea")
#'
#' @exportClass FgseaResult
setClass(
  "FgseaResult",
  representation(
    label        = "character",
    parent_label = "character",
    result_table = "DataFrame",
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
  rt      <- .df(object@result_table)
  n_sig   <- if (nrow(rt) > 0L) sum(rt$padj < 0.05, na.rm = TRUE) else 0L
  top_row <- if (nrow(rt) > 0L) rt[1L, ] else NULL

  cli::cli_rule("FgseaResult")
  .cli_aligned("Label        : {object@label}")
  .cli_aligned("Parent run   : {object@parent_label}")
  .cli_aligned("Created      : {format(object@created_at, '%Y-%m-%d %H:%M:%S')}")
  .cli_aligned("Rank by      : {object@rank_by}")
  gene_sets <- if (is.na(object@go_category)) "user-supplied" else
    object@go_category
  .cli_aligned("Gene sets    : {gene_sets}")
  .cli_aligned("Genes ranked : {object@params$n_genes_ranked}  ({object@params$n_genes_dropped_na} dropped - NA rank)")
  .cli_aligned("Pathways     : {nrow(rt)} tested  |  {n_sig} significant (padj < 0.05)")

  if (!is.null(top_row))
    .cli_aligned(
      "Top pathway  : {top_row$pathway}  (NES = {round(top_row$NES, 3)}, padj = {signif(top_row$padj, 3)})"
    )

  method <- object@params$method %||% "fgseaSimple"
  if (!is.null(object@params$nperm))
    method <- sprintf("%s, nperm = %d", method, object@params$nperm)
  .cli_aligned(
    "Params       : minSize = {object@params$minSize}  maxSize = {object@params$maxSize}  method = {method}"
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
#'     A named delta computed as the up-seed metric minus the down-seed metric
#'     from \code{cross_gene_summary}. Positive values indicate a larger
#'     static predicted neighbourhood metric for the up-labelled seeds;
#'     negative values indicate the opposite. Node, edge and gene counts
#'     include the seeds themselves, so a direction with more seeds counts
#'     more; compare \code{n_up_resolved} and \code{n_down_resolved}.
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
#' @section Gene sets:
#' \code{pathways} may be a named list of Ensembl gene ID vectors or a
#' GSEABase \code{GeneSetCollection}, for example a collection from the
#' msigdb package. Gene symbols in a \code{GeneSetCollection} are mapped to
#' Ensembl gene IDs with the Ensembl 107 annotation
#' (\code{\link{isonet_gene_names}}). The gene IDs of a collection without
#' an identifier type, such as one read with \code{GSEABase::getGmt()}, are
#' taken as Ensembl gene IDs when they all look like Ensembl gene IDs, and as
#' gene symbols otherwise. Gene sets are matched to the IDs of the ranked
#' genes (the analysed gene IDs, or the names of \code{stats}): IDs that
#' are ranked genes are kept, and other Ensembl gene IDs match the ranked
#' genes with or without a version suffix. Gene sets with no ranked gene are
#' an error.
#'
#' When \code{pathways = NULL} (default), the GO gene sets of MSigDB
#' release 2023.1 (collection C5) are fetched with the msigdb package from
#' ExperimentHub and kept for the R session. \code{go_category} selects the
#' subcollection(s): \code{"GO:BP"} (biological process), \code{"GO:MF"}
#' (molecular function), \code{"GO:CC"} (cellular component), or
#' \code{"ALL"} for all three.
#' fgsea estimates p-values by random sampling, so set the random seed before
#' this call.
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
#' @param pathways   Named list of character vectors (Ensembl gene IDs per
#'   pathway), a \code{GeneSetCollection}, or \code{NULL} (default) to fetch
#'   MSigDB GO gene sets with the msigdb package (see Gene sets).
#' @param go_category Character scalar.  GO subcollection to fetch when
#'   \code{pathways = NULL}.  One of \code{"GO:BP"}, \code{"GO:MF"},
#'   \code{"GO:CC"}, or \code{"ALL"}.  Default \code{"GO:BP"}.
#' @param minSize    Integer.  Minimum gene set size after filtering to genes
#'   present in the ranked list.  Default \code{15}.
#' @param maxSize    Integer.  Maximum gene set size.  Default \code{500}.
#' @param nperm      Integer or \code{NULL} (default).  With \code{NULL},
#'   fgsea estimates p-values with \code{fgseaMultilevel}, which also
#'   estimates very small p-values (down to its \code{eps}, which can be
#'   passed through \code{...}). A number runs \code{fgseaSimple} with that
#'   many permutations, which cannot give a p-value below
#'   \code{1 / (nperm + 1)}.
#' @param ...        Additional arguments passed to \code{fgsea::fgsea()}.
#'
#' @return The \linkS4class{IsoNetAnalyzeR} session object, invisibly, with the
#'   \linkS4class{FgseaResult} stored under \code{label} inside the analysis.
#'   Access via \code{get_gsea_result(obj, analysis, label)}.
#'
#' @examples
#' seeds <- read.csv(system.file("extdata", "gse141432_t2d_vs_lean_seeds.csv",
#'                               package = "IsoNetAnalyzeR"))
#' # Six genes for speed; the vignette analyses all 65. With few ranked genes,
#' # minSize must be small.
#' genes <- c("CLTCL1", "WAS", "VAV1", "HLA-DMA", "ICOS", "PRKACB")
#' obj <- IsoNetAnalyzeR("example")
#' obj <- run_differential_pulldown(
#'   obj, "t2d", seeds[seeds$gene_name %in% genes, ],
#'   order = 1, confirm = FALSE, verbose = FALSE
#' )
#' # Two small gene sets of gene symbols, for illustration. MSigDB collections
#' # are available as GeneSetCollections from the msigdb package; see the
#' # vignette.
#' gene_sets <- GSEABase::GeneSetCollection(list(
#'   GSEABase::GeneSet(c("WAS", "VAV1", "ICOS"), setName = "set_A",
#'                     geneIdType = GSEABase::SymbolIdentifier()),
#'   GSEABase::GeneSet(c("CLTCL1", "PRKACB", "HLA-DMA"), setName = "set_B",
#'                     geneIdType = GSEABase::SymbolIdentifier())
#' ))
#' set.seed(42)
#' obj <- run_fgsea(obj, "t2d", label = "example_gsea", pathways = gene_sets,
#'                  minSize = 3)
#' get_gsea_table(obj, "t2d", "example_gsea")
#'
#' # Custom ranking: a weighted node and edge difference per gene
#' cgs <- get_cross_gene_summary(obj, "t2d")
#' my_stats <- setNames(
#'   (cgs$up_nodes - cgs$down_nodes) * 0.7 +
#'     (cgs$up_edges - cgs$down_edges) * 0.3,
#'   cgs$gene_id
#' )
#' obj <- run_fgsea(obj, "t2d", label = "weighted", stats = my_stats,
#'                  pathways = gene_sets, minSize = 3)
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
                      nperm        = NULL,
                      ...) {

  obj_name <- .obj_name(substitute(obj))

  # -- Input validation --------------------------------------------------------
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  diff_res <- .resolve_diff_analysis(obj, analysis)

  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(trimws(label)))
    stop("`label` must be a non-empty character scalar.", call. = FALSE)

  go_category <- match.arg(go_category)
  minSize <- .validate_positive_integer(minSize, "minSize")
  maxSize <- .validate_positive_integer(maxSize, "maxSize")
  if (!is.null(nperm)) nperm <- .validate_positive_integer(nperm, "nperm")
  method <- if (is.null(nperm)) "fgseaMultilevel" else "fgseaSimple"
  if (maxSize < minSize)
    stop("`maxSize` must be greater than or equal to `minSize`.", call. = FALSE)

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

    cgs        <- .df(diff_res@cross_gene_summary)
    stats_vec  <- .compute_rank_preset(cgs, rank_by)

    message(sprintf(
      "Ranking genes by preset '%s'  [%s].",
      rank_by, paste(.FGSEA_RANK_PRESETS[[rank_by]], collapse = " - ")
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
  # Gene sets are matched to the IDs of the ranked genes
  ranked_genes <- names(stats_vec)
  if (is.null(pathways)) {
    pathways   <- .fetch_go_pathways(go_category, ranked_genes)
    go_cat_str <- sprintf("%s (MSigDB %s)", go_category, .MSIGDB_VERSION)
  } else {
    pathways   <- .as_pathway_list(pathways, ranked_genes)
    go_cat_str <- NA_character_
    message(sprintf("Using user-supplied gene sets (%d).", length(pathways)))
  }

  # -- Run fgsea ---------------------------------------------------------------
  message(sprintf(
    "Running fgsea (%s): %d genes ranked, %d pathway(s), minSize = %d, maxSize = %d ...",
    method, length(stats_vec), length(pathways), minSize, maxSize
  ))

  # fgsea runs fgseaSimple whenever an nperm argument is passed. Two direct
  # calls: do.call() would put the whole gene sets into warning calls
  result_tbl <- if (is.null(nperm)) {
    fgsea::fgsea(pathways = pathways, stats = stats_vec, minSize = minSize,
                 maxSize = maxSize, ...)
  } else {
    fgsea::fgsea(pathways = pathways, stats = stats_vec, minSize = minSize,
                 maxSize = maxSize, nperm = nperm, ...)
  }

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
    result_table = .as_DataFrame(result_tbl),
    ranked_genes = stats_vec,
    rank_by      = rank_label,
    go_category  = go_cat_str,
    params       = list(
      minSize            = minSize,
      maxSize            = maxSize,
      nperm              = nperm,
      method             = method,
      n_genes_ranked     = length(stats_vec),
      n_genes_dropped_na = n_dropped
    ),
    created_at   = Sys.time()
  )

  # -- Store and return --------------------------------------------------------
  diff_res@gsea_results[[label]] <- gsea_result
  obj@analyses[[analysis]]       <- diff_res

  message(
    "Access results with:\n",
    sprintf("  get_gsea_result(%s, \"%s\", \"%s\")", obj_name, analysis, label)
  )

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

  cols <- .FGSEA_RANK_PRESETS[[rank_by]]
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
#  Internal: gene sets from MSigDB (msigdb package) or the user
# ==============================================================================

#' MSigDB release fetched when \code{pathways = NULL}
#' @noRd
.MSIGDB_VERSION <- "2023.1"

#' MSigDB GO gene sets, loaded once per R session
#'
#' The C5 GO:BP, GO:MF and GO:CC gene sets (gene symbols) of the human MSigDB
#' collection, with the subcollection of each set. Loading the collection
#' from ExperimentHub takes most of a minute; only its GO sets are kept.
#' @noRd
.msigdb_go_sets <- function() {
  if (is.null(.annotation_cache$msigdb_go)) {
    gsc <- msigdb::getMsigdb(org = "hs", id = "SYM",
                             version = .MSIGDB_VERSION)
    types <- lapply(gsc, GSEABase::collectionType)
    category <- vapply(types, GSEABase::bcCategory, character(1L))
    subcategory <- vapply(types, function(x) {
      sub <- GSEABase::bcSubCategory(x)
      if (length(sub)) sub else NA_character_
    }, character(1L))
    go <- category == "c5" & subcategory %in% c("GO:BP", "GO:MF", "GO:CC")
    .annotation_cache$msigdb_go <- list(sets = gsc[go],
                                        subcategory = subcategory[go])
  }
  .annotation_cache$msigdb_go
}

#' Fetch GO gene sets from MSigDB via the msigdb package
#'
#' Takes the requested C5 GO subcollection(s) of the human MSigDB collection
#' (gene symbols, see \code{.msigdb_go_sets()}) and maps the symbols to
#' Ensembl gene IDs. The result is kept for the R session, per subcollection
#' and set of network genes.
#'
#' @param go_category Character scalar. One of "GO:BP", "GO:MF", "GO:CC",
#'   or "ALL".
#' @param genes Character vector of Ensembl gene IDs that symbols may map to.
#' @return Named list of character vectors (Ensembl gene IDs per pathway).
#' @noRd
.fetch_go_pathways <- function(go_category, genes) {

  for (pkg in c("msigdb", "GSEABase"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required to fetch GO gene sets ", pkg),
           "automatically. Install it with BiocManager::install('", pkg,
           "'), or supply your own gene sets via the `pathways` argument.",
           call. = FALSE)

  message(sprintf(
    "Fetching GO gene sets (C5, %s) from MSigDB %s via the msigdb package ...",
    go_category, .MSIGDB_VERSION
  ))
  # The converted sets are kept per subcollection for the network's genes
  cached <- .annotation_cache$go_pathways[[go_category]]
  if (!identical(cached$genes, genes)) {
    go <- .msigdb_go_sets()
    subcollections <- if (go_category == "ALL")
      c("GO:BP", "GO:MF", "GO:CC") else go_category
    sets <- GSEABase::geneIds(go$sets[go$subcategory %in% subcollections])
    cached <- list(genes = genes,
                   pathways = .network_gene_sets(sets, genes, symbols = TRUE))
    .annotation_cache$go_pathways[[go_category]] <- cached
  }
  pathways <- cached$pathways
  if (length(pathways) == 0L)
    stop(sprintf("No %s gene set of MSigDB %s contains a network gene.",
                 go_category, .MSIGDB_VERSION), call. = FALSE)

  message(sprintf("  %d gene set(s) retrieved.", length(pathways)))
  pathways
}

#' Gene sets as a named list of target gene IDs
#'
#' Accepts a named list of character vectors or a GSEABase
#' \code{GeneSetCollection} whose gene IDs are Ensembl gene IDs or gene
#' symbols. IDs of a collection without an identifier type
#' (\code{NullIdentifier}, as from \code{GSEABase::getGmt()}) are taken as
#' Ensembl gene IDs when every non-blank one looks like one, and as gene
#' symbols otherwise; a named list holds IDs of the target genes. See
#' \code{.network_gene_sets()} for the mapping. Gene sets in which no gene is
#' a target gene are an error.
#'
#' @param pathways Named list or \code{GeneSetCollection}.
#' @param genes Character vector of the target gene IDs: the ranked genes
#'   (\code{run_fgsea()}) or the network genes (\code{run_component_ora()}).
#' @return Named list of character vectors.
#' @noRd
.as_pathway_list <- function(pathways, genes) {
  symbols <- FALSE
  if (is(pathways, "GeneSetCollection")) {
    if (!requireNamespace("GSEABase", quietly = TRUE))
      stop("Package 'GSEABase' is required for GeneSetCollection input.",
           call. = FALSE)
    if (!length(pathways))
      stop("The GeneSetCollection contains no gene sets.", call. = FALSE)
    id_types <- unique(vapply(
      lapply(pathways, GSEABase::geneIdType), function(x) class(x)[1L],
      character(1L)
    ))
    # IDs without a declared type (as from GSEABase::getGmt()) are Ensembl
    # gene IDs when every non-blank one looks like one, and symbols otherwise
    untyped <- identical(id_types, "NullIdentifier")
    all_ids <- unlist(GSEABase::geneIds(pathways), use.names = FALSE)
    ensembl <- identical(id_types, "ENSEMBLIdentifier") ||
      (untyped &&
         all(grepl("^ENS[A-Z]*G[0-9]+", all_ids[nzchar(all_ids)])))
    if (!ensembl && !untyped && !identical(id_types, "SymbolIdentifier"))
      stop("GeneSetCollection gene IDs must be gene symbols or Ensembl gene ",
           "IDs; found ", paste(id_types, collapse = ", "), ". Convert them ",
           "with GSEABase::mapIdentifiers().", call. = FALSE)
    symbols <- !ensembl
  } else if (!is.list(pathways) || is.null(names(pathways))) {
    stop("`pathways` must be a named list of character vectors or a ",
         "GeneSetCollection.", call. = FALSE)
  }
  collection <- is(pathways, "GeneSetCollection")
  if (collection) pathways <- GSEABase::geneIds(pathways)
  ids <- .network_gene_sets(pathways, genes, symbols)
  if (!any(unlist(ids, use.names = FALSE) %in% genes)) {
    if (collection)
      stop("No gene of the gene sets is a network gene; gene sets must hold ",
           "gene symbols or the network's Ensembl gene IDs.", call. = FALSE)
    stop("No gene of the named list is a gene of the analysis; a named list ",
         "must hold the same gene IDs (for the bundled network, Ensembl gene ",
         "IDs). For gene symbols, pass a GeneSetCollection.", call. = FALSE)
  }
  ids
}

#' Gene sets mapped to the target gene IDs
#'
#' Gene symbols are mapped with the Ensembl 107 annotation
#' (\code{isonet_gene_names()}) to every target gene with that symbol. An ID
#' that is a target gene is kept; another Ensembl gene ID matches every
#' target gene with the same ID without version or GENCODE suffix (as for
#' labels), returned as the target has it; other IDs are kept. Gene sets left
#' empty are dropped.
#'
#' @param ids Named list of character vectors.
#' @param genes Character vector of the target gene IDs.
#' @param symbols Logical; \code{ids} are gene symbols.
#' @return Named list of character vectors.
#' @noRd
.network_gene_sets <- function(ids, genes, symbols) {
  # Factor-valued sets (e.g. from split() on a factor) would index by level
  ids <- lapply(ids, as.character)
  stable <- .strip_gene_version(genes)
  if (symbols) {
    # Each target gene under its symbol, so a symbol maps to all its genes
    ann <- isonet_gene_names()
    ann <- ann[!is.na(ann$gene_name), ]
    hit <- stable %in% ann$gene_id
    by_symbol <- split(genes[hit],
                       ann$gene_name[match(stable[hit], ann$gene_id)])
    ids <- lapply(ids, function(x) unique(unlist(by_symbol[x],
                                                 use.names = FALSE)))
  } else {
    # Exact target IDs stay; other Ensembl gene IDs match every target gene
    # with the same ID without version or suffix; the rest stay as given
    ids <- lapply(ids, function(x) {
      exact <- x %in% genes
      cut <- !exact & grepl("^ENS[A-Z]*G[0-9]+", x)
      key <- .strip_gene_version(x[cut])
      unique(c(x[exact], genes[stable %in% key], x[cut][!key %in% stable],
               x[!exact & !cut]))
    })
  }
  ids[lengths(ids) > 0L]
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
#' seeds <- read.csv(system.file("extdata", "gse141432_t2d_vs_lean_seeds.csv",
#'                               package = "IsoNetAnalyzeR"))
#' # Six genes for speed; the vignette analyses all 65. With few ranked genes,
#' # minSize must be small.
#' genes <- c("CLTCL1", "WAS", "VAV1", "HLA-DMA", "ICOS", "PRKACB")
#' obj <- IsoNetAnalyzeR("example")
#' obj <- run_differential_pulldown(
#'   obj, "t2d", seeds[seeds$gene_name %in% genes, ],
#'   order = 1, confirm = FALSE, verbose = FALSE
#' )
#' # Two small gene sets of gene symbols, for illustration. MSigDB collections
#' # are available as GeneSetCollections from the msigdb package; see the
#' # vignette.
#' gene_sets <- GSEABase::GeneSetCollection(list(
#'   GSEABase::GeneSet(c("WAS", "VAV1", "ICOS"), setName = "set_A",
#'                     geneIdType = GSEABase::SymbolIdentifier()),
#'   GSEABase::GeneSet(c("CLTCL1", "PRKACB", "HLA-DMA"), setName = "set_B",
#'                     geneIdType = GSEABase::SymbolIdentifier())
#' ))
#' set.seed(42)
#' obj <- run_fgsea(obj, "t2d", label = "example_gsea", pathways = gene_sets,
#'                  minSize = 3)
#' get_gsea_result(obj, "t2d", "example_gsea")
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
