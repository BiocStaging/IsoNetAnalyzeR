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
#' @slot result_table \code{DataFrame}.  Output of \code{fgsea::fgsea()} with
#'   standard columns: \code{pathway}, \code{pval}, \code{padj}, \code{ES},
#'   \code{NES}, \code{nMoreExtreme}, \code{size}, \code{leadingEdge}.
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
#'   \code{maxSize}, \code{nperm}, \code{n_genes_ranked},
#'   \code{n_genes_dropped_na}.
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

  .cli_aligned(
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
#'     A named delta computed as the up-seed metric minus the down-seed metric
#'     from \code{cross_gene_summary}. Positive values indicate a larger
#'     static predicted neighbourhood metric for the up-labelled seeds;
#'     negative values indicate the opposite.
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
#' (\code{\link{isonet_gene_names}}).
#'
#' When \code{pathways = NULL} (default), the GO gene sets of MSigDB
#' release 2023.1 (collection C5) are fetched with the msigdb package from
#' ExperimentHub. \code{go_category} selects the subcollection(s):
#' \code{"GO:BP"} (biological process), \code{"GO:MF"} (molecular function),
#' \code{"GO:CC"} (cellular component), or \code{"ALL"} for all three.
#' fgsea estimates p-values by permutation, so set the random seed before
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
#' @param nperm      Integer.  Number of permutations for p-value estimation.
#'   Default \code{1000}.
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
                      nperm        = 1000L,
                      ...) {

  obj_name <- .obj_name(substitute(obj))

  # -- Package checks ----------------------------------------------------------
  if (!requireNamespace("fgsea",    quietly = TRUE))
    stop("Package 'fgsea' is required. Install with:\n",
         "  BiocManager::install('fgsea')", call. = FALSE)

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
  nperm   <- .validate_positive_integer(nperm, "nperm")
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
  network_genes <- unique(igraph::V(obj@network)$gene)
  if (is.null(pathways)) {
    pathways   <- .fetch_go_pathways(go_category, network_genes)
    go_cat_str <- sprintf("%s (MSigDB %s)", go_category, .MSIGDB_VERSION)
  } else {
    pathways   <- .as_pathway_list(pathways, network_genes)
    go_cat_str <- NA_character_
    message(sprintf("Using user-supplied gene sets (%d).", length(pathways)))
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
    result_table = .as_DataFrame(result_tbl),
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
#  Internal: gene sets from MSigDB (msigdb package) or the user
# ==============================================================================

#' MSigDB release fetched when \code{pathways = NULL}
#' @noRd
.MSIGDB_VERSION <- "2023.1"

#' Fetch GO gene sets from MSigDB via the msigdb package
#'
#' Loads the human MSigDB collection (gene symbols) from ExperimentHub,
#' keeps the requested C5 GO subcollection(s) and maps the symbols to
#' Ensembl gene IDs.
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
  gsc <- msigdb::getMsigdb(org = "hs", id = "SYM", version = .MSIGDB_VERSION)

  subcollections <- if (go_category == "ALL")
    c("GO:BP", "GO:MF", "GO:CC") else go_category
  types <- lapply(gsc, GSEABase::collectionType)
  category <- vapply(types, GSEABase::bcCategory, character(1L))
  subcategory <- vapply(types, function(x) {
    sub <- GSEABase::bcSubCategory(x)
    if (length(sub)) sub else NA_character_
  }, character(1L))
  go_sets <- gsc[category == "c5" & subcategory %in% subcollections]

  pathways <- .as_pathway_list(go_sets, genes)
  if (length(pathways) == 0L)
    stop(sprintf("No GO gene sets were found for subcollection '%s'.",
                 go_category), call. = FALSE)

  message(sprintf("  %d gene set(s) retrieved.", length(pathways)))
  pathways
}

#' Gene sets as a named list of Ensembl gene IDs
#'
#' Accepts a named list of character vectors or a GSEABase
#' \code{GeneSetCollection} whose gene IDs are Ensembl gene IDs or gene
#' symbols. Symbols are mapped to Ensembl gene IDs with the Ensembl 107
#' annotation (\code{isonet_gene_names()}), restricted to \code{genes} so
#' that a symbol shared by alternative-haplotype genes maps to the network
#' gene. Gene sets left empty are dropped.
#'
#' @param pathways Named list or \code{GeneSetCollection}.
#' @param genes Character vector of Ensembl gene IDs present in the network.
#' @return Named list of character vectors.
#' @noRd
.as_pathway_list <- function(pathways, genes) {
  if (is(pathways, "GeneSetCollection")) {
    if (!requireNamespace("GSEABase", quietly = TRUE))
      stop("Package 'GSEABase' is required for GeneSetCollection input.",
           call. = FALSE)
    id_types <- unique(vapply(
      lapply(pathways, GSEABase::geneIdType), function(x) class(x)[1L],
      character(1L)
    ))
    ids <- GSEABase::geneIds(pathways)
    if (identical(id_types, "SymbolIdentifier")) {
      symbols <- isonet_gene_names()
      symbols <- symbols[!is.na(symbols$gene_name) &
                           symbols$gene_id %in% genes, ]
      by_symbol <- split(symbols$gene_id, symbols$gene_name)
      ids <- lapply(ids, function(x) unique(unlist(by_symbol[x],
                                                   use.names = FALSE)))
    } else if (!identical(id_types, "ENSEMBLIdentifier")) {
      stop("GeneSetCollection gene IDs must be gene symbols or Ensembl gene ",
           "IDs; found ", paste(id_types, collapse = ", "), ". Convert them ",
           "with GSEABase::mapIdentifiers().", call. = FALSE)
    }
    ids <- ids[lengths(ids) > 0L]
    return(ids)
  }
  if (!is.list(pathways) || is.null(names(pathways)))
    stop("`pathways` must be a named list of character vectors or a ",
         "GeneSetCollection.", call. = FALSE)
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
