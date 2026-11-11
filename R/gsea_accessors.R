## gsea_accessors.R
## Convenience accessors that return ready-to-use data frames from GSEA and
## ORA results stored inside IsoNetAnalyzeR session objects.
##
## These complement the object-level accessors (get_gsea_result,
## get_component_ora_result) which return the full S4 result objects.  Here
## every function returns a plain data frame so results can be immediately
## inspected, filtered, or passed to ggplot2 without unpacking slots.
##
## Exported:
##   get_gsea_table()   — result table from run_fgsea(), with optional filters
##   get_gsea_stats()   — ranked gene vector from run_fgsea() as a data frame
##   get_ora_table()    — result table from run_component_ora(), with optional
##                        component_id, padj_cutoff, and top_n filters

#' @include isonetanalyzer_class.R
NULL
#' @include diff_pulldown_analysis_class.R
NULL
#' @include component_analysis.R
NULL
#' @include run_fgsea.R
NULL
#' @include run_component_ora.R
NULL


# ==============================================================================
#  get_gsea_table
# ==============================================================================

#' Get the GSEA result table from a differential pulldown analysis
#'
#' Returns the \code{fgsea::fgsea()} result table stored in a
#' \linkS4class{FgseaResult} object as a \code{DataFrame}, with an optional
#' significance filter.
#'
#' @param obj        An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis   Character scalar.  Label of the
#'   \linkS4class{DifferentialPulldownAnalysis}.
#' @param gsea_label Character scalar.  Label of the
#'   \linkS4class{FgseaResult} stored inside the analysis.
#' @param padj_cutoff Numeric scalar or \code{NULL}.  When supplied, only rows
#'   with \code{padj <= padj_cutoff} are returned.  Default \code{NULL}
#'   (return all pathways).
#'
#' @return A \code{DataFrame} with columns \code{pathway}, \code{pval},
#'   \code{padj}, \code{ES}, \code{NES}, \code{nMoreExtreme}, \code{size},
#'   \code{leadingEdge}, sorted by \code{padj} ascending.
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
#'
#' # All tested gene sets, and those passing an adjusted p-value cutoff
#' get_gsea_table(obj, "t2d", "example_gsea")
#' get_gsea_table(obj, "t2d", "example_gsea", padj_cutoff = 0.05)
#'
#' @export
setGeneric("get_gsea_table",
           function(obj, analysis, gsea_label, padj_cutoff = NULL)
             standardGeneric("get_gsea_table"))

#' @describeIn get_gsea_table Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_gsea_table", "IsoNetAnalyzeR",
          function(obj, analysis, gsea_label, padj_cutoff = NULL) {

            # -- Validate padj_cutoff ------------------------------------------
            if (!is.null(padj_cutoff)) {
              padj_cutoff <- .validate_cutoff(padj_cutoff, "padj_cutoff")
            }

            # -- Retrieve result -----------------------------------------------
            res <- get_gsea_result(obj, analysis, gsea_label)
            tbl <- .df(res@result_table)

            if (nrow(tbl) == 0L) {
              message("GSEA result table is empty.")
              return(.as_DataFrame(tbl))
            }

            # -- Apply filter --------------------------------------------------
            if (!is.null(padj_cutoff)) {
              tbl <- tbl[!is.na(tbl$padj) & tbl$padj <= padj_cutoff, ,
                         drop = FALSE]

              if (nrow(tbl) == 0L)
                message(sprintf(
                  "No pathways with padj <= %.4g in '%s'.", padj_cutoff, gsea_label
                ))
            }

            rownames(tbl) <- NULL
            .as_DataFrame(tbl)
          })


# ==============================================================================
#  get_gsea_stats
# ==============================================================================

#' Get the ranked gene statistic from a GSEA result as a table
#'
#' Returns the named numeric vector of rank statistics that was passed to
#' \code{fgsea::fgsea()} as a \code{DataFrame}, making it easy to inspect
#' which genes drove the enrichment, join with other tables, or plot a ranking
#' barcode.
#'
#' @param obj        An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis   Character scalar.  Label of the
#'   \linkS4class{DifferentialPulldownAnalysis}.
#' @param gsea_label Character scalar.  Label of the
#'   \linkS4class{FgseaResult}.
#'
#' @return A \code{DataFrame} with columns:
#'   \describe{
#'     \item{\code{gene_id}}{Gene identifier (names from the ranked vector).}
#'     \item{\code{stat}}{Rank statistic value (e.g. delta_nodes).}
#'     \item{\code{rank}}{Integer rank (1 = highest statistic).}
#'   }
#'   Sorted by \code{stat} descending (i.e. rank 1 first).
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
#'
#' # The ranked statistic used by fgsea, highest first
#' head(get_gsea_stats(obj, "t2d", "example_gsea"))
#'
#' @export
setGeneric("get_gsea_stats",
           function(obj, analysis, gsea_label)
             standardGeneric("get_gsea_stats"))

#' @describeIn get_gsea_stats Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_gsea_stats", "IsoNetAnalyzeR",
          function(obj, analysis, gsea_label) {

            res       <- get_gsea_result(obj, analysis, gsea_label)
            stats_vec <- res@ranked_genes   # already sorted descending

            .as_DataFrame(data.frame(
              gene_id = names(stats_vec),
              stat    = unname(stats_vec),
              rank    = seq_along(stats_vec),
              stringsAsFactors = FALSE
            ))
          })


# ==============================================================================
#  get_ora_table
# ==============================================================================

#' Get the ORA result table from a component analysis
#'
#' Returns the \code{fgsea::fora()} result table stored in a
#' \linkS4class{ComponentOraResult} object as a \code{DataFrame}.  Results
#' can be filtered by component, significance threshold, and/or restricted to
#' the top \emph{N} terms per component.  All filters compose — they are
#' applied in the order: \code{component_id} → \code{padj_cutoff} →
#' \code{top_n}.
#'
#' @param obj                An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis           Character scalar.  Label of the parent
#'   \linkS4class{DifferentialPulldownAnalysis}.
#' @param component_analysis Character scalar.  Label of the
#'   \linkS4class{ComponentAnalysis}.
#' @param ora_label          Character scalar.  Label of the
#'   \linkS4class{ComponentOraResult}.
#' @param component_id       Integer scalar or \code{NULL}.  When supplied,
#'   only rows for that component are returned.  Default \code{NULL} (all
#'   components).
#' @param padj_cutoff        Numeric scalar or \code{NULL}.  When supplied,
#'   only rows with \code{padj <= padj_cutoff} are returned.  Default
#'   \code{NULL} (no filter).
#' @param top_n              Integer scalar or \code{NULL}.  When supplied,
#'   only the top \emph{N} terms per component (sorted by \code{padj}
#'   ascending, ties broken by \code{pval}) are returned.  Applied after
#'   \code{padj_cutoff}.  Default \code{NULL} (all terms).
#'
#' @return A \code{DataFrame} with columns \code{component_id}, \code{pathway},
#'   \code{pval}, \code{padj}, \code{overlap}, \code{pathway_size},
#'   \code{n_genes_tested}, \code{overlap_genes}.  Sorted by
#'   \code{component_id} then \code{padj} ascending.
#'
#' @examples
#' seeds <- read.csv(system.file("extdata", "gse141432_t2d_vs_lean_seeds.csv",
#'                               package = "IsoNetAnalyzeR"))
#' genes <- c("CLTCL1", "WAS", "VAV1", "HLA-DMA", "ICOS", "PRKACB")
#' obj <- IsoNetAnalyzeR("example")
#' obj <- run_differential_pulldown(
#'   obj, "t2d", seeds[seeds$gene_name %in% genes, ],
#'   order = 1, confirm = FALSE, verbose = FALSE
#' )
#' obj <- run_component_analysis(obj, "t2d", label = "comp")
#' # Two small gene sets of gene symbols, for illustration. MSigDB collections
#' # are available as GeneSetCollections from the msigdb package; see the
#' # vignette.
#' gene_sets <- GSEABase::GeneSetCollection(list(
#'   GSEABase::GeneSet(c("WAS", "VAV1", "ICOS"), setName = "set_A",
#'                     geneIdType = GSEABase::SymbolIdentifier()),
#'   GSEABase::GeneSet(c("CLTCL1", "PRKACB", "HLA-DMA"), setName = "set_B",
#'                     geneIdType = GSEABase::SymbolIdentifier())
#' ))
#' obj <- run_component_ora(obj, "t2d", "comp", label = "example",
#'                          pathways = gene_sets, minSize = 2)
#'
#' # Top two gene sets per component
#' get_ora_table(obj, "t2d", "comp", "example", top_n = 2)
#'
#' # Significant gene sets of component 1
#' get_ora_table(obj, "t2d", "comp", "example", component_id = 1,
#'               padj_cutoff = 0.05)
#'
#' @export
setGeneric("get_ora_table",
           function(obj, analysis, component_analysis, ora_label,
                    component_id  = NULL,
                    padj_cutoff   = NULL,
                    top_n         = NULL)
             standardGeneric("get_ora_table"))

#' @describeIn get_ora_table Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_ora_table", "IsoNetAnalyzeR",
          function(obj, analysis, component_analysis, ora_label,
                   component_id  = NULL,
                   padj_cutoff   = NULL,
                   top_n         = NULL) {

            # -- Validate arguments --------------------------------------------
            if (!is.null(component_id)) {
              component_id <- .validate_positive_integer(
                component_id, "component_id"
              )
            }

            if (!is.null(padj_cutoff)) {
              padj_cutoff <- .validate_cutoff(padj_cutoff, "padj_cutoff")
            }

            if (!is.null(top_n)) {
              top_n <- .validate_positive_integer(top_n, "top_n")
            }

            # -- Retrieve result table -----------------------------------------
            ora <- get_component_ora_result(obj, analysis, component_analysis,
                                            ora_label)
            tbl <- .df(ora@result_table)

            if (nrow(tbl) == 0L) {
              message("ORA result table is empty.")
              return(.as_DataFrame(tbl))
            }

            # -- Filter 1: component_id ----------------------------------------
            if (!is.null(component_id)) {
              available_ids <- sort(unique(tbl$component_id))

              if (!component_id %in% available_ids)
                stop(sprintf(
                  "component_id %d not found in ORA results. Available: %s",
                  component_id,
                  paste(available_ids, collapse = ", ")
                ), call. = FALSE)

              tbl <- tbl[tbl$component_id == component_id, , drop = FALSE]
            }

            # -- Filter 2: padj_cutoff -----------------------------------------
            if (!is.null(padj_cutoff)) {
              tbl <- tbl[!is.na(tbl$padj) & tbl$padj <= padj_cutoff, ,
                         drop = FALSE]

              if (nrow(tbl) == 0L) {
                message(sprintf(
                  "No terms with padj <= %.4g after filtering.", padj_cutoff
                ))
                return(.as_DataFrame(tbl))
              }
            }

            # -- Filter 3: top_n per component ---------------------------------
            if (!is.null(top_n)) {
              # Table is already sorted by component_id then padj, so we just
              # take the first top_n rows within each component_id group
              tbl <- do.call(rbind, lapply(
                split(tbl, tbl$component_id),
                function(grp) head(grp, top_n)
              ))
              rownames(tbl) <- NULL
            }

            rownames(tbl) <- NULL
            .as_DataFrame(tbl)
          })
