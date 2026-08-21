## run_component_ora.R
## Over-representation analysis (ORA) for connected components from a
## ComponentAnalysis result.
##
## Each component is tested independently against GO gene sets (or a
## user-supplied pathway list) using fgsea::fora().  The background universe
## is all genes present in the full IsoNetAnalyzeR network (most principled
## choice for protein-interaction network ORA).
##
## Node-to-gene mapping uses the reliable `gene` vertex attribute on the
## igraph network object.  Nodes are isoform IDs; the `gene` attribute holds
## the Ensembl gene ID.
##
## Results are stored inside the parent ComponentAnalysis object under
## @ora_results[[label]] as a single long-format data frame covering all
## components.  Users should access and filter results via get_ora_table(),
## which supports component_id, padj_cutoff, and top_n filters.
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: ComponentOraResult class, run_component_ora,
##           get_component_ora_result, show(ComponentOraResult)

#' @importFrom stats p.adjust setNames
#' @importFrom utils head
NULL

#' @include isonetanalyzer_class.R
NULL
#' @include diff_pulldown_analysis_class.R
NULL
#' @include component_analysis.R
NULL


# ==============================================================================
#  S4 class: ComponentOraResult
# ==============================================================================

#' ComponentOraResult - GO over-representation analysis on connected components
#'
#' Produced by \code{\link{run_component_ora}}.  Stores the ORA result table
#' across all components and all run parameters.  Stored inside
#' \code{ComponentAnalysis@@ora_results[[label]]}.
#'
#' @slot label              Character scalar.  User-supplied label, used as the
#'   \code{ora_label} argument in
#'   \code{get_component_ora_result(obj, analysis, component_analysis, ora_label)}.
#' @slot parent_label       Character scalar.  Label of the parent
#'   \linkS4class{ComponentAnalysis}.
#' @slot grandparent_label  Character scalar.  Label of the parent
#'   \linkS4class{DifferentialPulldownAnalysis}.
#' @slot result_table       Data frame.  One row per (component, pathway)
#'   combination.  Columns:
#'   \code{component_id} (integer),
#'   \code{pathway} (character),
#'   \code{pval} (numeric),
#'   \code{padj} (numeric; BH-adjusted within each component),
#'   \code{overlap} (integer; genes in both the component and the pathway),
#'   \code{pathway_size} (integer; pathway size after filtering to universe),
#'   \code{n_genes_tested} (integer; number of genes tested for this
#'     component),
#'   \code{overlap_genes} (list column of character vectors; the overlapping
#'     gene IDs).
#'   Sorted by \code{component_id} then \code{padj} ascending.
#' @slot universe           Character vector.  The background gene set used.
#' @slot nodes_used         Character scalar.  \code{"all"} or \code{"seeds"}.
#' @slot go_category        Character scalar.  GO subcollection tested, or
#'   \code{NA} when a user-supplied pathway list was used.
#' @slot params             Named list.  Run parameters: \code{minSize},
#'   \code{maxSize}, \code{n_components_tested},
#'   \code{n_components_skipped} (too few genes to test),
#'   \code{universe_size}.
#' @slot created_at         \code{POSIXct} scalar.  Timestamp of the run.
#'
#' @exportClass ComponentOraResult
setClass(
  "ComponentOraResult",
  representation(
    label             = "character",
    parent_label      = "character",
    grandparent_label = "character",
    result_table      = "data.frame",
    universe          = "character",
    nodes_used        = "character",
    go_category       = "character",
    params            = "list",
    created_at        = "POSIXct"
  )
)


# ==============================================================================
#  show()
# ==============================================================================

#' @describeIn ComponentOraResult Print a concise summary.
#' @param object A \linkS4class{ComponentOraResult} object.
setMethod("show", "ComponentOraResult", function(object) {
  rt <- object@result_table

  n_comp_tested <- object@params$n_components_tested
  n_comp_skip   <- object@params$n_components_skipped
  n_sig         <- if (nrow(rt) > 0L) sum(rt$padj < 0.05, na.rm = TRUE) else 0L

  # Top hit per component
  top_per_comp <- if (nrow(rt) > 0L) {
    rt_sorted <- rt[order(rt$component_id, rt$padj), ]
    rt_sorted[!duplicated(rt_sorted$component_id), ]
  } else NULL

  cli::cli_rule("ComponentOraResult")
  cli::cli_text("Label           : {object@label}")
  cli::cli_text("Parent comp. run: {object@parent_label}")
  cli::cli_text("Diff. analysis  : {object@grandparent_label}")
  cli::cli_text("Created         : {format(object@created_at, '%Y-%m-%d %H:%M:%S')}")
  cli::cli_text("GO category     : {object@go_category}")
  cli::cli_text("Nodes used      : {object@nodes_used}")
  cli::cli_text("Universe size   : {object@params$universe_size} genes")
  cli::cli_text(
    "Components      : {n_comp_tested} tested  |  {n_comp_skip} skipped (too few genes)"
  )
  cli::cli_text(
    "Results         : {nrow(rt)} term(s) across all components  |  {n_sig} significant (padj < 0.05)"
  )

  if (!is.null(top_per_comp) && nrow(top_per_comp) > 0L) {
    cli::cli_text("Top term per component:")
    for (i in seq_len(min(5L, nrow(top_per_comp)))) {
      r <- top_per_comp[i, ]
      cli::cli_text(
        "  Comp {r$component_id}: {r$pathway}  (padj = {signif(r$padj, 3)}, overlap = {r$overlap})"
      )
    }
    if (nrow(top_per_comp) > 5L)
      cli::cli_text("  ... ({nrow(top_per_comp) - 5L} more components)")
  }

  cli::cli_text(
    "Params          : minSize = {object@params$minSize}  maxSize = {object@params$maxSize}"
  )
  cli::cli_rule()
  invisible(object)
})


# ==============================================================================
#  Exported: run_component_ora
# ==============================================================================

#' Run GO over-representation analysis on connected components
#'
#' Tests each connected component from a \linkS4class{ComponentAnalysis} run
#' against GO gene sets (or a user-supplied pathway list) using
#' \code{fgsea::fora()}.  The result annotates each component with the GO
#' terms most significantly over-represented among its member genes, providing
#' a biological interpretation of what each protein community is doing.
#'
#' @section Node-to-gene mapping:
#' Nodes in the component network are isoform IDs.  The reliable \code{gene}
#' vertex attribute on the IsoNetAnalyzeR network object is used to map each
#' isoform to its Ensembl gene ID.  This mapping is applied automatically;
#' no user input is required.
#'
#' @section Background universe:
#' The default background is all genes represented in the full IsoNetAnalyzeR
#' network.  This is the most principled choice for protein-interaction network
#' ORA because it restricts the background to proteins that could in principle
#' have appeared in any component.  Using all GO-annotated genes would inflate
#' apparent enrichment by including thousands of proteins absent from the
#' network.  You can override this with the \code{universe} argument.
#'
#' @section GO pathways:
#' When \code{pathways = NULL} (default), gene sets are fetched from MSigDB
#' via \code{msigdbr} using category \code{"C5"} (GO gene sets).
#' \code{go_category} selects the subcollection:
#' \describe{
#'   \item{\code{"GO:BP"}}{Biological process - \emph{what} the proteins do.
#'     Most informative for isoform-switch interpretation. \strong{Default.}}
#'   \item{\code{"GO:MF"}}{Molecular function - the biochemical activity.
#'     Useful for identifying domain-level functional differences between
#'     isoforms.}
#'   \item{\code{"GO:CC"}}{Cellular component - subcellular location.
#'     Less directly relevant to isoform switching; use as secondary context.}
#'   \item{\code{"ALL"}}{All three GO subcollections combined (excludes HPO).}
#' }
#'
#' @section Which nodes to test:
#' \describe{
#'   \item{\code{nodes = "all"} (default)}{All nodes in each component are
#'     tested.  Recommended when components are small (typical for
#'     differential pulldown results), as using all nodes maximises the gene
#'     set size and statistical power.}
#'   \item{\code{nodes = "seeds"}}{Only seed nodes (the experimentally
#'     regulated isoforms) are tested.  More conservative; useful when
#'     components are large and neighbours dominate.}
#' }
#'
#' @param obj                An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis           Character scalar.  Label of the parent
#'   \linkS4class{DifferentialPulldownAnalysis}.
#' @param component_analysis Character scalar.  Label of the
#'   \linkS4class{ComponentAnalysis} to annotate.
#' @param label              Character scalar.  Key under which the result is
#'   stored; retrieve it afterwards with
#'   \code{get_component_ora_result(obj, analysis, component_analysis, label)}.
#' @param pathways           Named list of character vectors (Ensembl gene IDs
#'   per pathway), or \code{NULL} (default) to fetch GO gene sets
#'   automatically via \code{msigdbr}.
#' @param go_category        Character scalar.  GO subcollection to fetch when
#'   \code{pathways = NULL}.  One of \code{"GO:BP"}, \code{"GO:MF"},
#'   \code{"GO:CC"}, or \code{"ALL"}.  Default \code{"GO:BP"}.
#' @param nodes              Character scalar.  Which nodes to test per
#'   component: \code{"all"} (default) or \code{"seeds"}.
#' @param universe           Character vector of Ensembl gene IDs to use as
#'   the ORA background, or \code{NULL} (default) to use all genes in the
#'   IsoNetAnalyzeR network.
#' @param minSize            Integer.  Minimum gene set size after filtering
#'   to the universe.  Default \code{10}.
#' @param maxSize            Integer.  Maximum gene set size.  Default
#'   \code{500}.
#'
#' @return The \linkS4class{IsoNetAnalyzeR} session object, invisibly, with
#'   the \linkS4class{ComponentOraResult} stored under \code{label} inside the
#'   component analysis.  Access via
#'   \code{get_component_ora_result(obj, analysis, component_analysis, label)}.
#'
#' @examples
#' \dontrun{
#' # Run component analysis first
#' obj <- run_component_analysis(obj, "my_run", label = "comp_both")
#'
#' # Annotate components with GO:BP terms (default)
#' obj <- run_component_ora(obj, "my_run", "comp_both", label = "go_bp")
#'
#' # Also test molecular function
#' obj <- run_component_ora(obj, "my_run", "comp_both", label = "go_mf",
#'                          go_category = "GO:MF")
#'
#' # show() summary of the result object
#' get_component_ora_result(obj, "my_run", "comp_both", "go_bp")
#'
#' # All significant terms for component 1
#' get_ora_table(obj, "my_run", "comp_both", "go_bp",
#'               component_id = 1, padj_cutoff = 0.05)
#'
#' # Top term per component across all components
#' get_ora_table(obj, "my_run", "comp_both", "go_bp", top_n = 1)
#' }
#'
#' @export
run_component_ora <- function(obj,
                              analysis,
                              component_analysis,
                              label,
                              pathways           = NULL,
                              go_category        = c("GO:BP", "GO:MF", "GO:CC", "ALL"),
                              nodes              = c("all", "seeds"),
                              universe           = NULL,
                              minSize            = 10L,
                              maxSize            = 500L) {

  # -- Package checks ----------------------------------------------------------
  if (!requireNamespace("fgsea", quietly = TRUE))
    stop("Package 'fgsea' is required. Install with:\n",
         "  BiocManager::install('fgsea')", call. = FALSE)

  # -- Input validation --------------------------------------------------------
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  diff_res <- .resolve_diff_analysis(obj, analysis)
  ca       <- .resolve_component_analysis(obj, analysis, component_analysis)

  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(trimws(label)))
    stop("`label` must be a non-empty character scalar.", call. = FALSE)

  go_category <- match.arg(go_category)
  nodes       <- match.arg(nodes)
  minSize <- .validate_positive_integer(minSize, "minSize")
  maxSize <- .validate_positive_integer(maxSize, "maxSize")
  if (maxSize < minSize)
    stop("`maxSize` must be greater than or equal to `minSize`.", call. = FALSE)

  if (nrow(ca@node_table) == 0L)
    stop(sprintf(
      "ComponentAnalysis '%s' has no nodes - nothing to test.", component_analysis
    ), call. = FALSE)

  # -- Build isoform -> gene ID mapping from the network -----------------------
  # The `gene` vertex attribute on the igraph network is the reliable mapping.
  message("Building isoform-to-gene mapping from network vertex attributes ...")
  iso_to_gene <- .build_iso_gene_map(obj@network)

  # -- Build universe ----------------------------------------------------------
  if (is.null(universe)) {
    universe <- .build_network_universe(obj@network, iso_to_gene)
    message(sprintf(
      "Universe: %d unique gene(s) represented in the full network.",
      length(universe)
    ))
  } else {
    if (!is.character(universe) || length(universe) == 0L)
      stop("`universe` must be a non-empty character vector of gene IDs.",
           call. = FALSE)
    message(sprintf("Using user-supplied universe: %d gene(s).", length(universe)))
  }

  # -- Fetch / validate pathways -----------------------------------------------
  if (is.null(pathways)) {
    pathways   <- .fetch_go_pathways(go_category)   # reuse helper from run_fgsea.R
    go_cat_str <- go_category
  } else {
    if (!is.list(pathways) || is.null(names(pathways)))
      stop("`pathways` must be a named list of character vectors.", call. = FALSE)
    go_cat_str <- NA_character_
    message(sprintf("Using user-supplied pathway list (%d gene set(s)).",
                    length(pathways)))
  }

  # -- Extract gene sets per component -----------------------------------------
  message(sprintf(
    "Extracting gene sets per component (nodes = '%s') ...", nodes
  ))

  component_genes <- .extract_component_genes(
    node_table  = ca@node_table,
    iso_to_gene = iso_to_gene,
    nodes       = nodes
  )

  # -- Run fora() per component ------------------------------------------------
  comp_ids     <- sort(unique(ca@node_table$component_id))
  n_tested     <- 0L
  n_skipped    <- 0L
  result_rows  <- list()

  message(sprintf(
    "Running fgsea::fora() on %d component(s) ...", length(comp_ids)
  ))

  for (cid in comp_ids) {
    genes_in_comp <- component_genes[[as.character(cid)]]

    # Intersect with universe (fora() requires query to be a subset of universe)
    genes_in_comp <- intersect(genes_in_comp, universe)

    if (length(genes_in_comp) == 0L) {
      message(sprintf(
        "  Component %d: 0 genes after universe intersection - skipping.", cid
      ))
      n_skipped <- n_skipped + 1L
      next
    }

    # fora() requires at least 1 gene; practical minimum for meaningful results
    # is minSize but we let fgsea handle that filtering internally
    fora_res <- tryCatch(
      fgsea::fora(
        pathways = pathways,
        genes    = genes_in_comp,
        universe = universe,
        minSize  = minSize,
        maxSize  = maxSize
      ),
      error = function(e) {
        message(sprintf(
          "  Component %d: fora() failed - %s", cid, conditionMessage(e)
        ))
        NULL
      }
    )

    if (is.null(fora_res) || nrow(fora_res) == 0L) {
      n_skipped <- n_skipped + 1L
      next
    }

    n_tested <- n_tested + 1L

    # BH correction within this component (fora returns raw pval only)
    fora_res$padj <- p.adjust(fora_res$pval, method = "BH")
    fora_res$padj <- round(fora_res$padj, 6L)
    fora_res$pval <- round(fora_res$pval, 6L)

    # Attach component ID and gene count
    fora_df <- data.frame(
      component_id    = cid,
      pathway         = fora_res$pathway,
      pval            = fora_res$pval,
      padj            = fora_res$padj,
      overlap         = fora_res$overlap,
      pathway_size    = fora_res$size,
      n_genes_tested  = length(genes_in_comp),
      stringsAsFactors = FALSE
    )

    # Preserve overlap_genes as a list column
    fora_df$overlap_genes <- fora_res$overlapGenes

    # Sort by padj within component
    fora_df <- fora_df[order(fora_df$padj, fora_df$pval), ]

    result_rows[[length(result_rows) + 1L]] <- fora_df
  }

  # -- Combine results ---------------------------------------------------------
  if (length(result_rows) == 0L) {
    message("No components produced ORA results - returning empty result.")
    result_tbl <- .empty_ora_table()
  } else {
    result_tbl <- do.call(rbind, result_rows)
    rownames(result_tbl) <- NULL
  }

  n_sig <- sum(result_tbl$padj < 0.05, na.rm = TRUE)
  message(sprintf(
    "Done. %d component(s) tested, %d skipped. %d significant term(s) at padj < 0.05.",
    n_tested, n_skipped, n_sig
  ))

  # -- Assemble result object --------------------------------------------------
  ora_result <- new(
    "ComponentOraResult",
    label             = label,
    parent_label      = component_analysis,
    grandparent_label = analysis,
    result_table      = result_tbl,
    universe          = universe,
    nodes_used        = nodes,
    go_category       = go_cat_str,
    params            = list(
      minSize              = minSize,
      maxSize              = maxSize,
      n_components_tested  = n_tested,
      n_components_skipped = n_skipped,
      universe_size        = length(universe)
    ),
    created_at        = Sys.time()
  )

  # -- Store and return --------------------------------------------------------
  # Path: obj -> analyses[[analysis]] -> component_analyses[[component_analysis]]
  #          -> ora_results[[label]]
  obj@analyses[[analysis]]@component_analyses[[component_analysis]]@ora_results[[label]] <-
    ora_result

  message(sprintf(paste0(
    "Done. Access results with:\n",
    "  get_component_ora_result(obj, \"%s\", \"%s\", \"%s\")"),
    analysis, component_analysis, label
  ))

  invisible(obj)
}


# ==============================================================================
#  Internal: build isoform -> gene mapping from network
# ==============================================================================

#' Build a named character vector mapping isoform IDs to gene IDs
#'
#' Uses the `gene` vertex attribute of the IsoNetAnalyzeR igraph network.
#' Vertices without a gene attribute (NA or empty) are silently dropped.
#'
#' @param network An igraph object with a `gene` vertex attribute.
#' @return Named character vector: names = isoform IDs, values = gene IDs.
#' @noRd
.build_iso_gene_map <- function(network) {
  vnames  <- igraph::V(network)$name
  vgenes  <- igraph::V(network)$gene

  if (is.null(vgenes))
    stop(
      "The network graph has no 'gene' vertex attribute. ",
      "Cannot map isoform IDs to gene IDs for ORA.",
      call. = FALSE
    )

  keep    <- !is.na(vgenes) & nchar(vgenes) > 0L
  mapping <- setNames(vgenes[keep], vnames[keep])
  mapping
}


# ==============================================================================
#  Internal: build network universe
# ==============================================================================

#' Extract all unique gene IDs from the network as the ORA background universe
#'
#' @param network    igraph object.
#' @param iso_to_gene Named character vector from .build_iso_gene_map().
#' @return Character vector of unique gene IDs.
#' @noRd
.build_network_universe <- function(network, iso_to_gene) {
  unique(iso_to_gene)
}


# ==============================================================================
#  Internal: extract gene sets per component
# ==============================================================================

#' Map component nodes to gene IDs for ORA input
#'
#' For each component, converts the isoform node IDs in the node_table to
#' gene IDs via iso_to_gene, optionally restricting to seed nodes only.
#' Returns a named list (component_id -> character vector of unique gene IDs).
#'
#' @param node_table  Data frame. ComponentAnalysis@@node_table.
#' @param iso_to_gene Named character vector from .build_iso_gene_map().
#' @param nodes       Character scalar. "all" or "seeds".
#' @return Named list of character vectors.
#' @noRd
.extract_component_genes <- function(node_table, iso_to_gene, nodes) {

  if (nodes == "seeds")
    node_table <- node_table[node_table$is_seed, , drop = FALSE]

  comp_ids <- sort(unique(node_table$component_id))

  result <- lapply(comp_ids, function(cid) {
    nodes_in_comp <- node_table$node[node_table$component_id == cid]

    # At gene level, node IDs are already gene IDs — mapping would produce
    # all NAs. Detect this by checking whether any node matches a known
    # isoform ID in the mapping.
    if (any(nodes_in_comp %in% names(iso_to_gene))) {
      # Isoform level: translate isoform IDs -> gene IDs
      genes <- iso_to_gene[nodes_in_comp]
      genes <- genes[!is.na(genes)]
      unique(as.character(genes))
    } else {
      # Gene level: nodes are already gene IDs
      unique(nodes_in_comp)
    }
  })

  names(result) <- as.character(comp_ids)
  result
}


# ==============================================================================
#  Internal: empty result table constructor
# ==============================================================================

#' @noRd
.empty_ora_table <- function() {
  df <- data.frame(
    component_id   = integer(0),
    pathway        = character(0),
    pval           = numeric(0),
    padj           = numeric(0),
    overlap        = integer(0),
    pathway_size   = integer(0),
    n_genes_tested = integer(0),
    stringsAsFactors = FALSE
  )
  df$overlap_genes <- list()
  df
}


# ==============================================================================
#  Exported: get_component_ora_result
# ==============================================================================

#' Retrieve a stored component ORA result
#'
#' @param obj                An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis           Character scalar.  Label of the parent
#'   \linkS4class{DifferentialPulldownAnalysis}.
#' @param component_analysis Character scalar.  Label of the
#'   \linkS4class{ComponentAnalysis}.
#' @param ora_label          Character scalar.  Label of the
#'   \linkS4class{ComponentOraResult}.
#'
#' @return A \linkS4class{ComponentOraResult} object.
#'
#' @examples
#' \dontrun{
#' ora <- get_component_ora_result(obj, "my_run", "comp_both", "go_bp")
#' ora  # show() summary
#'
#' # Significant terms for component 2
#' res <- ora@@result_table
#' res[res$component_id == 2 & res$padj < 0.05,
#'     c("pathway", "padj", "overlap", "pathway_size")]
#'
#' # Top GO term per component (one row per component)
#' res[!duplicated(res$component_id),
#'     c("component_id", "pathway", "padj", "overlap")]
#'
#' # Inspect overlapping genes for the top term of component 1
#' top1 <- res[res$component_id == 1, ]
#' top1$overlap_genes[[1]]
#' }
#'
#' @export
setGeneric("get_component_ora_result",
           function(obj, analysis, component_analysis, ora_label)
             standardGeneric("get_component_ora_result"))

#' @describeIn get_component_ora_result Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_component_ora_result", "IsoNetAnalyzeR",
          function(obj, analysis, component_analysis, ora_label) {

            ca <- .resolve_component_analysis(obj, analysis, component_analysis)

            if (!ora_label %in% names(ca@ora_results))
              stop(sprintf(
                "No ORA result labelled '%s' in component analysis '%s'. Available: %s",
                ora_label, component_analysis,
                paste(names(ca@ora_results), collapse = ", ")
              ), call. = FALSE)

            ca@ora_results[[ora_label]]
          })
