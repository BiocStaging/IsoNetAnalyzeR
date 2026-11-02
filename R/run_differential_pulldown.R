## run_differential_pulldown.R
## User-facing function to run a differential virtual pulldown analysis.
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: run_differential_pulldown

#' @include isonetanalyzer_class.R
NULL
#' @include diff_pulldown_gene_class.R
NULL
#' @include diff_pulldown_analysis_class.R
NULL

# ==============================================================================
#  Exported function
# ==============================================================================

#' Run a differential virtual pulldown analysis
#'
#' Reads a long-format table of caller-selected isoform seeds (one isoform per
#' row), runs a virtual pulldown for the up- and down-labelled seeds of each
#' gene separately, computes per-direction and per-gene statistics, and
#' adds a \linkS4class{DifferentialPulldownAnalysis} object to the
#' \linkS4class{IsoNetAnalyzeR} session object.
#'
#' @section Input table format:
#' A CSV, Excel (\code{.xlsx} / \code{.xls}), or already-loaded data frame in
#' \strong{long format} with exactly the following columns (any order):
#' \describe{
#'   \item{gene_id}{Gene identifier (e.g. Ensembl gene ID or gene symbol).}
#'   \item{isoform_id}{Isoform / transcript identifier matching vertex names
#'     in the session network.}
#'   \item{direction}{Must be exactly \code{"up"} or \code{"down"} (case
#'     insensitive; trimmed).}
#' }
#' Extra columns are retained in \code{@@input_table}. If \code{seed_padj_col}
#' is supplied, that column is used for retained-seed validation.
#'
#' @section Retained-seed support:
#' With \code{seed_padj_col} supplied, decay filtering and identifier resolution
#' occur first, then unmapped rows are removed. Each retained gene must have
#' both up and down seeds and at least one retained adjusted p-value below
#' \code{seed_fdr}. Other retained seeds may be non-significant or have missing
#' adjusted p-values. No maximum seed count is imposed. Supply only transcripts
#' passing the study's gene and effect-size filters; this function does not
#' perform RNA testing or choose an effect floor. Audits are stored in
#' \code{obj@@metadata$analysis_seed_support[[label]]}.
#' With \code{seed_padj_col = NULL}, the generic three-column exploratory
#' contract is unchanged and does not assert statistical support.
#'
#' @section Parallel execution:
#' Parallelism is opt-in and requires the \pkg{future} and \pkg{future.apply}
#' packages (\code{install.packages(c("future", "future.apply"))}).
#' When \code{parallel = TRUE} the per-gene pulldowns are dispatched via
#' \code{future.apply::future_lapply()}, which delegates all scheduling to
#' whatever \pkg{future} \emph{plan} the user has set before calling this
#' function.
#'
#' Typical setup:
#' \preformatted{
#' library(future)
#'
#' # Recommended: works on all platforms including Windows
#' plan(multisession, workers = 4)
#'
#' # Unix/macOS only: lower overhead via forking
#' plan(multicore, workers = 4)
#'
#' # Back to sequential when done
#' plan(sequential)
#' }
#'
#' If no plan has been set explicitly, \pkg{future} defaults to
#' \code{sequential}, so passing \code{parallel = TRUE} without configuring a
#' plan is safe but will not speed anything up.  An informative message is
#' emitted in that situation.
#'
#' The network graph embedded in \code{obj} is treated as read-only inside
#' each worker; no mutable state is shared between workers.
#'
#' @section Progress reporting:
#' A progress bar is always shown during the per-gene loop via \pkg{progressr}.
#' The bar works in both sequential and parallel modes.  Users can customise
#' the handler before calling this function, e.g.:
#' \preformatted{
#' progressr::handlers("cli")          # rich bar (default when cli is installed)
#' progressr::handlers("txtprogressbar") # plain fallback
#' }
#' If \code{progressr::handlers(global = TRUE)} has been set in the session,
#' that handler is used automatically; otherwise \code{with_progress()} applies
#' the default handler locally for this call only.
#'
#' @param obj      An \linkS4class{IsoNetAnalyzeR} session object.
#' @param label    Character scalar. Key under which to store the result;
#'   retrieve it afterwards with \code{get_analysis(obj, label)}.
#' @param input    One of:
#'   \itemize{
#'     \item A character scalar - path to a CSV or Excel file.
#'     \item A data frame already loaded into R.
#'   }
#' @param order    Integer scalar. Number of BFS hops for the virtual pulldown.
#' @param cutoff   Numeric scalar in \[0, 1\]. Internal-edge fraction cutoff.
#'   Defaults to \code{obj@@metadata[["network_cutoff"]]} when \code{NULL}.
#' @param parallel Logical. When \code{TRUE}, per-gene pulldowns are
#'   dispatched via \code{future.apply::future_lapply()} using the active
#'   \pkg{future} plan.  Requires the \pkg{future} and \pkg{future.apply}
#'   packages to be installed.  Default \code{FALSE}.
#' @param confirm  Logical. When \code{TRUE} (the default in interactive
#'   sessions), a network coverage report is printed before the pulldowns
#'   begin and the user is asked to confirm before proceeding.  Set to
#'   \code{FALSE} to skip the prompt in scripts or non-interactive contexts.
#'   The coverage report is always printed regardless of this flag; only the
#'   interactive pause is suppressed.
#' @param resolve_duplicates Logical. Resolve versioned and protein-identical
#'   redundant transcript identifiers before analysis. Same-direction
#'   representative collisions are collapsed and genes with opposite-direction
#'   collisions are excluded. Default \code{TRUE}; use \code{FALSE} to match
#'   identifiers exactly.
#' @param verbose  Logical. Print per-gene progress messages. Default
#'   \code{TRUE}.
#' @param seed_padj_col Optional name of the caller-supplied adjusted p-value
#'   column. When non-NULL, enforce significant retained support and both
#'   mapped directions after biological eligibility and identifier resolution.
#' @param seed_fdr Numeric scalar in (0, 1]. A retained seed is significant
#'   when its supplied adjusted p-value is strictly below this threshold.
#'   Default 0.05. No additional p-value adjustment is performed.
#'
#' @return An updated \linkS4class{IsoNetAnalyzeR} object (invisibly); access
#'   the result via \code{get_analysis(obj, label)}.
#'
#' @examples
#' seeds <- read.csv(system.file("extdata", "gse141432_t2d_vs_lean_seeds.csv",
#'                               package = "IsoNetAnalyzeR"))
#' cltcl1 <- seeds[seeds$gene_name == "CLTCL1", ]
#' obj <- IsoNetAnalyzeR("example")
#'
#' # Exploratory run: gene_id, isoform_id and direction only
#' obj <- run_differential_pulldown(obj, "explore", cltcl1, order = 1,
#'                                  confirm = FALSE, verbose = FALSE)
#'
#' # Require at least one significant retained seed per gene
#' obj <- run_differential_pulldown(
#'   obj, "supported", cltcl1, order = 1, confirm = FALSE, verbose = FALSE,
#'   seed_padj_col = "within_gene_padj", seed_fdr = 0.05
#' )
#' get_metadata(obj, "analysis_seed_support")$supported$gene_audit
#' get_cross_gene_summary(obj, "supported")
#'
#' @export
run_differential_pulldown <- function(obj,
                                      label,
                                      input,
                                      order,
                                      cutoff   = NULL,
                                      parallel = FALSE,
                                      confirm  = interactive(),
                                      verbose  = TRUE,
                                      resolve_duplicates = TRUE,
                                      seed_padj_col = NULL,
                                      seed_fdr = 0.05) {

  obj_name <- .obj_name(substitute(obj))

  # -- Input validation --------------------------------------------------------
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(trimws(label)))
    stop("`label` must be a single non-empty character string.", call. = FALSE)

  if (label %in% names(obj@analyses))
    warning(sprintf(
      "An analysis labelled '%s' already exists and will be overwritten.", label
    ))

  order <- .validate_positive_integer(order, "order")

  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel))
    stop("`parallel` must be a single non-NA logical value.", call. = FALSE)

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("`verbose` must be a single non-NA logical value.", call. = FALSE)

  if (!is.logical(resolve_duplicates) || length(resolve_duplicates) != 1L ||
      is.na(resolve_duplicates))
    stop("`resolve_duplicates` must be a single non-NA logical value.",
         call. = FALSE)

  # -- Check for future packages when parallel = TRUE -------------------------
  if (parallel) {
    if (!requireNamespace("future",       quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE))
      stop(
        "Parallel execution requires the 'future' and 'future.apply' packages.\n",
        "Install them with:\n",
        "  install.packages(c('future', 'future.apply'))\n",
        "Or run with parallel = FALSE for sequential execution.",
        call. = FALSE
      )

    # Inform the user if no parallel plan is active so the flag is not
    # silently a no-op.
    if (inherits(future::plan(), "sequential"))
      message(
        "parallel = TRUE but the active future plan is 'sequential'.\n",
        "Set a parallel plan before calling this function, e.g.:\n",
        "  future::plan(future::multisession, workers = 4)"
      )
  }

  # -- Cutoff resolution -------------------------------------------------------
  if (is.null(cutoff)) {
    cutoff <- obj@metadata[["network_cutoff"]]
    if (is.null(cutoff))
      stop(
        "`cutoff` not supplied and `obj@metadata$network_cutoff` is not set. ",
        "Please provide a numeric cutoff explicitly.",
        call. = FALSE
      )
    cutoff <- .validate_cutoff(cutoff, "obj@metadata$network_cutoff")
    if (verbose)
      message(sprintf("Using stored network_cutoff: %.4f", cutoff))
  } else {
    cutoff <- .validate_cutoff(cutoff)
    if (verbose)
      message(sprintf("Using user-supplied cutoff: %.4f", cutoff))
  }

  if (!is.logical(confirm) || length(confirm) != 1L || is.na(confirm))
    stop("`confirm` must be a single non-NA logical value.", call. = FALSE)

  # -- Load and validate input table -------------------------------------------
  tbl      <- .load_diff_input(input)
  tbl      <- .validate_diff_table(
    tbl, obj@network, resolve_duplicates = resolve_duplicates
  )
  duplicate_resolution <- attr(tbl, "dedup_resolution")
  biotype_exclusions <- attr(tbl, "biotype_exclusions")
  attr(tbl, "biotype_exclusions") <- NULL
  obj@metadata$analysis_biotype_exclusions[[label]] <- biotype_exclusions
  submitted_tbl <- if (is.null(duplicate_resolution)) tbl else
    duplicate_resolution$submitted_table
  attr(tbl, "dedup_resolution") <- NULL
  if (!is.null(seed_padj_col)) {
    support <- .filter_supported_seed_sets(tbl, obj@network, seed_padj_col, seed_fdr)
    obj@metadata$analysis_seed_support[[label]] <- support[c("gene_audit", "row_audit", "padj_column", "alpha")]
    tbl <- support$input
    if (!nrow(tbl)) stop(
      "No genes retain a significant mapped seed together with both up and down seeds.",
      call. = FALSE
    )
  }
  gene_ids <- unique(tbl$gene_id)

  # -- Network coverage report -------------------------------------------------
  # Coverage is computed against the full network vertex set.  Since
  # virtual_pulldown() now retains isolated vertices (isoforms with no
  # interact==1 edges) rather than deleting them, an isoform is "found" if it
  # exists anywhere in obj@network - even if it has no predicted interactions.
  # That case will produce a single-node subgraph (zero interactions), which
  # is a valid and interpretable result.
  network_vertices <- igraph::V(obj@network)$name
  pred_edge_ids <- which(igraph::edge_attr(obj@network, "interact") == 1L)
  interacting_vertices <- if (length(pred_edge_ids)) {
    unique(as.vector(igraph::ends(obj@network, pred_edge_ids, names = TRUE)))
  } else {
    character(0)
  }
  coverage <- .report_network_coverage(
    tbl, network_vertices, gene_ids, confirm, interacting_vertices,
    submitted_tbl = submitted_tbl,
    resolution_audit = if (is.null(duplicate_resolution)) NULL else
      duplicate_resolution$audit
  )

  if (verbose)
    message(sprintf(
      "Input table: %d rows, %d gene(s), %d unique isoform(s).  Mode: %s.",
      nrow(tbl),
      length(gene_ids),
      length(unique(tbl$isoform_id)),
      if (parallel) sprintf("parallel [%s]", class(future::plan())[[1L]])
      else "sequential"
    ))

  # -- Per-gene pulldown -------------------------------------------------------
  # Both branches call .process_one_gene(); only the dispatcher differs.
  # future_lapply() behaves identically to lapply() - the active future plan
  # decides whether work actually runs in parallel.
  #
  # Progress is reported via progressr::progressor(), which works transparently
  # in both sequential and parallel modes.  In parallel mode, worker signals are
  # relayed back to the main process automatically by the progressr framework.
  lapply_fn <- if (parallel) future.apply::future_lapply else base::lapply

  progressr::with_progress({
    p <- progressr::progressor(steps  = length(gene_ids))

    gene_results_list <- lapply_fn(
      gene_ids,
      FUN = function(gid) {
        p(message = sprintf(
          "Processing genes [%d%%] | %s",
          round(100 * which(gene_ids == gid) / length(gene_ids)),
          gid
        ))
        .process_one_gene(
          gid, tbl, obj, order, cutoff, verbose,
          network_vertices = network_vertices,
          submitted_tbl = submitted_tbl
        )
      }
    )
  })

  # -- Collect results; handle per-gene errors ---------------------------------
  # .process_one_gene() never stops: failures are returned as a
  # list(error = <character>) sentinel so that one bad gene does not abort
  # the whole run.
  gene_results <- stats::setNames(vector("list", length(gene_ids)), gene_ids)
  gene_errors  <- character(0)

  for (i in seq_along(gene_ids)) {
    gid <- gene_ids[[i]]
    res <- gene_results_list[[i]]

    if (is.list(res) && identical(names(res), "error")) {
      gene_errors[[gid]] <- res$error
      warning(sprintf(
        "Gene '%s': processing failed - %s", gid, res$error
      ), call. = FALSE)
    } else {
      gene_results[[gid]] <- res
    }
  }

  # Drop NULL slots (failed genes) before assembly
  gene_results <- Filter(Negate(is.null), gene_results)

  if (length(gene_results) == 0L && length(gene_errors) > 0L)
    stop(
      paste0(
        "Every gene failed during differential pulldown processing:\n",
        paste(sprintf("- %s: %s", names(gene_errors), gene_errors),
              collapse = "\n")
      ),
      call. = FALSE
    )

  # -- Assemble top-level result -----------------------------------------------
  result <- .new_diff_pulldown_analysis(
    label        = label,
    input_table  = tbl,
    gene_results = gene_results,
    order        = order,
    cutoff       = cutoff
  )

  if (verbose)
    message(sprintf("Done. DifferentialPulldownAnalysis '%s' created.", label))

  obj@analyses[[label]] <- result
  if (is.null(obj@metadata[["analysis_network_coverage"]]))
    obj@metadata[["analysis_network_coverage"]] <- list()
  obj@metadata[["analysis_network_coverage"]][[label]] <- coverage
  if (!is.null(duplicate_resolution)) {
    if (is.null(obj@metadata[["analysis_duplicate_resolution"]]))
      obj@metadata[["analysis_duplicate_resolution"]] <- list()
    obj@metadata[["analysis_duplicate_resolution"]][[label]] <-
      duplicate_resolution
  }
  message(sprintf("Done. Access the analysis with:\n  get_analysis(%s, \"%s\")",
                  obj_name, label))
  invisible(obj)
}


# ==============================================================================
#  Internal: per-gene unit of work
# ==============================================================================

#' Process a single gene in the differential pulldown loop
#'
#' Contains the per-gene body shared by the sequential and parallel code paths.
#' Errors are caught and returned as a \code{list(error = <character>)}
#' sentinel so that a failure in one gene does not abort a parallel job.
#'
#' @param gid              Character scalar. Gene ID to process.
#' @param tbl              Validated long-format input data frame.
#' @param obj              \linkS4class{IsoNetAnalyzeR} session object (read-only).
#' @param order            Integer. BFS hops.
#' @param cutoff           Numeric. Internal-edge fraction cutoff.
#' @param verbose          Logical. Emit per-gene progress messages.
#' @param network_vertices Character vector of all vertex names in the network.
#'
#' @return A \linkS4class{GeneDiffPulldown} object, or a
#'   \code{list(error = <character>)} sentinel on failure.
#' @noRd
.process_one_gene <- function(gid, tbl, obj, order, cutoff, verbose,
                              network_vertices = character(0),
                              submitted_tbl = tbl) {
  tryCatch({
    if (verbose) message(sprintf("  Processing gene: %s", gid))

    gene_tbl      <- tbl[tbl$gene_id == gid, , drop = FALSE]
    submitted_gene_tbl <- submitted_tbl[
      submitted_tbl$gene_id == gid, , drop = FALSE
    ]
    up_isoforms   <- gene_tbl$isoform_id[gene_tbl$direction == "up"]
    down_isoforms <- gene_tbl$isoform_id[gene_tbl$direction == "down"]

    up_graph   <- .safe_pulldown(obj, up_isoforms,   order, cutoff, gid, "up",   verbose)
    down_graph <- .safe_pulldown(obj, down_isoforms, order, cutoff, gid, "down", verbose)

    up_stats   <- .safe_analyze(up_graph)
    down_stats <- .safe_analyze(down_graph)

    .new_gene_diff_pulldown(
      gene_id          = gid,
      up_isoforms      = up_isoforms,
      down_isoforms    = down_isoforms,
      up_graph         = up_graph,
      down_graph       = down_graph,
      up_stats         = up_stats,
      down_stats       = down_stats,
      network_vertices = network_vertices,
      submitted_counts = c(
        up = length(unique(submitted_gene_tbl$isoform_id[
          submitted_gene_tbl$direction == "up"
        ])),
        down = length(unique(submitted_gene_tbl$isoform_id[
          submitted_gene_tbl$direction == "down"
        ]))
      )
    )
  },
  error = function(e) list(error = conditionMessage(e)))
}


# ==============================================================================
#  Internal: network coverage report
# ==============================================================================

#' Print a network coverage report and optionally ask the user to confirm
#'
#' Called once before the per-gene loop.  Computes how many isoforms and genes
#' from the input table are unmapped, mapped without a predicted edge, or
#' mapped with at least one \code{interact == 1} edge. Prints a formatted
#' summary via \pkg{cli} and
#' - when \code{confirm = TRUE} in an interactive session - asks the user
#' whether to proceed.
#'
#' @param tbl              Validated long-format input data frame.
#' @param network_vertices Character vector of all network vertex names.
#' @param gene_ids         Character vector of unique gene IDs in \code{tbl}.
#' @param confirm          Logical. Prompt the user before proceeding.
#' @param interacting_vertices Character vector of vertices incident to at
#'   least one edge with \code{interact == 1}.
#' @param min_coverage_warn Mapping fraction below which a warning is emitted.
#'
#' @return Invisibly, a typed list containing overall, per-direction, and gene
#'   coverage summaries. Stops if the user declines to proceed.
#' @noRd
.report_network_coverage <- function(tbl, network_vertices, gene_ids, confirm,
                                     interacting_vertices = character(0),
                                     min_coverage_warn = 0.10,
                                     submitted_tbl = tbl,
                                     resolution_audit = NULL) {

  if (is.null(resolution_audit)) {
    resolution_rows <- data.frame(
      submitted_id = submitted_tbl$isoform_id,
      resolved_id = submitted_tbl$isoform_id,
      gene_id = submitted_tbl$gene_id,
      direction = submitted_tbl$direction,
      resolution = "direct",
      action = "retained",
      stringsAsFactors = FALSE
    )
  } else {
    resolution_rows <- resolution_audit
  }

  classify <- function(rows, direction) {
    rows <- rows[!duplicated(rows$submitted_id), , drop = FALSE]
    ids <- rows$resolved_id
    mapped <- ids %in% network_vertices
    with_edge <- mapped & ids %in% interacting_vertices
    n_input <- nrow(rows)
    n_mapped <- sum(mapped)
    data.frame(
      direction = as.character(direction),
      n_input = as.integer(n_input),
      n_direct = as.integer(sum(rows$resolution == "direct")),
      n_version_normalized = as.integer(sum(
        rows$resolution == "version_normalized"
      )),
      n_remapped_from_duplicate = as.integer(sum(
        rows$resolution == "duplicate_remapped"
      )),
      n_resolution_incompatible = as.integer(sum(rows$resolution %in% c(
        "map_target_unavailable", "map_network_incompatible"
      ))),
      n_resolved_representatives = as.integer(length(unique(ids))),
      n_unmapped = as.integer(sum(!mapped)),
      n_mapped_no_predicted_edge = as.integer(sum(mapped & !with_edge)),
      n_mapped_with_predicted_edge = as.integer(sum(with_edge)),
      n_mapped = as.integer(n_mapped),
      pct_mapped = as.numeric(if (n_input) 100 * n_mapped / n_input else NA_real_),
      stringsAsFactors = FALSE
    )
  }

  # ---- isoform-level coverage ------------------------------------------------
  overall <- classify(resolution_rows, "overall")
  directions <- unique(resolution_rows$direction)
  by_direction <- do.call(rbind, lapply(directions, function(direction) {
    classify(
      resolution_rows[resolution_rows$direction == direction, , drop = FALSE],
      direction
    )
  }))
  rownames(by_direction) <- NULL

  # ---- gene-level coverage ---------------------------------------------------
  # A gene is "covered" if at least one of its isoforms is in the network
  submitted_gene_ids <- unique(submitted_tbl$gene_id)
  gene_covered <- vapply(submitted_gene_ids, function(g) {
    isos <- resolution_rows$resolved_id[resolution_rows$gene_id == g]
    any(isos %in% network_vertices)
  }, logical(1L))

  n_genes_total   <- length(submitted_gene_ids)
  n_genes_covered <- sum(gene_covered)
  pct_genes       <- if (n_genes_total > 0L)
    100 * n_genes_covered / n_genes_total else 0

  # ---- print report ----------------------------------------------------------
  cli::cli_rule("Network Coverage Report")
  .cli_aligned("")

  .cli_aligned("{.strong Isoforms}")
  .cli_aligned("  Unique submitted IDs            : {overall$n_input}")
  .cli_aligned("  Unique resolved representatives : {overall$n_resolved_representatives}")
  .cli_aligned("  Duplicate-remapped IDs          : {overall$n_remapped_from_duplicate}")
  .cli_aligned("  Version-normalized IDs          : {overall$n_version_normalized}")
  .cli_aligned("  Mapped, predicted edge          : {overall$n_mapped_with_predicted_edge}")
  .cli_aligned("  Mapped, no predicted edge       : {overall$n_mapped_no_predicted_edge}")
  .cli_aligned("  Unmapped                        : {overall$n_unmapped}")
  .cli_aligned("  Mapping                         : {overall$n_mapped} / {overall$n_input} ({round(overall$pct_mapped, 1)}%)")
  for (i in seq_len(nrow(by_direction))) {
    row <- by_direction[i, , drop = FALSE]
    .cli_aligned(
      "  {row$direction}: mapped with edge {row$n_mapped_with_predicted_edge}; mapped without edge {row$n_mapped_no_predicted_edge}; unmapped {row$n_unmapped} ({round(row$pct_mapped, 1)}% mapped)"
    )
  }
  .cli_aligned("")

  .cli_aligned("{.strong Genes}")
  .cli_aligned(
    "  Total in input                  : {n_genes_total}"
  )
  .cli_aligned(
    "  At least 1 isoform in network   : {n_genes_covered} / {n_genes_total} ({round(pct_genes, 1)}%)"
  )
  .cli_aligned("  No submitted isoforms mapped    : {n_genes_total - n_genes_covered}")
  .cli_aligned(
    "  Note: isoforms found in network with no interact==1 edges will produce a"
  )
  .cli_aligned(
    "        zero-interaction subgraph (single node), which is a valid result."
  )
  .cli_aligned("")
  cli::cli_rule()

  warn_rows <- rbind(overall, by_direction)
  warn_rows <- warn_rows[warn_rows$n_input > 0L &
                           warn_rows$pct_mapped < 100 * min_coverage_warn,
                         , drop = FALSE]
  for (i in seq_len(nrow(warn_rows))) {
    row <- warn_rows[i, , drop = FALSE]
    input_examples <- unique(resolution_rows$submitted_id[
      if (row$direction == "overall") rep(TRUE, nrow(resolution_rows)) else
        resolution_rows$direction == row$direction
    ])
    cli::cli_warn(c(
      "!" = "Only {round(row$pct_mapped, 1)}% of {row$direction} input isoforms were found in the network.",
      "i" = "This often indicates an annotation-release, version-suffix, or identifier-type mismatch.",
      "i" = "Network IDs look like: {paste(utils::head(network_vertices, 2L), collapse = ', ')}",
      "i" = "Input IDs look like: {paste(utils::head(input_examples, 2L), collapse = ', ')}"
    ))
  }

  # ---- optional confirmation prompt ------------------------------------------
  if (confirm && interactive()) {
    answer <- readline("Proceed with the analysis? [y/N] ")
    if (!tolower(trimws(answer)) %in% c("y", "yes"))
      stop("Analysis cancelled by user.", call. = FALSE)
  }

  gene_summary <- data.frame(
    n_genes_input = as.integer(n_genes_total),
    n_genes_analyzed = as.integer(length(gene_ids)),
    n_genes_dropped_opposite_collision = as.integer(
      length(setdiff(submitted_gene_ids, gene_ids))
    ),
    n_genes_with_mapped_isoform = as.integer(n_genes_covered),
    n_genes_without_mapped_isoform = as.integer(n_genes_total - n_genes_covered),
    pct_genes_with_mapped_isoform = as.numeric(pct_genes),
    stringsAsFactors = FALSE
  )

  invisible(list(
    overall = overall,
    by_direction = by_direction,
    genes = gene_summary,
    warning_threshold = as.numeric(min_coverage_warn)
  ))
}

# Require evidence on mapped rows remaining after biotype and ID resolution.
.filter_supported_seed_sets <- function(tbl, network, padj_column, alpha) {
  if (!is.character(padj_column) || length(padj_column) != 1L ||
      is.na(padj_column) || !nzchar(padj_column) || !padj_column %in% names(tbl))
    stop("`seed_padj_col` must name an adjusted p-value column in the input.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha > 1)
    stop("`seed_fdr` must be a finite number in (0, 1].", call. = FALSE)
  p <- tbl[[padj_column]]
  if (!is.numeric(p) || any(!is.na(p) & (!is.finite(p) | p < 0 | p > 1)))
    stop("The seed adjusted p-value column must be numeric in [0, 1] or NA.", call. = FALSE)
  mapped <- tbl$isoform_id %in% igraph::V(network)$name
  significant <- !is.na(p) & p < alpha
  genes <- unique(tbl$gene_id)
  audit <- do.call(rbind, lapply(genes, function(g) {
    rows <- tbl$gene_id == g & mapped
    up <- any(rows & tbl$direction == "up")
    down <- any(rows & tbl$direction == "down")
    sig <- sum(rows & significant)
    data.frame(gene_id = g, has_up = up, has_down = down,
      n_significant_retained_candidates = sig,
      retained = up && down && sig > 0L,
      reason = if (!up || !down) "not_bidirectional_after_mapping"
        else if (!sig) "no_significant_retained_seed" else "retained",
      stringsAsFactors = FALSE)
  }))
  keep <- mapped & tbl$gene_id %in% audit$gene_id[audit$retained]
  row_audit <- data.frame(gene_id = tbl$gene_id, isoform_id = tbl$isoform_id,
    direction = tbl$direction, adjusted_pvalue = p, mapped = mapped,
    significant = significant, retained = keep,
    reason = ifelse(!mapped, "transcript_not_in_session_network",
      audit$reason[match(tbl$gene_id, audit$gene_id)]), stringsAsFactors = FALSE)
  list(input = tbl[keep, , drop = FALSE], gene_audit = audit,
       row_audit = row_audit, padj_column = padj_column, alpha = alpha)
}

#' Load the differential input from a file path or data frame
#' @param input Character scalar (file path) or data frame.
#' @return data.frame.
#' @noRd
.load_diff_input <- function(input) {

  if (is.data.frame(input))
    return(input)

  if (!is.character(input) || length(input) != 1L)
    stop(
      "`input` must be a file path (character scalar) or a data frame.",
      call. = FALSE
    )

  if (!file.exists(input))
    stop(sprintf("File not found: '%s'", input), call. = FALSE)

  ext <- tolower(tools::file_ext(input))

  if (ext == "csv")
    return(utils::read.csv(input, stringsAsFactors = FALSE))

  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop(
        "Package 'readxl' is required to read Excel files: ",
        "install.packages('readxl').",
        call. = FALSE
      )
    return(as.data.frame(readxl::read_excel(input)))
  }

  stop(sprintf(
    "Unsupported file extension '.%s'. Supply a CSV or Excel (.xlsx/.xls) file.",
    ext
  ), call. = FALSE)
}


# ==============================================================================
#  Internal: input validation
# ==============================================================================

#' Validate and normalise the differential input table
#'
#' Checks required columns exist, trims whitespace, lower-cases direction,
#' drops rows with blank gene/isoform IDs, and warns about unrecognised
#' direction values.
#'
#' @param tbl data.frame.
#' @param network Optional igraph used for identifier and gene checks.
#' @param resolve_duplicates Logical; apply packaged duplicate resolution.
#' @param duplicate_map Optional validated map used by tests/custom workflows.
#' @return Cleaned data.frame.
#' @noRd
.validate_diff_table <- function(tbl, network = NULL,
                                 resolve_duplicates = FALSE,
                                 duplicate_map = NULL) {

  required <- c("gene_id", "isoform_id", "direction")
  missing  <- setdiff(required, names(tbl))

  if (length(missing))
    stop(sprintf(
      "Input table is missing required column(s): %s",
      paste(missing, collapse = ", ")
    ), call. = FALSE)

  if (!is.logical(resolve_duplicates) || length(resolve_duplicates) != 1L ||
      is.na(resolve_duplicates))
    stop("`resolve_duplicates` must be a single non-NA logical value.",
         call. = FALSE)

  tbl$gene_id    <- trimws(as.character(tbl$gene_id))
  tbl$isoform_id <- trimws(as.character(tbl$isoform_id))
  tbl$direction  <- tolower(trimws(as.character(tbl$direction)))

  blank_gene    <- is.na(tbl$gene_id) | !nzchar(tbl$gene_id)
  blank_isoform <- is.na(tbl$isoform_id) | !nzchar(tbl$isoform_id)

  if (any(blank_gene | blank_isoform)) {
    n_drop <- sum(blank_gene | blank_isoform)
    warning(sprintf(
      "Dropping %d row(s) with missing or empty gene_id or isoform_id.", n_drop
    ), call. = FALSE)
    tbl <- tbl[!(blank_gene | blank_isoform), , drop = FALSE]
  }

  bad_dirs <- is.na(tbl$direction) | !tbl$direction %in% c("up", "down")
  if (any(bad_dirs))
    stop(sprintf(
      "Unrecognised direction value(s): %s. Must be 'up' or 'down'.",
      paste(unique(ifelse(is.na(tbl$direction[bad_dirs]), "<missing>",
                          tbl$direction[bad_dirs])), collapse = ", ")
    ), call. = FALSE)

  if (nrow(tbl) == 0L)
    stop("Input table is empty after cleaning.", call. = FALSE)

  duplicate_rows <- duplicated(tbl[required])
  if (any(duplicate_rows)) {
    warning(sprintf(
      "Removing %d exact duplicate input row(s).", sum(duplicate_rows)
    ), call. = FALSE)
    tbl <- tbl[!duplicate_rows, , drop = FALSE]
  }

  isoform_directions <- tapply(tbl$direction, tbl$isoform_id,
                               function(x) length(unique(x)))
  direction_conflicts <- names(isoform_directions)[isoform_directions > 1L]
  if (length(direction_conflicts))
    stop(sprintf(
      "Each isoform must have one direction; conflicting isoform(s): %s.",
      paste(direction_conflicts, collapse = ", ")
    ), call. = FALSE)

  isoform_genes <- tapply(tbl$gene_id, tbl$isoform_id,
                          function(x) length(unique(x)))
  gene_conflicts <- names(isoform_genes)[isoform_genes > 1L]
  if (length(gene_conflicts))
    stop(sprintf(
      "Each isoform must belong to one input gene; conflicting isoform(s): %s.",
      paste(gene_conflicts, collapse = ", ")
    ), call. = FALSE)

  excluded_ids <- .excluded_decay_ids(network)
  excluded <- sub("\\.[0-9]+$", "", tbl$isoform_id) %in% excluded_ids
  biotype_exclusions <- tbl[excluded, , drop = FALSE]
  if (any(excluded)) {
    message(sprintf("Excluded %d submitted seed row(s) by transcript-biotype policy.", sum(excluded)))
    tbl <- tbl[!excluded, , drop = FALSE]
    if (!nrow(tbl)) stop("No seeds remain after transcript-biotype exclusion.", call. = FALSE)
  }

  if (resolve_duplicates) {
    submitted_tbl <- tbl
    resolved_tbl <- resolve_duplicate_isoforms(
      tbl, network = network, duplicate_map = duplicate_map
    )
    audit <- attr(resolved_tbl, "dedup_audit")
    attr(resolved_tbl, "dedup_audit") <- NULL
    audit$gene_id <- submitted_tbl$gene_id
    audit$direction <- submitted_tbl$direction
    audit$action <- "retained"
    tbl <- resolved_tbl

    collision_key <- paste(tbl$gene_id, tbl$isoform_id, sep = "\r")
    submitted_per_key <- tapply(
      audit$submitted_id, collision_key,
      function(x) length(unique(x))
    )
    collision_keys <- names(submitted_per_key)[submitted_per_key > 1L]
    drop_genes <- character(0)
    collapse_rows <- integer(0)

    for (key in collision_keys) {
      rows <- which(collision_key == key)
      directions <- unique(tbl$direction[rows])
      representative <- unique(tbl$isoform_id[rows])[[1L]]
      gene <- unique(tbl$gene_id[rows])[[1L]]
      submitted_description <- paste(
        paste0(audit$submitted_id[rows], " (", tbl$direction[rows], ")"),
        collapse = ", "
      )
      if (length(directions) > 1L) {
        drop_genes <- c(drop_genes, gene)
        warning(
          sprintf("Dropping gene '%s': protein-identical submitted transcripts ",
                  gene),
          sprintf("%s resolve to representative '%s' with opposite directions.",
                  submitted_description, representative),
          call. = FALSE
        )
      } else {
        collapse_rows <- c(collapse_rows, rows[-1L])
        audit$action[rows[-1L]] <- "collapsed_same_direction"
        message(
          sprintf("Collapsed %d submitted isoform IDs into one protein-identical ",
                  length(unique(audit$submitted_id[rows]))),
          sprintf("network representative '%s' for gene '%s'.",
                  representative, gene)
        )
      }
    }

    if (length(drop_genes)) {
      drop_genes <- unique(drop_genes)
      audit$action[audit$gene_id %in% drop_genes] <-
        "dropped_opposite_direction"
    }
    keep <- !(tbl$gene_id %in% drop_genes)
    if (length(collapse_rows)) keep[collapse_rows] <- FALSE
    tbl <- tbl[keep, , drop = FALSE]

    if (nrow(tbl) == 0L) {
      stop(
        "No analyzable rows remain after protein-identical collision ",
        "handling. All submitted genes contained opposite-direction ",
        "collisions.",
        call. = FALSE
      )
    }

    resolution_summary <- data.frame(
      n_submitted_rows = nrow(submitted_tbl),
      n_resolved_rows = nrow(tbl),
      n_duplicate_remapped = sum(audit$resolution == "duplicate_remapped"),
      n_version_normalized = sum(audit$resolution == "version_normalized"),
      n_same_direction_collapsed = sum(
        audit$action == "collapsed_same_direction"
      ),
      n_opposite_direction_genes_dropped = length(unique(
        audit$gene_id[audit$action == "dropped_opposite_direction"]
      )),
      stringsAsFactors = FALSE
    )
    attr(tbl, "dedup_resolution") <- list(
      submitted_table = submitted_tbl,
      audit = audit,
      summary = resolution_summary
    )
  }

  if (!is.null(network)) {
    network_gene <- stats::setNames(
      as.character(igraph::V(network)$gene), igraph::V(network)$name
    )
    mapped <- tbl$isoform_id %in% names(network_gene)
    mismatched <- mapped & tbl$gene_id != unname(network_gene[tbl$isoform_id])
    if (any(mismatched))
      stop(sprintf(
        paste0(
          "Input transcript-to-gene assignment disagrees with the network for: %s. ",
          "Check identifiers and annotation release."
        ),
        paste(unique(tbl$isoform_id[mismatched]), collapse = ", ")
      ), call. = FALSE)
  }

  if (nrow(biotype_exclusions)) attr(tbl, "biotype_exclusions") <- biotype_exclusions
  tbl
}


# ==============================================================================
#  Internal: safe pulldown wrapper
# ==============================================================================

#' Run virtual_pulldown, returning NULL instead of stopping on failure
#'
#' @param obj       IsoNetAnalyzeR object.
#' @param isoforms  Character vector of seed isoform IDs.
#' @param order     Integer.
#' @param cutoff    Numeric.
#' @param gene_id   Character scalar - used in warning messages only.
#' @param direction Character scalar - "up" or "down".
#' @param verbose   Logical.
#'
#' @return igraph subgraph, or NULL on failure.
#' @noRd
.safe_pulldown <- function(obj, isoforms, order, cutoff,
                           gene_id, direction, verbose) {
  if (length(isoforms) == 0L) {
    if (verbose)
      message(sprintf(
        "    [%s] No %s-regulated isoforms - skipping.", gene_id, direction
      ))
    return(NULL)
  }

  tryCatch(
    virtual_pulldown(
      obj           = obj,
      seed_vertices = isoforms,
      order         = order,
      cutoff        = cutoff,
      resolve_duplicates = FALSE
    ),
    error = function(e) {
      warning(sprintf(
        "Gene '%s' (%s): virtual_pulldown failed - %s",
        gene_id, direction, conditionMessage(e)
      ), call. = FALSE)
      NULL
    }
  )
}


# ==============================================================================
#  Internal: safe analyze_subgraph wrapper
# ==============================================================================

#' Run analyze_subgraph, returning a zero-row data frame on NULL input
#'
#' @param g igraph or NULL.
#' @return data.frame (may have zero rows).
#' @noRd
.safe_analyze <- function(g) {
  if (is.null(g))
    return(data.frame(
      isoform            = character(0),
      is_seed            = logical(0),
      degree             = integer(0),
      clustering         = numeric(0),
      distance_from_seed = numeric(0),
      mean_pi            = numeric(0),
      pi_entropy         = numeric(0),
      stringsAsFactors   = FALSE
    ))
  analyze_subgraph(g)
}
