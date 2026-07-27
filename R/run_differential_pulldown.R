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
#' Reads a long-format table of differentially expressed isoforms (one isoform
#' per row), runs a virtual pulldown for the up- and down-regulated isoforms
#' of each gene separately, computes per-direction and per-gene statistics, and
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
#' Extra columns are silently retained in \code{@@input_table} but are not
#' used by the analysis.
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
#' @param verbose  Logical. Print per-gene progress messages. Default
#'   \code{TRUE}.
#'
#' @return An updated \linkS4class{IsoNetAnalyzeR} object (invisibly); access
#'   the result via \code{get_analysis(obj, label)}.
#'
#' @examples
#' \dontrun{
#' obj <- IsoNetAnalyzeR("demo")
#'
#' # --- Sequential (default, no extra packages needed) -----------------------
#' obj <- run_differential_pulldown(
#'   obj, "condition_A_vs_B", "de_isoforms.csv", order = 1
#' )
#'
#' # --- Parallel (requires future + future.apply) ----------------------------
#' library(future)
#' plan(multisession, workers = 4)   # works on all platforms
#'
#' obj <- run_differential_pulldown(
#'   obj, "condition_A_vs_B", "de_isoforms.csv",
#'   order    = 1,
#'   parallel = TRUE
#' )
#'
#' plan(sequential)   # restore after the run
#'
#' # --- Inspect results ------------------------------------------------------
#' get_analysis(obj, "condition_A_vs_B")
#' get_gene_result(obj, "condition_A_vs_B", "ENSG00000012048")
#' plot_network(get_gene_result(obj, "condition_A_vs_B", "ENSG00000012048"))
#' }
#'
#' @export
run_differential_pulldown <- function(obj,
                                      label,
                                      input,
                                      order,
                                      cutoff   = NULL,
                                      parallel = FALSE,
                                      confirm  = interactive(),
                                      verbose  = TRUE) {

  # -- Input validation --------------------------------------------------------
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.character(label) || length(label) != 1L || nchar(trimws(label)) == 0L)
    stop("`label` must be a single non-empty character string.", call. = FALSE)

  if (label %in% names(obj@analyses))
    warning(sprintf(
      "An analysis labelled '%s' already exists and will be overwritten.", label
    ))

  if (!is.numeric(order) || length(order) != 1L || order < 1L)
    stop("`order` must be a positive integer scalar.", call. = FALSE)
  order <- as.integer(order)

  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel))
    stop("`parallel` must be a single non-NA logical value.", call. = FALSE)

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
    if (is.null(cutoff) || is.na(cutoff))
      stop(
        "`cutoff` not supplied and `obj@metadata$network_cutoff` is not set. ",
        "Please provide a numeric cutoff explicitly.",
        call. = FALSE
      )
    if (verbose)
      message(sprintf("Using stored network_cutoff: %.4f", cutoff))
  } else {
    if (!is.numeric(cutoff) || length(cutoff) != 1L ||
        cutoff < 0 || cutoff > 1)
      stop("`cutoff` must be a single numeric value in [0, 1].", call. = FALSE)
    if (verbose)
      message(sprintf("Using user-supplied cutoff: %.4f", cutoff))
  }

  if (!is.logical(confirm) || length(confirm) != 1L || is.na(confirm))
    stop("`confirm` must be a single non-NA logical value.", call. = FALSE)

  # -- Load and validate input table -------------------------------------------
  tbl      <- .load_diff_input(input)
  tbl      <- .validate_diff_table(tbl)
  gene_ids <- unique(tbl$gene_id)

  # -- Network coverage report -------------------------------------------------
  # Coverage is computed against the full network vertex set.  Since
  # virtual_pulldown() now retains isolated vertices (isoforms with no
  # interact==1 edges) rather than deleting them, an isoform is "found" if it
  # exists anywhere in obj@network - even if it has no predicted interactions.
  # That case will produce a single-node subgraph (zero interactions), which
  # is a valid and interpretable result.
  network_vertices <- igraph::V(obj@network)$name
  .report_network_coverage(tbl, network_vertices, gene_ids, confirm)

  if (verbose)
    message(sprintf(
      "Input table: %d rows, %d gene(s), %d unique isoform(s).  Mode: %s.",
      nrow(tbl),
      length(gene_ids),
      length(unique(tbl$isoform_id)),
      if (parallel) paste0("parallel [", class(future::plan())[[1L]], "]")
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
        .process_one_gene(gid, tbl, obj, order, cutoff, verbose,
                          network_vertices = network_vertices)
      }
    )
  })

  # -- Collect results; handle per-gene errors ---------------------------------
  # .process_one_gene() never stops: failures are returned as a
  # list(error = <character>) sentinel so that one bad gene does not abort
  # the whole run.
  gene_results <- stats::setNames(vector("list", length(gene_ids)), gene_ids)

  for (i in seq_along(gene_ids)) {
    gid <- gene_ids[[i]]
    res <- gene_results_list[[i]]

    if (is.list(res) && identical(names(res), "error")) {
      warning(sprintf(
        "Gene '%s': processing failed - %s", gid, res$error
      ), call. = FALSE)
    } else {
      gene_results[[gid]] <- res
    }
  }

  # Drop NULL slots (failed genes) before assembly
  gene_results <- Filter(Negate(is.null), gene_results)

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
  message(sprintf("Done. Access the analysis with:\n  get_analysis(obj, \"%s\")", label))
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
                              network_vertices = character(0)) {
  tryCatch({
    if (verbose) message(sprintf("  Processing gene: %s", gid))

    gene_tbl      <- tbl[tbl$gene_id == gid, , drop = FALSE]
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
      network_vertices = network_vertices
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
#' from the input table are represented in the \emph{interaction-filtered}
#' subgraph (edges with \code{interact == 1}), which is the graph actually used
#' by \code{virtual_pulldown()}.  Prints a formatted summary via \pkg{cli} and
#' - when \code{confirm = TRUE} in an interactive session - asks the user
#' whether to proceed.
#'
#' @param tbl              Validated long-format input data frame.
#' @param network_vertices Character vector of vertex names in the
#'   \emph{interaction-filtered} subgraph (passed in from the main function
#'   body, which constructs it once to avoid repeating the edge filter).
#' @param gene_ids         Character vector of unique gene IDs in \code{tbl}.
#' @param confirm          Logical. Prompt the user before proceeding.
#'
#' @return Invisibly \code{NULL}.  Stops with an informative message if the
#'   user declines to proceed.
#' @noRd
.report_network_coverage <- function(tbl, network_vertices, gene_ids, confirm) {

  # ---- isoform-level coverage ------------------------------------------------
  all_isoforms     <- unique(tbl$isoform_id)
  found_isoforms   <- all_isoforms[all_isoforms %in% network_vertices]
  n_iso_total      <- length(all_isoforms)
  n_iso_found      <- length(found_isoforms)
  pct_iso          <- if (n_iso_total > 0L) 100 * n_iso_found / n_iso_total else 0

  up_isoforms_all   <- unique(tbl$isoform_id[tbl$direction == "up"])
  down_isoforms_all <- unique(tbl$isoform_id[tbl$direction == "down"])
  n_up_found_all   <- sum(up_isoforms_all   %in% network_vertices)
  n_down_found_all <- sum(down_isoforms_all %in% network_vertices)

  pct_up           <- if (length(up_isoforms_all)   > 0L)
    100 * n_up_found_all   / length(up_isoforms_all)   else 0
  pct_down         <- if (length(down_isoforms_all) > 0L)
    100 * n_down_found_all / length(down_isoforms_all) else 0

  # ---- gene-level coverage ---------------------------------------------------
  # A gene is "covered" if at least one of its isoforms is in the network
  gene_covered <- vapply(gene_ids, function(g) {
    isos <- tbl$isoform_id[tbl$gene_id == g]
    any(isos %in% network_vertices)
  }, logical(1L))

  n_genes_total   <- length(gene_ids)
  n_genes_covered <- sum(gene_covered)
  pct_genes       <- if (n_genes_total > 0L)
    100 * n_genes_covered / n_genes_total else 0

  # ---- print report ----------------------------------------------------------
  cli::cli_rule("Network Coverage Report")
  cli::cli_text("")

  cli::cli_text("{.strong Isoforms}")
  cli::cli_text(
    "  Total in input   : {n_iso_total}"
  )
  cli::cli_text(
    "  Found in network : {n_iso_found} / {n_iso_total} ({round(pct_iso, 1)}%)"
  )
  cli::cli_text(
    "  Up-regulated     : {n_up_found_all} / {length(up_isoforms_all)} found ({round(pct_up, 1)}%)"
  )
  cli::cli_text(
    "  Down-regulated   : {n_down_found_all} / {length(down_isoforms_all)} found ({round(pct_down, 1)}%)"
  )
  cli::cli_text("")

  cli::cli_text("{.strong Genes}")
  cli::cli_text(
    "  Total in input   : {n_genes_total}"
  )
  cli::cli_text(
    "  At least 1 isoform in network : {n_genes_covered} / {n_genes_total} ({round(pct_genes, 1)}%)"
  )
  cli::cli_text(
    "  No isoforms in network        : {n_genes_total - n_genes_covered} gene(s) - these will produce no subgraph"
  )
  cli::cli_text(
    "  Note: isoforms found in network with no interact==1 edges will produce a"
  )
  cli::cli_text(
    "        zero-interaction subgraph (single node), which is a valid result."
  )
  cli::cli_text("")
  cli::cli_rule()

  # ---- optional confirmation prompt ------------------------------------------
  if (confirm && interactive()) {
    answer <- readline("Proceed with the analysis? [y/N] ")
    if (!tolower(trimws(answer)) %in% c("y", "yes"))
      stop("Analysis cancelled by user.", call. = FALSE)
  }

  invisible(NULL)
}

#' Load the differential input from a file path or data frame
#'
#' Supports CSV and Excel (\code{.xlsx} / \code{.xls}) files.  Returns a data
#' frame; does not validate column names (that is done separately).
#'
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
#' @return Cleaned data.frame.
#' @noRd
.validate_diff_table <- function(tbl) {

  required <- c("gene_id", "isoform_id", "direction")
  missing  <- setdiff(required, names(tbl))

  if (length(missing))
    stop(sprintf(
      "Input table is missing required column(s): %s",
      paste(missing, collapse = ", ")
    ), call. = FALSE)

  tbl$gene_id    <- trimws(as.character(tbl$gene_id))
  tbl$isoform_id <- trimws(as.character(tbl$isoform_id))
  tbl$direction  <- tolower(trimws(as.character(tbl$direction)))

  blank_gene    <- nchar(tbl$gene_id)    == 0L
  blank_isoform <- nchar(tbl$isoform_id) == 0L

  if (any(blank_gene | blank_isoform)) {
    n_drop <- sum(blank_gene | blank_isoform)
    warning(sprintf(
      "Dropping %d row(s) with empty gene_id or isoform_id.", n_drop
    ), call. = FALSE)
    tbl <- tbl[!(blank_gene | blank_isoform), , drop = FALSE]
  }

  bad_dirs <- !tbl$direction %in% c("up", "down")
  if (any(bad_dirs))
    stop(sprintf(
      "Unrecognised direction value(s): %s. Must be 'up' or 'down'.",
      paste(unique(tbl$direction[bad_dirs]), collapse = ", ")
    ), call. = FALSE)

  if (nrow(tbl) == 0L)
    stop("Input table is empty after cleaning.", call. = FALSE)

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
      cutoff        = cutoff
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
