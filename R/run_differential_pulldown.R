## run_differential_pulldown.R
## User-facing function to run a differential virtual pulldown analysis.
##
## Internal helpers (prefixed with ".") are not exported.
## Exported: run_differential_pulldown


# ==============================================================================
#  Exported function
# ==============================================================================

#' Run a differential virtual pulldown analysis
#'
#' Reads a long-format table of differentially expressed isoforms (one isoform
#' per row), runs a virtual pulldown for the up- and down-regulated isoforms
#' of each gene separately, computes per-direction and per-gene statistics, and
#' returns a \linkS4class{DifferentialPulldownAnalysis} object.  The result is
#' optionally stored in \code{obj@@analyses} when \code{store = TRUE}.
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
#' @param obj     An \linkS4class{IsoNetAnalyzeR} session object.
#' @param input   One of:
#'   \itemize{
#'     \item A character scalar — path to a CSV or Excel file.
#'     \item A data frame already loaded into R.
#'   }
#' @param order   Integer scalar. Number of BFS hops for the virtual pulldown.
#' @param cutoff  Numeric scalar in \[0, 1\].  Internal-edge fraction cutoff.
#'   Defaults to \code{obj@@metadata[["network_cutoff"]]} when \code{NULL}.
#' @param label   Character scalar. Key used when \code{store = TRUE}.
#'   Defaults to \code{"diff_pulldown_<timestamp>"}.
#' @param verbose Logical. Print per-gene progress messages.  Default
#'   \code{TRUE}.
#'
#' @return An updated \linkS4class{IsoNetAnalyzeR}
#'   object (invisibly); access the result via
#'   \code{obj@@analyses[[label]]}.
#'
#' @examples
#' \dontrun{
#' obj <- IsoNetAnalyzeR("demo")
#'
#' # From a CSV file
#' obj <- run_differential_pulldown(obj,
#'                                  "condition_A_vs_B",
#'                                  "de_isoforms.csv",
#'                                  order = 1)
#'
#' # From an already-loaded data frame, stored in the session object
#' obj <- run_differential_pulldown(
#'   obj, "condition_A_vs_B", de_df, order = 1
#' )
#' obj@@analyses$condition_A_vs_B
#'
#' # Inspect a single gene
#' obj@@analyses$condition_A_vs_B@@gene_results[["ENSG00000012048"]]
#'
#' # Plot both directions for one gene
#' plot_network(obj@@analyses$condition_A_vs_B@@gene_results[["ENSG00000012048"]])
#' }
#'
#' @export
run_differential_pulldown <- function(obj,
                                      label,
                                      input,
                                      order,
                                      cutoff  = NULL,
                                      verbose = TRUE) {

  # ── Input validation — session object ───────────────────────────────────────
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.numeric(order) || length(order) != 1L || order < 1L)
    stop("`order` must be a positive integer scalar.", call. = FALSE)
  order <- as.integer(order)

  # ── Cutoff resolution ────────────────────────────────────────────────────────
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

  # ── Load and validate input table ────────────────────────────────────────────
  tbl <- .load_diff_input(input)
  tbl <- .validate_diff_table(tbl)

  if (verbose)
    message(sprintf(
      "Input table: %d rows, %d gene(s), %d unique isoform(s).",
      nrow(tbl),
      length(unique(tbl$gene_id)),
      length(unique(tbl$isoform_id))
    ))

  # ── Per-gene pulldown loop ───────────────────────────────────────────────────
  gene_ids     <- unique(tbl$gene_id)
  gene_results <- stats::setNames(
    vector("list", length(gene_ids)),
    gene_ids
  )

  for (gid in gene_ids) {
    if (verbose) message(sprintf("  Processing gene: %s", gid))

    gene_tbl      <- tbl[tbl$gene_id == gid, , drop = FALSE]
    up_isoforms   <- gene_tbl$isoform_id[gene_tbl$direction == "up"]
    down_isoforms <- gene_tbl$isoform_id[gene_tbl$direction == "down"]

    up_graph   <- .safe_pulldown(obj, up_isoforms,   order, cutoff, gid, "up",   verbose)
    down_graph <- .safe_pulldown(obj, down_isoforms, order, cutoff, gid, "down", verbose)

    up_stats   <- .safe_analyze(up_graph)
    down_stats <- .safe_analyze(down_graph)

    gene_results[[gid]] <- .new_gene_diff_pulldown(
      gene_id       = gid,
      up_isoforms   = up_isoforms,
      down_isoforms = down_isoforms,
      up_graph      = up_graph,
      down_graph    = down_graph,
      up_stats      = up_stats,
      down_stats    = down_stats
    )
  }

  # ── Assemble top-level result ────────────────────────────────────────────────
  if (is.null(label))
    label <- paste0("diff_pulldown_", format(Sys.time(), "%Y%m%d_%H%M%S"))

  result <- .new_diff_pulldown_analysis(
    label        = label,
    input_table  = tbl,
    gene_results = gene_results,
    order        = order,
    cutoff       = cutoff
  )

  if (verbose)
    message(sprintf("Done. DifferentialPulldownAnalysis '%s' created.", label))

  if (!is.character(label) || length(label) != 1L || nchar(trimws(label)) == 0L)
    stop("`label` must be a single non-empty character string.", call. = FALSE)

  if (label %in% names(obj@analyses))
    warning(sprintf(
      "An analysis labelled '%s' already exists and will be overwritten.", label
    ))

  obj@analyses[[label]] <- result
  message(sprintf("Analysis '%s' stored in obj@analyses.", label))
  invisible(obj)
}


# ==============================================================================
#  Internal: input loading
# ==============================================================================

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

  if (ext == "csv") {
    return(utils::read.csv(input, stringsAsFactors = FALSE))
  }

  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop(
        "Package 'readxl' is required to read Excel files: ",
        "install.packages('readxl').",
        call. = FALSE
      )
    return(as.data.frame(readxl::read_excel(input)))
  }

  stop(
    sprintf(
      "Unsupported file extension '.%s'. Supply a CSV or Excel (.xlsx/.xls) file.",
      ext
    ),
    call. = FALSE
  )
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

  # Normalise strings
  tbl$gene_id    <- trimws(as.character(tbl$gene_id))
  tbl$isoform_id <- trimws(as.character(tbl$isoform_id))
  tbl$direction  <- tolower(trimws(as.character(tbl$direction)))

  # Drop rows with empty identifiers, with a warning
  blank_gene    <- nchar(tbl$gene_id)    == 0L
  blank_isoform <- nchar(tbl$isoform_id) == 0L

  if (any(blank_gene | blank_isoform)) {
    n_drop <- sum(blank_gene | blank_isoform)
    warning(sprintf(
      "Dropping %d row(s) with empty gene_id or isoform_id.", n_drop
    ), call. = FALSE)
    tbl <- tbl[!(blank_gene | blank_isoform), , drop = FALSE]
  }

  # Check direction values
  valid_dirs   <- c("up", "down")
  bad_dirs     <- !tbl$direction %in% valid_dirs
  if (any(bad_dirs)) {
    bad_vals <- unique(tbl$direction[bad_dirs])
    stop(sprintf(
      "Unrecognised direction value(s): %s. Must be 'up' or 'down'.",
      paste(bad_vals, collapse = ", ")
    ), call. = FALSE)
  }

  if (nrow(tbl) == 0L)
    stop("Input table is empty after cleaning.", call. = FALSE)

  tbl
}


# ==============================================================================
#  Internal: safe pulldown wrapper
# ==============================================================================

#' Run virtual_pulldown, returning NULL instead of stopping on failure
#'
#' Catches errors that arise when seed isoforms are absent from the network or
#' no nodes survive the cutoff, and emits a warning instead so the per-gene
#' loop continues.
#'
#' @param obj       IsoNetAnalyzeR object.
#' @param isoforms  Character vector of seed isoform IDs.
#' @param order     Integer.
#' @param cutoff    Numeric.
#' @param gene_id   Character scalar — used in the warning message only.
#' @param direction Character scalar — "up" or "down", for warning messages.
#' @param verbose   Logical.
#'
#' @return igraph subgraph, or NULL on failure.
#' @noRd
.safe_pulldown <- function(obj, isoforms, order, cutoff,
                           gene_id, direction, verbose) {

  if (length(isoforms) == 0L) {
    if (verbose)
      message(sprintf(
        "    [%s] No %s-regulated isoforms — skipping.", gene_id, direction
      ))
    return(NULL)
  }

  tryCatch(
    virtual_pulldown(
      obj          = obj,
      seed_vertices = isoforms,
      order        = order,
      cutoff       = cutoff
    ),
    error = function(e) {
      warning(sprintf(
        "Gene '%s' (%s): virtual_pulldown failed — %s",
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
  if (is.null(g)) {
    return(data.frame(
      isoform            = character(0),
      degree             = integer(0),
      clustering         = numeric(0),
      distance_from_seed = numeric(0),
      stringsAsFactors   = FALSE
    ))
  }
  analyze_subgraph(g)
}
