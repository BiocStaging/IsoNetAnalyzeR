# ==============================================================================
#  lookup_virtual_pulldown — single quick look-up (supports multiple seeds)
# ==============================================================================

#' Perform a virtual pulldown look-up with one or more seeds
#'
#' Runs a virtual pulldown for one or more seed isoforms or genes and returns
#' the result as a \linkS4class{VirtualPulldownLookup} object.  All supplied
#' seeds are treated as order-0 anchors in a single pulldown.  The network is
#' taken directly from \code{obj}; the cutoff defaults to the value stored in
#' \code{obj@@metadata[["network_cutoff"]]} but can be overridden to explore
#' alternative thresholds.  The result is optionally stored inside the session
#' object (under \code{obj@@lookups}) when \code{store = TRUE}.
#'
#' @param obj    An \linkS4class{IsoNetAnalyzeR} session object.
#' @param seed   Character vector. One or more isoform vertex names or gene
#'   names.  When multiple seeds are given they are all used as order-0 anchors
#'   in the same virtual pulldown.
#' @param order  Integer. Number of interaction hops.
#' @param cutoff Numeric in \[0, 1\]. Internal-edge fraction threshold. Falls
#'   back to \code{obj@@metadata[["network_cutoff"]]} when \code{NULL}.
#'   Supply an explicit value to explore alternative thresholds without
#'   modifying the session object.
#' @param level  Character scalar. \code{"isoform"} (default) or \code{"gene"}.
#' @param store  Logical. If \code{TRUE} the result is stored in
#'   \code{obj@@lookups} and the updated \code{obj} is returned invisibly.  If
#'   \code{FALSE} (default) only the \linkS4class{VirtualPulldownLookup} object
#'   is returned.
#' @param label  Character scalar. Key under which to store the look-up when
#'   \code{store = TRUE}.  Defaults to \code{"<seed1>_AND_<seed2>_<timestamp>"}.
#'
#' @return When \code{store = FALSE}: a \linkS4class{VirtualPulldownLookup}
#'   object.  When \code{store = TRUE}: the updated
#'   \linkS4class{IsoNetAnalyzeR} object (invisibly); access the look-up via
#'   \code{obj@@lookups[[label]]}.
#'
#' @examples
#' \dontrun{
#' obj <- IsoNetAnalyzeR("demo")
#'
#' # Single seed — uses stored network_cutoff
#' lkp <- lookup_virtual_pulldown(obj, seed = "ENST00000379236", order = 1)
#' lkp
#'
#' # Multiple seeds — all treated as order-0 anchors in one pulldown
#' lkp <- lookup_virtual_pulldown(
#'   obj,
#'   seed  = c("ENST00000379236", "ENST00000445118"),
#'   order = 1
#' )
#'
#' # Custom cutoff — explore a stricter threshold without changing the session
#' lkp <- lookup_virtual_pulldown(
#'   obj,
#'   seed   = "ENST00000379236",
#'   order  = 1,
#'   cutoff = 0.4
#' )
#'
#' # Gene-level seeds — all isoforms of each gene become order-0 anchors
#' lkp <- lookup_virtual_pulldown(
#'   obj,
#'   seed  = c("BRCA1", "TP53"),
#'   order = 1,
#'   level = "gene"
#' )
#'
#' # Persistent look-up — result stored in obj
#' obj <- lookup_virtual_pulldown(
#'   obj,
#'   seed  = c("ENST00000379236", "ENST00000445118"),
#'   order = 1,
#'   store = TRUE,
#'   label = "combined_lkp"
#' )
#' obj@@lookups$combined_lkp
#' }
#'
#' @export
lookup_virtual_pulldown <- function(obj,
                                    seed,
                                    order,
                                    cutoff = NULL,
                                    level  = c("isoform", "gene"),
                                    store  = FALSE,
                                    label  = NULL) {

  level <- match.arg(level)

  # ── Input validation ──────────────────────────────────────────────────────────
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.character(seed) || length(seed) == 0L || any(nchar(trimws(seed)) == 0L))
    stop("`seed` must be a non-empty character vector with no blank strings.",
         call. = FALSE)

  if (!is.numeric(order) || length(order) != 1L || order < 1L)
    stop("`order` must be a positive integer scalar.", call. = FALSE)

  order <- as.integer(order)

  # ── Cutoff resolution ─────────────────────────────────────────────────────────
  # Resolved here so we can report the value in the progress message and store
  # it on the lookup object — virtual_pulldown() will resolve again internally
  # but that is a cheap operation.
  if (is.null(cutoff)) {
    cutoff <- obj@metadata[["network_cutoff"]]
    if (is.null(cutoff) || is.na(cutoff))
      stop(
        "`cutoff` not supplied and `obj@metadata$network_cutoff` is not set. ",
        "Please provide a numeric cutoff explicitly.",
        call. = FALSE
      )
    message(sprintf("Using stored network_cutoff: %.4f", cutoff))
  } else {
    if (!is.numeric(cutoff) || length(cutoff) != 1L ||
        cutoff < 0 || cutoff > 1)
      stop("`cutoff` must be a single numeric value in [0, 1].", call. = FALSE)
    message(sprintf("Using user-supplied cutoff: %.4f", cutoff))
  }

  # ── Seed resolution ───────────────────────────────────────────────────────────
  resolved <- .resolve_seeds(obj@network, seed, level)

  missing_seeds <- names(resolved)[vapply(resolved, length, integer(1L)) == 0L]
  if (length(missing_seeds) > 0L)
    stop(sprintf(
      "The following seed(s) were not found in the network: %s",
      paste(missing_seeds, collapse = ", ")
    ), call. = FALSE)

  # Flatten all resolved vertices into one combined seed vector —
  # all seeds are treated as order-0 anchors in a single pulldown
  sv <- unique(unlist(resolved, use.names = FALSE))

  # ── Run pulldown ──────────────────────────────────────────────────────────────
  message(sprintf(
    "Running virtual pulldown: %d seed vertex/vertices, order = %d, cutoff = %.4f",
    length(sv), order, cutoff
  ))

  subgraph   <- virtual_pulldown(obj = obj, seed_vertices = sv, order = order,
                                 cutoff = cutoff)
  node_stats <- analyze_subgraph(subgraph)

  lookup_obj <- .new_vp_lookup(
    seed       = seed,
    level      = level,
    subgraph   = subgraph,
    node_stats = node_stats,
    cutoff     = cutoff,
    params     = list(order = order, level = level)
  )

  # ── Return or store ───────────────────────────────────────────────────────────
  if (!store)
    return(lookup_obj)

  # Build a safe default label from all seeds
  if (is.null(label))
    label <- paste0(
      paste(gsub("[^A-Za-z0-9_]", "_", seed), collapse = "_AND_"), "_",
      format(Sys.time(), "%Y%m%d_%H%M%S")
    )

  if (!is.character(label) || length(label) != 1L || nchar(trimws(label)) == 0L)
    stop("`label` must be a single non-empty character string.", call. = FALSE)

  if (label %in% names(obj@lookups))
    warning(sprintf(
      "A look-up labelled '%s' already exists and will be overwritten.", label
    ))

  obj@lookups[[label]] <- lookup_obj
  message(sprintf("Look-up '%s' stored in obj@lookups.", label))
  invisible(obj)
}


# ==============================================================================
#  Internal: seed resolution helper
# ==============================================================================

#' Resolve seeds to vertex names in the network
#'
#' At isoform level each seed is returned as-is (wrapped in a list).
#' At gene level the network's \code{gene} vertex attribute is used to
#' collect all isoforms belonging to each gene.
#'
#' @param network igraph object.
#' @param seeds   Character vector of seeds.
#' @param level   \code{"isoform"} or \code{"gene"}.
#'
#' @return Named list (same names as \code{seeds}), each element a character
#'   vector of isoform vertex names.
#' @noRd
.resolve_seeds <- function(network, seeds, level) {
  if (level == "isoform") {
    return(stats::setNames(as.list(seeds), seeds))
  }

  # Gene level: look up via vertex attribute
  if (is.null(igraph::V(network)$gene))
    stop(
      "The network does not carry a 'gene' vertex attribute, ",
      "which is required for level = 'gene'.",
      call. = FALSE
    )

  gene_attr <- igraph::V(network)$gene
  v_names   <- igraph::V(network)$name

  stats::setNames(
    lapply(seeds, function(g) {
      v_names[!is.na(gene_attr) & gene_attr == g]
    }),
    seeds
  )
}
