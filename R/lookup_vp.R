#' @include isonetanalyzer_class.R
NULL
#' @include lookup_class.R
NULL

# ==============================================================================
#  lookup_virtual_pulldown - single quick look-up (supports multiple seeds)
# ==============================================================================

#' Perform a virtual pulldown look-up with one or more seeds
#'
#' Runs a virtual pulldown for one or more seed isoforms or genes and adds
#' the result as a \linkS4class{VirtualPulldownLookup} object to the
#' \linkS4class{IsoNetAnalyzeR} object.  All supplied seeds are treated as
#' order-0 anchors in a single pulldown.  The network is taken directly from
#' \code{obj}; the cutoff defaults to the value stored in
#' \code{obj@@metadata[["network_cutoff"]]} but can be overridden to explore
#' alternative thresholds.
#'
#' @param obj    An \linkS4class{IsoNetAnalyzeR} session object.
#' @param label  Character scalar. Key under which to store the look-up.
#' @param seed   Character vector. One or more isoform vertex names or gene
#'   names.  When multiple seeds are given they are all used as order-0 anchors
#'   in the same virtual pulldown.
#' @param order  Integer. Number of interaction hops.
#' @param cutoff Numeric in \[0, 1\]. Internal-edge fraction threshold. Falls
#'   back to \code{obj@@metadata[["network_cutoff"]]} when \code{NULL}.
#'   Supply an explicit value to explore alternative thresholds without
#'   modifying the session object.
#' @param level  Character scalar. \code{"isoform"} (default) or \code{"gene"}.
#' @param resolve_duplicates Logical. For isoform-level lookups, resolve
#'   versioned and protein-identical redundant transcript IDs against the
#'   active network. Default \code{TRUE}.
#'
#' @return An updated \linkS4class{IsoNetAnalyzeR} object (invisibly);
#' access the look-up via \code{get_lookup(obj, label)}.
#'
#' @examples
#' obj <- IsoNetAnalyzeR("example")
#'
#' # One isoform seed, using the session's stored cutoff
#' obj <- lookup_virtual_pulldown(obj, "cltcl1_up", seed = "ENST00000621271",
#'                                order = 1)
#'
#' # Several seeds act as joint order-0 anchors; try a stricter cutoff
#' obj <- lookup_virtual_pulldown(
#'   obj, "cltcl1_both",
#'   seed = c("ENST00000621271", "ENST00000427926"),
#'   order = 1, cutoff = 0.4
#' )
#'
#' # All network isoforms of a gene as seeds
#' obj <- lookup_virtual_pulldown(obj, "cltcl1_gene", seed = "ENSG00000070371",
#'                                order = 1, level = "gene")
#' list_lookups(obj)
#' get_lookup(obj, "cltcl1_gene")
#'
#' @export
lookup_virtual_pulldown <- function(obj,
                                    label,
                                    seed,
                                    order,
                                    cutoff = NULL,
                                    level  = c("isoform", "gene"),
                                    resolve_duplicates = TRUE) {

  obj_name <- .obj_name(substitute(obj))

  level <- match.arg(level)

  # -- Input validation --------------------------------------------------------
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(trimws(label)))
    stop("`label` must be a single non-empty character string.", call. = FALSE)

  if (label %in% names(obj@lookups))
    warning(sprintf(
      "A look-up labelled '%s' already exists and will be overwritten.", label
    ))

  if (!is.character(seed) || length(seed) == 0L || anyNA(seed) ||
      any(!nzchar(trimws(seed))))
    stop("`seed` must be a non-empty character vector with no blank strings.",
         call. = FALSE)

  if (!is.logical(resolve_duplicates) || length(resolve_duplicates) != 1L ||
      is.na(resolve_duplicates))
    stop("`resolve_duplicates` must be a single non-NA logical value.",
         call. = FALSE)

  order <- .validate_positive_integer(order, "order")

  # -- Cutoff resolution -------------------------------------------------------
  # Resolved here so we can report the value in the progress message and store
  # it on the lookup object — virtual_pulldown() will resolve again internally
  # but that is a cheap operation.
  if (is.null(cutoff)) {
    cutoff <- obj@metadata[["network_cutoff"]]
    if (is.null(cutoff))
      stop(
        "`cutoff` not supplied and `obj@metadata$network_cutoff` is not set. ",
        "Please provide a numeric cutoff explicitly.",
        call. = FALSE
      )
    cutoff <- .validate_cutoff(cutoff, "obj@metadata$network_cutoff")
    message(sprintf("Using stored network_cutoff: %.4f", cutoff))
  } else {
    cutoff <- .validate_cutoff(cutoff)
    message(sprintf("Using user-supplied cutoff: %.4f", cutoff))
  }

  # -- Seed resolution ---------------------------------------------------------
  resolution_audit <- NULL
  if (identical(level, "isoform") && resolve_duplicates) {
    seed <- resolve_duplicate_isoforms(seed, network = obj@network)
    resolution_audit <- attr(seed, "dedup_audit")
    .warn_duplicate_resolution(resolution_audit)
    seed <- unique(as.character(seed))
  }
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

  # -- Run pulldown ------------------------------------------------------------
  message(sprintf(
    "Running virtual pulldown: %d seed vertex/vertices, order = %d, cutoff = %.4f",
    length(sv), order, cutoff
  ))

  subgraph   <- virtual_pulldown(obj = obj, seed_vertices = sv, order = order,
                                 cutoff = cutoff, resolve_duplicates = FALSE)
  node_stats <- analyze_subgraph(subgraph)

  lookup_obj <- .new_vp_lookup(
    seed       = seed,
    level      = level,
    subgraph   = subgraph,
    node_stats = node_stats,
    cutoff     = cutoff,
    params     = list(
      order = order, level = level,
      duplicate_resolution = resolution_audit
    )
  )

  obj@lookups[[label]] <- lookup_obj
  message(
    sprintf("Look-up '%s' stored. Access with: get_lookup(%s, \"%s\")",
            label, obj_name, label))
  invisible(obj)
}


# ==============================================================================
#  Internal: seed resolution helper
# ==============================================================================

#' Resolve seeds to vertex names in the network
#'
#' At isoform level each seed is checked against network membership and returned
#' as a one-element character vector, or \code{character(0)} when absent.
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
    network_vertices <- igraph::V(network)$name
    return(stats::setNames(
      lapply(seeds, function(seed) {
        if (seed %in% network_vertices) seed else character(0)
      }),
      seeds
    ))
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
