## VirtualPulldownLookup-class.R
## S4 class definition for a single virtual pulldown look-up.
##
## Exported: VirtualPulldownLookup (constructor), show method


# ==============================================================================
#  Class definition
# ==============================================================================

#' VirtualPulldownLookup — single quick virtual pulldown look-up
#'
#' Stores the result of a one-shot virtual pulldown query at either isoform or
#' gene level, without running a full analysis across many seeds.  Can
#' optionally be stored inside an \linkS4class{IsoNetAnalyzeR} session object.
#'
#' @slot seed       Character vector. The seed isoform(s) or gene(s) queried.
#'   Multiple seeds are stored as a vector when all were used as order-0
#'   anchors in a single pulldown.
#' @slot level      Character scalar. Either \code{"isoform"} or \code{"gene"}.
#' @slot subgraph   An \code{igraph} object representing the pulldown subnetwork.
#' @slot node_stats Data frame. Per-node statistics for the subgraph
#'   (output of \code{\link{analyze_subgraph}}).
#' @slot cutoff     Numeric scalar. The internal-edge fraction cutoff used
#'   for this look-up. Stored explicitly so it is always visible on the object
#'   independent of the session metadata.
#' @slot params     Named list. Additional parameters used: \code{order},
#'   \code{level}.
#' @slot created_at \code{POSIXct} scalar. Timestamp of the look-up.
#'
#' @exportClass VirtualPulldownLookup
setClass(
  "VirtualPulldownLookup",
  representation(
    seed       = "character",
    level      = "character",
    subgraph   = "ANY",        # igraph
    node_stats = "data.frame",
    cutoff     = "numeric",
    params     = "list",
    created_at = "POSIXct"
  )
)


# ==============================================================================
#  Validator
# ==============================================================================

setValidity("VirtualPulldownLookup", function(object) {
  msgs <- character()

  if (!is.character(object@seed) ||
      length(object@seed) == 0L ||
      any(nchar(trimws(object@seed)) == 0L))
    msgs <- c(msgs, "@seed must be a non-empty character vector with no blank strings.")

  if (!object@level %in% c("isoform", "gene"))
    msgs <- c(msgs, "@level must be 'isoform' or 'gene'.")

  if (!inherits(object@subgraph, "igraph"))
    msgs <- c(msgs, "@subgraph must be an igraph object.")

  if (length(object@cutoff) != 1L || is.na(object@cutoff) ||
      object@cutoff < 0 || object@cutoff > 1)
    msgs <- c(msgs, "@cutoff must be a single numeric value in [0, 1].")

  if (length(object@created_at) != 1L)
    msgs <- c(msgs, "@created_at must be a length-1 POSIXct.")

  if (length(msgs)) msgs else TRUE
})


# ==============================================================================
#  Constructor (internal — called by lookup_virtual_pulldown)
# ==============================================================================

#' Build a VirtualPulldownLookup object
#'
#' Internal constructor called by \code{\link{lookup_virtual_pulldown}}.
#'
#' @param seed       Character vector. The queried seed(s).
#' @param level      Character scalar. \code{"isoform"} or \code{"gene"}.
#' @param subgraph   igraph subgraph.
#' @param node_stats Data frame from \code{\link{analyze_subgraph}}.
#' @param cutoff     Numeric scalar. The cutoff used for this look-up.
#' @param params     Named list of additional run parameters.
#'
#' @return A \linkS4class{VirtualPulldownLookup} object.
#' @noRd
.new_vp_lookup <- function(seed, level, subgraph, node_stats, cutoff, params) {
  obj <- new(
    "VirtualPulldownLookup",
    seed       = as.character(seed),
    level      = as.character(level),
    subgraph   = subgraph,
    node_stats = node_stats,
    cutoff     = as.numeric(cutoff),
    params     = params,
    created_at = Sys.time()
  )
  validObject(obj)
  obj
}


# ==============================================================================
#  show()
# ==============================================================================

#' @describeIn VirtualPulldownLookup Print a concise summary.
setMethod("show", "VirtualPulldownLookup", function(object) {

  # Format seed display: single seed inline, multiple seeds on separate lines
  seed_str <- if (length(object@seed) == 1L) {
    object@seed
  } else {
    paste0(
      "\n    ", paste(object@seed, collapse = "\n    ")
    )
  }

  cat(
    "── VirtualPulldownLookup ────────────────────────────────────────\n",
    sprintf("  Seed(s)    : %s  [%s level]\n", seed_str, object@level),
    sprintf("  Created    : %s\n", format(object@created_at, "%Y-%m-%d %H:%M:%S")),
    sprintf("  Order      : %s\n", object@params[["order"]]  %||% "?"),
    sprintf("  Cutoff     : %.4f\n", object@cutoff),
    sprintf("  Subgraph   : %d vertices, %d edges\n",
            igraph::vcount(object@subgraph),
            igraph::ecount(object@subgraph)),
    sprintf("  Node stats : %d row(s), %d col(s)\n",
            nrow(object@node_stats), ncol(object@node_stats)),
    sep = ""
  )
  cat("────────────────────────────────────────────────────────────────\n")
  invisible(object)
})
