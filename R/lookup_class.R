## VirtualPulldownLookup-class.R
## S4 class definition for a single virtual pulldown look-up.
##
## Exported: VirtualPulldownLookup (constructor), show method,
##           get_lookup_stats, get_lookup_subgraph


# ==============================================================================
#  Class definition
# ==============================================================================

#' VirtualPulldownLookup - single quick virtual pulldown look-up
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
#  Constructor (internal - called by lookup_virtual_pulldown)
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
#' @param object An \linkS4class{VirtualPulldownLookup} session object.
setMethod("show", "VirtualPulldownLookup", function(object) {

  # Format seed display: single seed inline, multiple seeds on separate lines
  cli::cli_rule("VirtualPulldownLookup")

  if (length(object@seed) == 1L) {
    cli::cli_text("Seed(s)    : {object@seed}  [{object@level} level]")
  } else {
    cli::cli_text("Seed(s)    : [{object@level} level]")
    for (s in object@seed) {
      cli::cli_text("    {s}")
    }
  }

  cli::cli_text("Created    : {format(object@created_at, '%Y-%m-%d %H:%M:%S')}")
  cli::cli_text("Order      : {object@params[['order']] %||% '?'}")
  cli::cli_text("Cutoff     : {round(object@cutoff, 4)}")
  cli::cli_text("Subgraph   : {igraph::vcount(object@subgraph)} vertices, {igraph::ecount(object@subgraph)} edges")
  cli::cli_text("Node stats : {nrow(object@node_stats)} row(s), {ncol(object@node_stats)} col(s)")
  cli::cli_rule()
  invisible(object)
})


# ==============================================================================
#  get_lookup_stats
# ==============================================================================

#' Get the per-node statistics from a virtual pulldown look-up
#'
#' Returns the \code{node_stats} data frame stored in a
#' \linkS4class{VirtualPulldownLookup} object directly, without having to
#' unpack the S4 slot manually.  The data frame contains one row per node in
#' the pulldown subgraph and is the output of \code{\link{analyze_subgraph}}.
#'
#' @param obj   An \linkS4class{IsoNetAnalyzeR} session object.
#' @param label Character scalar.  Label of the
#'   \linkS4class{VirtualPulldownLookup}; use \code{list_lookups(obj)} to see
#'   available labels.
#'
#' @return A data frame with one row per node in the subgraph and columns
#'   produced by \code{\link{analyze_subgraph}}: \code{isoform},
#'   \code{is_seed}, \code{degree}, \code{clustering},
#'   \code{distance_from_seed}, \code{mean_pi}, \code{pi_entropy}.
#'
#' @examples
#' \dontrun{
#' # Run a look-up first
#' obj <- lookup_virtual_pulldown(obj, "my_lkp",
#'                                seed = "ENST00000379236", order = 1)
#'
#' # Retrieve the node stats data frame
#' stats <- get_lookup_stats(obj, "my_lkp")
#' head(stats)
#'
#' # Seed nodes only
#' stats[stats$is_seed, ]
#' }
#'
#' @export
setGeneric("get_lookup_stats",
           function(obj, label) standardGeneric("get_lookup_stats"))

#' @describeIn get_lookup_stats Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_lookup_stats", "IsoNetAnalyzeR",
          function(obj, label) {
            get_lookup(obj, label)@node_stats
          })


# ==============================================================================
#  get_lookup_subgraph
# ==============================================================================

#' Get the subgraph from a virtual pulldown look-up
#'
#' Returns the \code{igraph} subgraph stored in a
#' \linkS4class{VirtualPulldownLookup} object directly, without having to
#' unpack the S4 slot manually.  The subgraph can be passed directly to
#' \code{plot_network()}, \code{igraph} functions, or any other graph-level
#' utility.
#'
#' @param obj   An \linkS4class{IsoNetAnalyzeR} session object.
#' @param label Character scalar.  Label of the
#'   \linkS4class{VirtualPulldownLookup}; use \code{list_lookups(obj)} to see
#'   available labels.
#'
#' @return An \code{igraph} object representing the pulldown subnetwork.
#'
#' @examples
#' \dontrun{
#' # Run a look-up first
#' obj <- lookup_virtual_pulldown(obj, "my_lkp",
#'                                seed = "ENST00000379236", order = 1)
#'
#' # Retrieve the igraph subgraph
#' g <- get_lookup_subgraph(obj, "my_lkp")
#' igraph::vcount(g)
#' igraph::ecount(g)
#'
#' # Pass directly to plot_network
#' plot_network(g)
#' }
#'
#' @export
setGeneric("get_lookup_subgraph",
           function(obj, label) standardGeneric("get_lookup_subgraph"))

#' @describeIn get_lookup_subgraph Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_lookup_subgraph", "IsoNetAnalyzeR",
          function(obj, label) {
            get_lookup(obj, label)@subgraph
          })
