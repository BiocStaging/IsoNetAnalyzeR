## IsoNetAnalyzeR-class.R
## S4 class definition for the main IsoNetAnalyzeR session object.
##
## Exported: IsoNetAnalyzeR (constructor), show method,
##           get_network_diagnostics

setOldClass("igraph")
#' @importFrom methods is new setClass setGeneric setMethod setValidity validObject
#' @importFrom rlang .data
NULL

# ==============================================================================
#  Class definition
# ==============================================================================

#' IsoNetAnalyzeR - main session object
#'
#' The central object of the \pkg{IsoNetAnalyzeR} package. It holds the
#' isoform interaction network, session metadata, a named list of full
#' analyses (e.g. virtual pulldown runs), and a named list of quick look-ups.
#'
#' @slot name        Character scalar. A user-defined label for the session.
#' @slot network     An \code{igraph} object representing the isoform
#'   interaction network. Populated from the package's internal data by
#'   default; can be replaced with a user-supplied graph in the future.
#' @slot created_at  \code{POSIXct} scalar. Timestamp recorded when the object
#'   was constructed.
#' @slot metadata    Named list of arbitrary session-level metadata (e.g.
#'   \code{pi_threshold}, \code{network_cutoff} carried over from internal
#'   data, or user-supplied fields).
#' @slot analyses    Named list of \linkS4class{DifferentialPulldownAnalysis}
#'   objects (and potentially other analysis classes added in the future).
#' @slot lookups     Named list of \linkS4class{VirtualPulldownLookup} objects
#'   (and potentially other look-up classes added in the future).
#' @slot network_diagnostics Either \code{NULL} (before diagnostics have been
#'   run) or a \linkS4class{NetworkDiagnostics} object produced by
#'   \code{\link{run_network_diagnostics}}.
#'
#' @exportClass IsoNetAnalyzeR
setClass(
  "IsoNetAnalyzeR",
  representation(
    name                = "character",
    network             = "ANY",        # igraph; ANY avoids a hard igraph S4 import
    created_at          = "POSIXct",
    metadata            = "list",
    analyses            = "list",
    lookups             = "list",
    network_diagnostics = "ANY"         # NULL or NetworkDiagnostics
  )
)


# ==============================================================================
#  Validator
# ==============================================================================

setValidity("IsoNetAnalyzeR", function(object) {
  msgs <- character()

  if (length(object@name) != 1L || nchar(trimws(object@name)) == 0L)
    msgs <- c(msgs, "@name must be a non-empty character scalar.")

  if (!inherits(object@network, "igraph"))
    msgs <- c(msgs, "@network must be an igraph object.")

  if (length(object@created_at) != 1L)
    msgs <- c(msgs, "@created_at must be a length-1 POSIXct.")

  # Every element of @analyses must carry the right S4 class
  allowed <- c("DifferentialPulldownAnalysis")

  bad_analyses <- vapply(
    object@analyses,
    function(x) !any(vapply(allowed, function(cls) is(x, cls), logical(1))),
    logical(1)
  )

  if (any(bad_analyses))
    msgs <- c(msgs, paste0(
      "@analyses contains unsupported object(s): ",
      paste(names(object@analyses)[bad_analyses], collapse = ", ")
    ))

  # Every element of @lookups must carry the right S4 class
  bad_lookups <- vapply(
    object@lookups,
    function(x) !is(x, "VirtualPulldownLookup"),
    logical(1)
  )
  if (any(bad_lookups))
    msgs <- c(msgs, paste0(
      "@lookups contains unsupported object(s): ",
      paste(names(object@lookups)[bad_lookups], collapse = ", ")
    ))

  # @network_diagnostics must be NULL or a NetworkDiagnostics object
  if (!is.null(object@network_diagnostics) &&
      !is(object@network_diagnostics, "NetworkDiagnostics"))
    msgs <- c(msgs,
              "@network_diagnostics must be NULL or a NetworkDiagnostics object.")

  if (length(msgs)) msgs else TRUE
})


# ==============================================================================
#  Constructor
# ==============================================================================

#' Create a new IsoNetAnalyzeR session object
#'
#' Initialises an \linkS4class{IsoNetAnalyzeR} object. When no \code{network}
#' is supplied the package's internal isoform interaction network
#' (\code{.network_data}) is loaded automatically, along with the associated
#' \code{.network_cutoff} and \code{.pi_threshold} values stored as metadata.
#'
#' @param name     Character scalar. A human-readable label for the session
#'   (e.g. \code{"my_experiment"}).  Defaults to \code{"IsoNetAnalyzeR"}.
#' @param network  An \code{igraph} object to use as the interaction network.
#'   Pass \code{NULL} (the default) to load the package's built-in network.
#' @param metadata Named list of additional metadata to attach to the session.
#'   This is merged with any automatically populated fields
#'   (\code{pi_threshold}, \code{network_cutoff}) - user-supplied values take
#'   precedence on name collision.
#'
#' @return An \linkS4class{IsoNetAnalyzeR} object.
#'
#' @examples
#' \dontrun{
#' # Load the built-in network
#' obj <- IsoNetAnalyzeR("my_experiment")
#'
#' # Supply a custom network
#' obj <- IsoNetAnalyzeR("custom", network = my_igraph)
#' }
#'
#' @export
IsoNetAnalyzeR <- function(name     = "IsoNetAnalyzeR",
                           network  = NULL,
                           metadata = list()) {

  # -- Network resolution ------------------------------------------------------
  auto_meta <- list()

  if (is.null(network)) {
    if (!exists(".network_data", envir = asNamespace(packageName()))) {
      stop(
        "No network supplied and the package's internal network (.network_data) ",
        "could not be found. Have you built the package with `devtools::load_all()`?",
        call. = FALSE
      )
    }
    network   <- get(".network_data",  envir = asNamespace(packageName()))
    auto_meta <- list(
      network_source  = "package_internal",
      pi_threshold    = tryCatch(
        get(".pi_threshold",   envir = asNamespace(packageName())),
        error = function(e) NA_real_
      ),
      network_cutoff  = tryCatch(
        get(".network_cutoff", envir = asNamespace(packageName())),
        error = function(e) NA_real_
      )
    )
    message(sprintf(
      "Loaded package-internal network: %d vertices, %d edges.",
      igraph::vcount(network), igraph::ecount(network)
    ))
  } else {
    if (!inherits(network, "igraph"))
      stop("`network` must be an igraph object.", call. = FALSE)
    auto_meta <- list(network_source = "user_supplied")
    message(sprintf(
      "Using user-supplied network: %d vertices, %d edges.",
      igraph::vcount(network), igraph::ecount(network)
    ))
  }

  # User metadata wins on collision
  merged_meta <- utils::modifyList(auto_meta, metadata)

  #-- Construct & validate -----------------------------------------------------
  obj <- new(
    "IsoNetAnalyzeR",
    name       = as.character(name),
    network    = network,
    created_at = Sys.time(),
    metadata            = merged_meta,
    analyses            = list(),
    lookups             = list(),
    network_diagnostics = NULL
  )

  validObject(obj)
  obj
}


# ==============================================================================
#  show()
# ==============================================================================

#' @describeIn IsoNetAnalyzeR Print a concise summary of an IsoNetAnalyzeR
#'   session object.
#' @param object An \linkS4class{IsoNetAnalyzeR} object.
setMethod("show", "IsoNetAnalyzeR", function(object) {
  cli::cli_rule("IsoNetAnalyzeR session")
  cli::cli_text("Name       : {object@name}")
  cli::cli_text("Created    : {format(object@created_at, '%Y-%m-%d %H:%M:%S')}")
  cli::cli_text("Network    : {igraph::vcount(object@network)} vertices / {igraph::ecount(object@network)} edges  [source: {object@metadata[['network_source']] %||% 'unknown'}]")
  cli::cli_text("Analyses   : {length(object@analyses)} stored  [{(.slot_names_str(object@analyses))}]")
  cli::cli_text("Lookups    : {length(object@lookups)} stored  [{(.slot_names_str(object@lookups))}]")
  diag_status <- if (is.null(object@network_diagnostics)) "not run" else
    format(object@network_diagnostics@created_at, '%Y-%m-%d %H:%M:%S')
  cli::cli_text("Diagnostics: {diag_status}")

  if (length(object@metadata) > 0L) {
    extra <- setdiff(names(object@metadata), "network_source")
    if (length(extra)) {
      cli::cli_text("Metadata   : {paste(extra, collapse = ', ')}")
    }
  }
  cli::cli_rule()
  invisible(object)
})


# ==============================================================================
#  Internal helpers
# ==============================================================================

#' Null-coalescing operator (internal)
#' @noRd
`%||%` <- function(x, y) if (!is.null(x)) x else y

#' Format a named list's keys as a comma-separated string
#' @noRd
.slot_names_str <- function(lst) {
  if (length(lst) == 0L) return("none")
  paste(names(lst), collapse = ", ")
}

#' Return the package name from within the package namespace
#' @noRd
packageName <- function() "IsoNetAnalyzeR"


# ==============================================================================
#  Accessor generics and methods
# ==============================================================================

#' List stored analysis labels
#'
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
#' @return Character vector of analysis labels, or \code{character(0)} if none.
#' @export
setGeneric("list_analyses", function(obj) standardGeneric("list_analyses"))

#' @describeIn IsoNetAnalyzeR List labels of stored analyses.
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
setMethod("list_analyses", "IsoNetAnalyzeR", function(obj) {
  names(obj@analyses)
})


#' Retrieve a stored analysis result
#'
#' @param obj   An \linkS4class{IsoNetAnalyzeR} session object.
#' @param label Character scalar. The key under which the analysis was stored.
#'
#' @return A \linkS4class{DifferentialPulldownAnalysis} object.
#' @export
setGeneric("get_analysis", function(obj, label) standardGeneric("get_analysis"))

#' @describeIn IsoNetAnalyzeR Retrieve a stored analysis by label.
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
#' @param label Character scalar. The key under which the analysis was stored.
setMethod("get_analysis", "IsoNetAnalyzeR", function(obj, label) {
  if (!label %in% names(obj@analyses))
    stop(sprintf(
      "No analysis labelled '%s'. Available: %s",
      label, .slot_names_str(obj@analyses)
    ), call. = FALSE)
  obj@analyses[[label]]
})


#' List stored lookup labels
#'
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
#' @return Character vector of lookup labels, or \code{character(0)} if none.
#' @export
setGeneric("list_lookups", function(obj) standardGeneric("list_lookups"))

#' @describeIn IsoNetAnalyzeR List labels of stored lookups.
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
setMethod("list_lookups", "IsoNetAnalyzeR", function(obj) {
  names(obj@lookups)
})


#' Retrieve a stored lookup result
#'
#' @param obj   An \linkS4class{IsoNetAnalyzeR} session object.
#' @param label Character scalar. The key under which the lookup was stored.
#'
#' @return A \linkS4class{VirtualPulldownLookup} object.
#' @export
setGeneric("get_lookup", function(obj, label) standardGeneric("get_lookup"))

#' @describeIn IsoNetAnalyzeR Retrieve a stored lookup by label.
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
#' @param label Character scalar. The key under which the lookup was stored.
setMethod("get_lookup", "IsoNetAnalyzeR", function(obj, label) {
  if (!label %in% names(obj@lookups))
    stop(sprintf(
      "No lookup labelled '%s'. Available: %s",
      label, .slot_names_str(obj@lookups)
    ), call. = FALSE)
  obj@lookups[[label]]
})


# ==============================================================================
#  get_network_diagnostics
# ==============================================================================

#' Retrieve the network diagnostics result from an IsoNetAnalyzeR object
#'
#' Returns the \linkS4class{NetworkDiagnostics} object stored in the session,
#' or stops with an informative message if \code{\link{run_network_diagnostics}}
#' has not been called yet.
#'
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
#'
#' @return A \linkS4class{NetworkDiagnostics} object.
#'
#' @examples
#' \dontrun{
#' obj <- run_network_diagnostics(obj)
#' get_network_diagnostics(obj)
#' }
#'
#' @export
setGeneric("get_network_diagnostics",
           function(obj) standardGeneric("get_network_diagnostics"))

#' @describeIn get_network_diagnostics Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_network_diagnostics", "IsoNetAnalyzeR", function(obj) {
  if (is.null(obj@network_diagnostics))
    stop(
      "No network diagnostics found. Run `run_network_diagnostics(obj)` first.",
      call. = FALSE
    )
  obj@network_diagnostics
})
