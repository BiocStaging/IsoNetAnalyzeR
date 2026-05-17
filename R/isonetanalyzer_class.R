## IsoNetAnalyzeR-class.R
## S4 class definition for the main IsoNetAnalyzeR session object.
##
## Exported: IsoNetAnalyzeR (constructor), show method,
##           get_session_name, set_session_name,
##           list_metadata, get_metadata, set_metadata

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
#'
#' @exportClass IsoNetAnalyzeR
setClass(
  "IsoNetAnalyzeR",
  representation(
    name       = "character",
    network    = "ANY",        # igraph; ANY avoids a hard igraph S4 import
    created_at = "POSIXct",
    metadata   = "list",
    analyses   = "list",
    lookups    = "list"
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
    metadata   = merged_meta,
    analyses   = list(),
    lookups    = list()
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
#  Session name accessors
# ==============================================================================

#' Get the session name of an IsoNetAnalyzeR object
#'
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
#' @return Character scalar. The session name stored in \code{@@name}.
#'
#' @examples
#' \dontrun{
#' get_session_name(obj)
#' }
#'
#' @export
setGeneric("get_session_name", function(obj) standardGeneric("get_session_name"))

#' @describeIn get_session_name Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_session_name", "IsoNetAnalyzeR", function(obj) {
  obj@name
})


#' Set the session name of an IsoNetAnalyzeR object
#'
#' @param obj   An \linkS4class{IsoNetAnalyzeR} session object.
#' @param value Character scalar. The new session name.
#'
#' @return The updated \linkS4class{IsoNetAnalyzeR} object, invisibly.
#'
#' @examples
#' \dontrun{
#' obj <- set_session_name(obj, "final_experiment")
#' get_session_name(obj)
#' }
#'
#' @export
setGeneric("set_session_name",
           function(obj, value) standardGeneric("set_session_name"))

#' @describeIn set_session_name Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("set_session_name", "IsoNetAnalyzeR", function(obj, value) {
  if (!is.character(value) || length(value) != 1L || nchar(trimws(value)) == 0L)
    stop("`value` must be a single non-empty character string.", call. = FALSE)
  obj@name <- value
  invisible(obj)
})


# ==============================================================================
#  Metadata accessors
# ==============================================================================

# Fields managed automatically by the constructor and other package functions.
# Users must not overwrite these directly; dedicated functions exist for that.
.protected_metadata_keys <- c("network_source", "pi_threshold", "network_cutoff")


#' List the keys present in the metadata of an IsoNetAnalyzeR object
#'
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
#' @return Character vector of metadata key names, or \code{character(0)} if
#'   none are set.
#'
#' @examples
#' \dontrun{
#' list_metadata(obj)
#' }
#'
#' @export
setGeneric("list_metadata", function(obj) standardGeneric("list_metadata"))

#' @describeIn list_metadata Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("list_metadata", "IsoNetAnalyzeR", function(obj) {
  names(obj@metadata)
})


#' Get a single metadata field from an IsoNetAnalyzeR object
#'
#' Retrieves one named entry from the \code{@@metadata} slot.
#'
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
#' @param key Character scalar. Name of the metadata field to retrieve.
#'
#' @return The value stored under \code{key}.
#'
#' @examples
#' \dontrun{
#' get_metadata(obj, "pi_threshold")
#' get_metadata(obj, "author")
#' }
#'
#' @export
setGeneric("get_metadata", function(obj, key) standardGeneric("get_metadata"))

#' @describeIn get_metadata Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("get_metadata", "IsoNetAnalyzeR", function(obj, key) {
  if (!is.character(key) || length(key) != 1L || nchar(trimws(key)) == 0L)
    stop("`key` must be a single non-empty character string.", call. = FALSE)

  if (!key %in% names(obj@metadata))
    stop(sprintf(
      "No metadata field '%s'. Available: %s",
      key, .slot_names_str(obj@metadata)
    ), call. = FALSE)

  obj@metadata[[key]]
})


#' Set a metadata field on an IsoNetAnalyzeR object
#'
#' Adds or updates a single named entry in the \code{@@metadata} slot.
#' The fields \code{network_source}, \code{pi_threshold}, and
#' \code{network_cutoff} are managed internally by the package and cannot be
#' set with this function; attempting to do so raises an error.
#'
#' @param obj   An \linkS4class{IsoNetAnalyzeR} session object.
#' @param key   Character scalar. Name of the metadata field to set.
#' @param value The value to store.  Can be any R object, though scalar
#'   values (character, numeric, logical) are recommended for readability.
#'
#' @return The updated \linkS4class{IsoNetAnalyzeR} object, invisibly.
#'
#' @examples
#' \dontrun{
#' obj <- set_metadata(obj, "author", "Jane Doe")
#' obj <- set_metadata(obj, "notes", "Pilot run, March 2025")
#' print_metadata(obj)
#' }
#'
#' @export
setGeneric("set_metadata",
           function(obj, key, value) standardGeneric("set_metadata"))

#' @describeIn set_metadata Method for \linkS4class{IsoNetAnalyzeR}.
setMethod("set_metadata", "IsoNetAnalyzeR", function(obj, key, value) {
  if (!is.character(key) || length(key) != 1L || nchar(trimws(key)) == 0L)
    stop("`key` must be a single non-empty character string.", call. = FALSE)

  if (key %in% .protected_metadata_keys)
    stop(sprintf(
      "'%s' is a protected metadata field managed internally by the package and cannot be set via set_metadata(). Use the appropriate package function to modify it.",
      key
    ), call. = FALSE)

  obj@metadata[[key]] <- value
  invisible(obj)
})
