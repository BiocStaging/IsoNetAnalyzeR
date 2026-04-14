## IsoNetAnalyzeR-metadata.R
## Generic and method for printing the metadata slot of an IsoNetAnalyzeR
## session object.
##
## Exported: print_metadata

#' @include isonetanalyzer_class.R
NULL

# ==============================================================================
#  Generic
# ==============================================================================

#' Print the metadata of an IsoNetAnalyzeR session object
#'
#' Displays every key-value pair stored in the \code{@@metadata} slot in a
#' readable, aligned format.  Values that are themselves lists or vectors are
#' shown compactly on one line.
#'
#' @param object An \linkS4class{IsoNetAnalyzeR} session object.
#' @param ... Currently unused; reserved for future arguments.
#'
#' @return \code{object}, invisibly.
#'
#' @examples
#' \dontrun{
#' obj <- IsoNetAnalyzeR("demo")
#' print_metadata(obj)
#'
#' # Add custom metadata and inspect
#' obj@@metadata$author <- "Jane Doe"
#' obj@@metadata$notes  <- "Pilot run, March 2025"
#' print_metadata(obj)
#' }
#'
#' @export
setGeneric("print_metadata", function(object, ...) standardGeneric("print_metadata"))


# ==============================================================================
#  Method
# ==============================================================================

#' @describeIn print_metadata Print metadata for an \linkS4class{IsoNetAnalyzeR}
#'   object.
setMethod("print_metadata", "IsoNetAnalyzeR", function(object, ...) {

  meta <- object@metadata

  cat("── Metadata  [", object@name, "] ─────────────\n", sep = "")

  if (length(meta) == 0L) {
    cat("  (empty)\n")
    cat("──────────────────────────────────────────\n")
    return(invisible(object))
  }

  # Determine padding width from the longest key name
  keys      <- names(meta)
  pad_width <- max(nchar(keys))

  for (key in keys) {
    val <- meta[[key]]

    val_str <- if (is.null(val)) {
      "NULL"
    } else if (length(val) == 0L) {
      "(length-0)"
    } else if (is.list(val)) {
      # Nested list: show names only, or fallback to type
      if (!is.null(names(val))) {
        paste0("list(", paste(names(val), collapse = ", "), ")")
      } else {
        paste0("list [", length(val), " element(s)]")
      }
    } else if (length(val) == 1L) {
      # Scalar — format numerics neatly
      if (is.numeric(val)) formatC(val, format = "g", digits = 6) else as.character(val)
    } else {
      # Vector: first few elements + count
      shown   <- utils::head(val, 4L)
      preview <- if (is.numeric(shown)) {
        paste(formatC(shown, format = "g", digits = 4), collapse = ", ")
      } else {
        paste(shown, collapse = ", ")
      }
      if (length(val) > 4L)
        preview <- paste0(preview, ", … [", length(val), " total]")
      preview
    }

    cat(sprintf("  %-*s : %s\n", pad_width, key, val_str))
  }

  cat("──────────────────────────────────────────\n")
  invisible(object)
})
