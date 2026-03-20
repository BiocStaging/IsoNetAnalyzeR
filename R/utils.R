#' Resolve the graph argument
#'
#' Returns \code{g} unchanged if provided, otherwise falls back to the graph
#' stored in the package environment. Raises an error when neither
#' is available.
#'
#' @param g An igraph object or \code{NULL}.
#' @return An igraph object.
#' @noRd
.resolve_graph <- function(g) {
  if (!is.null(g)) {
    if (!inherits(g, "igraph")) {
      stop("`g` must be an igraph object, got: ", class(g)[1], call. = FALSE)
    }
    return(g)
  }

  # isoform_network_graph is available via R/sysdata.rda (internal = TRUE)
  isoform_network_graph
}
