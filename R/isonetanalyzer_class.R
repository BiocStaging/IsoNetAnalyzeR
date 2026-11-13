## IsoNetAnalyzeR-class.R
## S4 class definition for the main IsoNetAnalyzeR session object.
##
## Exported: IsoNetAnalyzeR (constructor), show method,
##           get_network_diagnostics

setOldClass("igraph")
#' @importFrom methods is new setClass setGeneric setMethod setValidity validObject
#' @importFrom rlang .data
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame
#' @importFrom patchwork plot_annotation
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
#' @return Objects of this class are created by \code{\link{IsoNetAnalyzeR}}.
#'
#' @examples
#' obj <- IsoNetAnalyzeR("example")
#' obj
#' list_analyses(obj)
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
#'   Repeated edges between the same two transcripts are merged into one
#'   (\code{interact} is 1 when any copy is 1; numeric attributes such as
#'   \code{pi} are averaged over the predicted copies, or over all copies
#'   when none is predicted) and self-loops are removed, because they would
#'   inflate the degrees used by the neighbourhood cutoff.
#' @param metadata Named list of additional metadata to attach to the session.
#'   This is merged with any automatically populated fields
#'   (\code{pi_threshold}, \code{network_cutoff}) - user-supplied values take
#'   precedence on name collision.
#'
#' @param include_nmd Include nonsense-mediated decay transcripts. Default FALSE.
#' @param include_nsd Include non-stop decay transcripts. Default FALSE.
#'   Exclusions apply to seeds and neighbours. Set both flags TRUE for an
#'   unfiltered network.
#' @param transcript_biotypes Optional transcript annotation: a
#'   \code{data.frame} or \code{DataFrame} with \code{isoform_id} and
#'   \code{transcript_biotype} columns, which overrides vertex annotations,
#'   or an \code{EnsDb} object. A network with Ensembl transcript IDs
#'   (\code{ENST}) is always annotated from Ensembl release 107 (the GENCODE
#'   v41 gene build of the bundled network), loaded once per R session from
#'   AnnotationHub record AH104864 and cached on disk by AnnotationHub, or
#'   from the \code{EnsDb} given here; a \code{transcript_biotype} vertex
#'   attribute and a table override it for their transcripts. The full
#'   annotation is needed so that seeds outside the network that are decay
#'   transcripts are excluded too. If it cannot be loaded (for example
#'   offline), a session is still created when both \code{include_nmd} and
#'   \code{include_nsd} are \code{TRUE} (without biotypes and gene symbols),
#'   or when every vertex has a biotype from the attribute or table (with a
#'   warning that decay transcripts outside the network cannot be blocked).
#'   Networks without Ensembl transcript IDs do not load the annotation
#'   when the session is created; gene-symbol labels of Ensembl gene IDs in
#'   plots and printed summaries can still load it.
#'   Unknown vertex biotypes are an error when either decay exclusion is
#'   enabled. Gene-level biotypes are not used.
#' @return An \linkS4class{IsoNetAnalyzeR} object. The exclusion audit is stored
#'   in metadata$network_biotype_filter. Create a new session to change the
#'   policy. When biotypes come from an EnsDb, the network also gets a
#'   \code{gene_name} vertex attribute with gene symbols, used for plot
#'   labels.
#'
#' @examples
#' # Session on the bundled network (NMD and NSD transcripts excluded)
#' obj <- IsoNetAnalyzeR("example")
#' obj
#'
#' # A user-supplied network needs vertex names, a `gene` vertex attribute,
#' # a binary `interact` edge attribute and known transcript biotypes
#' g <- igraph::graph_from_data_frame(
#'   data.frame(from = "tx1", to = "tx2", interact = 1L, pi = 0.9),
#'   directed = FALSE,
#'   vertices = data.frame(name = c("tx1", "tx2"), gene = c("g1", "g2"),
#'                         transcript_biotype = "protein_coding")
#' )
#' IsoNetAnalyzeR("custom", network = g, metadata = list(network_cutoff = 0))
#'
#' @export
IsoNetAnalyzeR <- function(name     = "IsoNetAnalyzeR",
                           network  = NULL,
                           metadata = list(),
                           include_nmd = FALSE,
                           include_nsd = FALSE,
                           transcript_biotypes = NULL) {

  user_supplied_network <- !is.null(network)

  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !nzchar(trimws(name)))
    stop("`name` must be a single non-empty character string.", call. = FALSE)

  if (!is.list(metadata) ||
      (length(metadata) > 0L && (is.null(names(metadata)) ||
       any(is.na(names(metadata))) || any(!nzchar(names(metadata))))))
    stop("`metadata` must be a named list.", call. = FALSE)

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
      "Loaded the bundled network: %s transcripts, %s edges.",
      .fmt_int(igraph::vcount(network)), .fmt_int(igraph::ecount(network))
    ))
  } else {
    if (!inherits(network, "igraph"))
      stop("`network` must be an igraph object.", call. = FALSE)
    auto_meta <- list(network_source = "user_supplied")
    message(sprintf(
      "Using the supplied network: %s transcripts, %s edges.",
      .fmt_int(igraph::vcount(network)), .fmt_int(igraph::ecount(network))
    ))
  }

  if (user_supplied_network) {
    .validate_network_contract(network)
    network <- .simplify_network(network)
  }

  # User metadata wins on collision
  merged_meta <- utils::modifyList(auto_meta, metadata)

  if (!is.null(merged_meta[["network_cutoff"]]) &&
      !(!user_supplied_network &&
        length(merged_meta[["network_cutoff"]]) == 1L &&
        is.na(merged_meta[["network_cutoff"]])))
    .validate_cutoff(merged_meta[["network_cutoff"]], "metadata$network_cutoff")

  filtered <- .filter_network_biotypes(
    network, include_nmd, include_nsd, transcript_biotypes
  )
  network <- filtered$network
  merged_meta$network_biotype_filter <- filtered$audit
  message(
    sprintf("Excluded %s NMD and %s NSD transcripts with their edges;\n",
            .fmt_int(filtered$audit$excluded_nmd),
            .fmt_int(filtered$audit$excluded_nsd)),
    sprintf("%s transcripts and %s edges remain.",
            .fmt_int(igraph::vcount(network)),
            .fmt_int(igraph::ecount(network)))
  )

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
  source <- object@metadata[["network_source"]] %||% "unknown"
  source <- switch(source, package_internal = "bundled",
                   user_supplied = "user-supplied", source)
  n_vertices <- .fmt_int(igraph::vcount(object@network))
  n_edges    <- .fmt_int(igraph::ecount(object@network))
  diag_status <- if (is.null(object@network_diagnostics)) "not run" else
    format(object@network_diagnostics@created_at, '%Y-%m-%d %H:%M:%S')

  cli::cli_rule("IsoNetAnalyzeR session")
  .cli_aligned("Name        : {object@name}")
  .cli_aligned("Created     : {format(object@created_at, '%Y-%m-%d %H:%M:%S')}")
  .cli_aligned("Network     : {n_vertices} transcripts, {n_edges} edges ({source})")
  .cli_aligned("Analyses    : {(.slot_names_str(object@analyses))}")
  .cli_aligned("Lookups     : {(.slot_names_str(object@lookups))}")
  .cli_aligned("Diagnostics : {diag_status}")

  extra <- setdiff(names(object@metadata), "network_source")
  if (length(extra)) {
    keys <- strwrap(paste(extra, collapse = ", "),
                    width = cli::console_width() - 15L)
    .cli_aligned("Metadata    : {keys[1L]}")
    for (line in keys[-1L]) .cli_aligned("              {line}")
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

#' Format counts with thousands separators
#' @noRd
.fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")

#' cli::cli_text() that keeps column alignment
#'
#' cli_text() collapses repeated spaces, which misaligns the padded labels of
#' printed summaries. Runs of spaces outside glue braces become non-breaking
#' spaces, which cli prints as ordinary spaces. A line break inside the text,
#' with the indentation that follows it, still reads as one space.
#' @noRd
.cli_aligned <- function(..., .envir = parent.frame()) {
  text  <- gsub("\n[ \t]*", " ", paste0(...))
  chars <- strsplit(text, "")[[1L]]
  depth <- cumsum(chars == "{") - cumsum(chars == "}")
  space <- chars == " " & depth == 0L
  n <- length(space)
  run <- space & (c(FALSE, space[-n]) | c(space[-1L], FALSE))
  chars[run] <- "\u00a0"
  cli::cli_text(paste(chars, collapse = ""), .envir = .envir)
}

#' Variable name of the session object, for "Access with" hints
#'
#' Takes \code{substitute(obj)} from the calling function and returns the
#' caller's variable name, or \code{"obj"} when the argument was an
#' expression.
#' @noRd
.obj_name <- function(expr) if (is.name(expr)) as.character(expr) else "obj"

#' Return the package name from within the package namespace
#' @noRd
packageName <- function() "IsoNetAnalyzeR"

#' Remove an Ensembl version suffix (ENST00000621271.4 -> ENST00000621271)
#' @noRd
.strip_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))

#' A gene ID for matching: cut at the first dot, so that GENCODE's
#' ENSG00000002586.20_PAR_Y matches ENSG00000002586 too
#' @noRd
.strip_gene_version <- function(x) sub("\\..*$", "", as.character(x))

#' Validate a positive whole-number scalar
#' @noRd
.validate_positive_integer <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x > .Machine$integer.max || x != floor(x))
    stop(sprintf("`%s` must be a finite positive integer scalar.", arg),
         call. = FALSE)
  as.integer(x)
}

#' Validate a probability-like cutoff
#' @noRd
.validate_cutoff <- function(x, arg = "cutoff") {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 0 || x > 1)
    stop(sprintf("`%s` must be a single finite numeric value in [0, 1].", arg),
         call. = FALSE)
  as.numeric(x)
}

#' Validate the internal network contract
#' @noRd
.validate_network_contract <- function(network) {
  if (!inherits(network, "igraph"))
    stop("`network` must be an igraph object.", call. = FALSE)

  if (igraph::is_directed(network))
    stop("`network` must be undirected.", call. = FALSE)

  if (igraph::vcount(network) == 0L)
    stop("`network` must contain at least one vertex.", call. = FALSE)

  vertex_names <- igraph::V(network)$name
  if (is.null(vertex_names) || length(vertex_names) != igraph::vcount(network) ||
      anyNA(vertex_names) || any(!nzchar(trimws(vertex_names))) ||
      anyDuplicated(vertex_names))
    stop("`network` vertex names must be unique, non-missing, and non-empty.",
         call. = FALSE)

  genes <- igraph::vertex_attr(network, "gene")
  if (is.null(genes) || length(genes) != igraph::vcount(network) ||
      anyNA(genes) || any(!nzchar(trimws(as.character(genes)))))
    stop("`network` must have a non-missing, non-empty `gene` attribute for every vertex.",
         call. = FALSE)

  if (!"interact" %in% igraph::edge_attr_names(network))
    stop("`network` must have a binary `interact` edge attribute.", call. = FALSE)

  interact <- igraph::edge_attr(network, "interact")
  if (!(is.numeric(interact) || is.logical(interact)) || anyNA(interact) ||
      any(!is.finite(as.numeric(interact))) ||
      any(!as.numeric(interact) %in% c(0, 1)))
    stop("`network` edge attribute `interact` must contain only non-missing binary values (0 or 1).",
         call. = FALSE)

  invisible(TRUE)
}

#' Merge repeated edges and remove self-loops
#'
#' Neighbourhood cutoffs compare vertex degrees, which repeated edges and
#' self-loops would inflate. One edge is kept per transcript pair:
#' \code{interact} is 1 when any copy is 1, other numeric attributes are
#' averaged and the rest keep their first value. A pair with a predicted copy
#' keeps only its predicted copies, so its \code{pi} is the mean of those,
#' as in inst/scripts/reduce_network.R.
#' @noRd
.simplify_network <- function(network) {
  if (igraph::is_simple(network)) return(network)
  loops <- igraph::which_loop(network)
  n_repeated <- sum(igraph::which_multiple(network) & !loops)
  ends <- igraph::as_edgelist(network, names = FALSE)
  pair <- paste(pmin(ends[, 1L], ends[, 2L]), pmax(ends[, 1L], ends[, 2L]))
  predicted <- as.numeric(igraph::edge_attr(network, "interact")) == 1
  network <- igraph::delete_edges(
    network, which(!predicted & pair %in% pair[predicted])
  )
  attrs <- igraph::edge_attr_names(network)
  combine <- lapply(attrs, function(a) {
    if (a == "interact") return(function(x) max(as.numeric(x)))
    if (is.numeric(igraph::edge_attr(network, a))) "mean" else "first"
  })
  names(combine) <- attrs
  network <- igraph::simplify(network, remove.multiple = TRUE,
                              remove.loops = TRUE, edge.attr.comb = combine)
  message(sprintf(
    "Merged %s repeated edge(s) and removed %s self-loop(s) from the network.",
    .fmt_int(n_repeated), .fmt_int(sum(loops))
  ))
  network
}

#' Coerce a tabular input to a data.frame
#'
#' Accepts a \code{data.frame} (including tibbles) or an S4Vectors
#' \code{DataFrame}; column names are kept as they are.
#' @noRd
.as_data_frame <- function(x, arg) {
  if (is(x, "DataFrame"))
    return(as.data.frame(x, optional = TRUE))
  if (is.data.frame(x))
    return(as.data.frame(x, stringsAsFactors = FALSE))
  stop(sprintf("`%s` must be a data.frame or DataFrame.", arg), call. = FALSE)
}

#' Is a tabular input (data.frame, tibble or DataFrame)?
#' @noRd
.is_table <- function(x) is.data.frame(x) || is(x, "DataFrame")

#' Store a table as an S4Vectors DataFrame
#' @noRd
.as_DataFrame <- function(x) {
  if (is(x, "DataFrame")) return(x)
  S4Vectors::DataFrame(x, check.names = FALSE)
}

#' A stored DataFrame slot as a data.frame, for internal computation
#' @noRd
.df <- function(x) {
  if (is(x, "DataFrame")) as.data.frame(x, optional = TRUE) else x
}


# ==============================================================================
#  Accessor generics and methods
# ==============================================================================

#' List stored analysis labels
#'
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
#' @return Character vector of analysis labels, or \code{character(0)} if none.
#' @examples
#' seeds <- read.csv(system.file("extdata", "gse141432_t2d_vs_lean_seeds.csv",
#'                               package = "IsoNetAnalyzeR"))
#' obj <- IsoNetAnalyzeR("example")
#' obj <- run_differential_pulldown(
#'   obj, "t2d", seeds[seeds$gene_name == "CLTCL1", ],
#'   order = 1, confirm = FALSE, verbose = FALSE
#' )
#' list_analyses(obj)
#'
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
#' @examples
#' seeds <- read.csv(system.file("extdata", "gse141432_t2d_vs_lean_seeds.csv",
#'                               package = "IsoNetAnalyzeR"))
#' obj <- IsoNetAnalyzeR("example")
#' obj <- run_differential_pulldown(
#'   obj, "t2d", seeds[seeds$gene_name == "CLTCL1", ],
#'   order = 1, confirm = FALSE, verbose = FALSE
#' )
#' list_analyses(obj)
#' get_analysis(obj, "t2d")
#'
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
#' @examples
#' obj <- IsoNetAnalyzeR("example")
#' obj <- lookup_virtual_pulldown(obj, "cltcl1_up", seed = "ENST00000621271",
#'                                order = 1)
#' list_lookups(obj)
#'
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
#' @examples
#' obj <- IsoNetAnalyzeR("example")
#' obj <- lookup_virtual_pulldown(obj, "cltcl1_up", seed = "ENST00000621271",
#'                                order = 1)
#' get_lookup(obj, "cltcl1_up")
#'
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
#' # A small deterministic network keeps this example fast; on the bundled
#' # network the diagnostics take about half a minute.
#' g <- igraph::make_graph("Zachary")
#' n <- igraph::vcount(g)
#' igraph::V(g)$name <- paste0("tx", seq_len(n))
#' igraph::V(g)$gene <- paste0("gene", ceiling(seq_len(n) / 2))
#' igraph::V(g)$transcript_biotype <- "protein_coding"
#' igraph::E(g)$interact <- 1L
#' toy <- IsoNetAnalyzeR("toy", network = g,
#'                       metadata = list(network_cutoff = 0.25))
#' toy <- run_network_diagnostics(toy)
#' get_network_diagnostics(toy)
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


# ==============================================================================
#  set_metadata
# ==============================================================================

#' Set a metadata field on an IsoNetAnalyzeR session object
#'
#' Adds or updates a single named entry in the \code{@@metadata} slot of an
#' \linkS4class{IsoNetAnalyzeR} session object.  If the key already exists its
#' value is silently overwritten; if it does not exist it is created.
#'
#' The following keys are set automatically by the constructor and should
#' not normally be overwritten:
#' \code{network_source}, \code{pi_threshold}, \code{network_cutoff}.
#'
#' @param obj   An \linkS4class{IsoNetAnalyzeR} session object.
#' @param key   Character scalar. The metadata field name.
#' @param value Any R object to store under \code{key}.
#'
#' @return An updated \linkS4class{IsoNetAnalyzeR} object (invisibly).
#'
#' @examples
#' obj <- IsoNetAnalyzeR("example")
#' obj <- set_metadata(obj, "author", "Jane Doe")
#' obj <- set_metadata(obj, "notes", "Pilot run")
#' get_metadata(obj, "author")
#'
#' @export
set_metadata <- function(obj, key, value) {

  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.character(key) || length(key) != 1L || nchar(trimws(key)) == 0L)
    stop("`key` must be a single non-empty character string.", call. = FALSE)

  obj@metadata[[key]] <- value
  invisible(obj)
}


# ==============================================================================
#  get_metadata
# ==============================================================================

#' Retrieve metadata from an IsoNetAnalyzeR session object
#'
#' Returns one named entry of the session metadata, or the complete metadata
#' list. The constructor and the analysis functions record provenance here:
#' \code{network_source}, \code{network_cutoff} and
#' \code{network_biotype_filter} (the transcript eligibility audit) are set
#' when the session is created, and \code{analysis_network_coverage},
#' \code{analysis_seed_support} and \code{analysis_duplicate_resolution} are
#' lists keyed by analysis label.
#'
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
#' @param key Character scalar naming one metadata entry, or \code{NULL}
#'   (default) to return all entries.
#'
#' @return The value stored under \code{key}, or the named list of all
#'   metadata when \code{key} is \code{NULL}.
#'
#' @seealso \code{\link{set_metadata}}, \code{\link{print_metadata}}
#'
#' @examples
#' obj <- IsoNetAnalyzeR("example")
#' names(get_metadata(obj))
#' get_metadata(obj, "network_cutoff")
#' get_metadata(obj, "network_biotype_filter")$excluded_nmd
#'
#' @export
get_metadata <- function(obj, key = NULL) {

  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (is.null(key))
    return(obj@metadata)

  if (!is.character(key) || length(key) != 1L || is.na(key) ||
      !nzchar(trimws(key)))
    stop("`key` must be a single non-empty character string.", call. = FALSE)

  if (!key %in% names(obj@metadata))
    stop(sprintf("No metadata entry '%s'. Available: %s",
                 key, .slot_names_str(obj@metadata)), call. = FALSE)

  obj@metadata[[key]]
}


# ==============================================================================
#  get_network
# ==============================================================================

#' Retrieve the interaction network of an IsoNetAnalyzeR session
#'
#' Returns the isoform interaction network used by the session, after
#' transcript eligibility filtering.
#'
#' @param obj An \linkS4class{IsoNetAnalyzeR} session object.
#'
#' @return An undirected \code{igraph} object. Vertices are transcripts with
#'   \code{name}, \code{gene} and \code{transcript_biotype} attributes; edges
#'   carry a binary \code{interact} attribute and, for the bundled network,
#'   the predicted interaction probability \code{pi}.
#'
#' @examples
#' obj <- IsoNetAnalyzeR("example")
#' net <- get_network(obj)
#' c(vertices = igraph::vcount(net), edges = igraph::ecount(net))
#'
#' @export
get_network <- function(obj) {

  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  obj@network
}
