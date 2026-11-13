## Transcript and gene annotation from Ensembl via AnnotationHub.
##
## The bundled network uses GENCODE v41 transcripts. Ensembl release 107 is
## the same gene build, and its EnsDb is AnnotationHub record AH104864.

.ISONET_ENSDB <- "AH104864"

.annotation_cache <- new.env(parent = emptyenv())

#' The Ensembl 107 EnsDb matching the bundled network
#'
#' Loads AnnotationHub record AH104864 once per R session. AnnotationHub
#' caches the download (about 440 MB) on disk.
#' @noRd
.isonet_ensdb <- function() {
  if (is.null(.annotation_cache$ensdb)) {
    message("Loading the Ensembl 107 (GENCODE v41) annotation from ",
            "AnnotationHub (", .ISONET_ENSDB, ") ...")
    edb <- tryCatch(
      AnnotationHub::AnnotationHub()[[.ISONET_ENSDB]],
      error = function(e) stop(
        "Could not load AnnotationHub record ", .ISONET_ENSDB, ": ",
        conditionMessage(e), call. = FALSE
      )
    )
    .annotation_cache$ensdb <- edb
  }
  .annotation_cache$ensdb
}

.check_ensdb <- function(edb) {
  if (!is(edb, "EnsDb"))
    stop("`edb` must be an EnsDb object.", call. = FALSE)
  edb
}

#' Transcript biotypes used for protein-network eligibility
#'
#' Returns the transcript biotype of every transcript in an Ensembl
#' annotation. By default this is Ensembl release 107, the gene build of
#' GENCODE v41 used by the bundled network, loaded from AnnotationHub
#' (record AH104864) and cached for the R session. These transcript
#' biotypes are not gene-level coding labels or evidence of stable protein
#' production.
#'
#' @param edb An \code{EnsDb} object, or \code{NULL} (default) for the
#'   Ensembl 107 annotation from AnnotationHub.
#' @return A \code{DataFrame} with \code{isoform_id} and
#'   \code{transcript_biotype} columns.
#' @examples
#' biotypes <- isonet_transcript_biotypes()
#' head(biotypes)
#' head(sort(table(biotypes$transcript_biotype), decreasing = TRUE))
#'
#' @export
isonet_transcript_biotypes <- function(edb = NULL) {
  if (is.null(edb)) {
    if (is.null(.annotation_cache$biotypes))
      .annotation_cache$biotypes <- isonet_transcript_biotypes(.isonet_ensdb())
    return(.annotation_cache$biotypes)
  }
  tx <- ensembldb::transcripts(.check_ensdb(edb),
                               columns = c("tx_id", "tx_biotype"),
                               return.type = "data.frame")
  S4Vectors::DataFrame(isoform_id = tx$tx_id,
                       transcript_biotype = tx$tx_biotype)
}

#' Gene symbols for Ensembl gene IDs
#'
#' Returns the gene name, usually the HGNC symbol, of every gene in an
#' Ensembl annotation; by default Ensembl release 107 from AnnotationHub
#' (record AH104864), as for \code{\link{isonet_transcript_biotypes}}.
#' Plots use it to label nodes, and it can annotate any result table keyed
#' by Ensembl gene ID.
#'
#' @inheritParams isonet_transcript_biotypes
#' @return A \code{DataFrame} with \code{gene_id} (Ensembl, without version)
#'   and \code{gene_name} columns. \code{gene_name} is \code{NA} for genes
#'   without a symbol.
#' @examples
#' gene_names <- isonet_gene_names()
#' gene_names[gene_names$gene_id == "ENSG00000070371", ]
#'
#' @export
isonet_gene_names <- function(edb = NULL) {
  if (is.null(edb)) {
    if (is.null(.annotation_cache$gene_names))
      .annotation_cache$gene_names <- isonet_gene_names(.isonet_ensdb())
    return(.annotation_cache$gene_names)
  }
  genes <- ensembldb::genes(.check_ensdb(edb),
                            columns = c("gene_id", "gene_name"),
                            return.type = "data.frame")
  name <- genes$gene_name
  name[!is.na(name) & !nzchar(name)] <- NA_character_
  S4Vectors::DataFrame(gene_id = genes$gene_id, gene_name = name)
}

# Gene symbols for Ensembl gene IDs (versioned or not); NA when unknown.
# Only Ensembl gene IDs are looked up, so networks with other identifiers
# never load the annotation; a failure to load it leaves the symbols NA and
# is not retried until the annotation has been loaded by another route.
.gene_symbols <- function(gene_ids) {
  gene_ids <- as.character(gene_ids)
  stable <- .strip_gene_version(gene_ids)
  out <- rep(NA_character_, length(gene_ids))
  ensembl <- grepl("^ENSG[0-9]+$", stable)
  if (!any(ensembl)) return(out)
  if (isTRUE(.annotation_cache$symbols_failed) &&
      is.null(.annotation_cache$ensdb))
    return(out)
  ann <- tryCatch(isonet_gene_names(), error = function(e) {
    .annotation_cache$symbols_failed <- TRUE
    NULL
  })
  if (is.null(ann)) return(out)
  out[ensembl] <- ann$gene_name[match(stable[ensembl], ann$gene_id)]
  out
}

# Gene symbol where known, otherwise the ID itself.
.gene_label <- function(gene_ids) {
  symbols <- .gene_symbols(gene_ids)
  ifelse(is.na(symbols), gene_ids, symbols)
}

# "SYMBOL (gene ID)" for titles and printed summaries, or the ID alone.
.gene_title <- function(gene_id) {
  symbol <- .gene_symbols(gene_id)
  if (is.na(symbol)) gene_id else sprintf("%s (%s)", symbol, gene_id)
}

.validate_biotype_table <- function(x) {
  required <- c("isoform_id", "transcript_biotype")
  x <- .as_data_frame(x, "transcript_biotypes")
  if (!all(required %in% names(x)))
    stop("`transcript_biotypes` must have isoform_id and transcript_biotype columns.", call. = FALSE)
  x <- x[required]
  x$isoform_id <- .strip_version(trimws(as.character(x$isoform_id)))
  x$transcript_biotype <- tolower(trimws(as.character(x$transcript_biotype)))
  if (anyNA(x$isoform_id) || any(!nzchar(x$isoform_id)) ||
      anyNA(x$transcript_biotype) || any(!nzchar(x$transcript_biotype)))
    stop("Transcript-biotype annotations must have non-missing IDs and biotypes.", call. = FALSE)
  x <- unique(x)
  if (anyDuplicated(x$isoform_id))
    stop("Conflicting transcript biotypes for the same normalized ID.", call. = FALSE)
  x
}

#' Annotate transcript biotypes and exclude decay transcripts
#'
#' Biotypes come, in increasing priority, from an EnsDb (the
#' \code{transcript_biotypes} argument, or the Ensembl 107 annotation when
#' the network has Ensembl transcript IDs), from a
#' \code{transcript_biotype} vertex attribute, and from a
#' \code{transcript_biotypes} table. Gene symbols from the EnsDb are added
#' as a \code{gene_name} vertex attribute when the network has none.
#' Networks without Ensembl transcript IDs do not load the Ensembl 107
#' annotation here; without any exclusion it is optional, and a failure to load
#' it leaves the vertices unannotated.
#' @noRd
.filter_network_biotypes <- function(network, include_nmd, include_nsd,
                                     transcript_biotypes = NULL) {
  for (arg in c("include_nmd", "include_nsd")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop(sprintf("`%s` must be TRUE or FALSE.", arg), call. = FALSE)
  }
  ids <- igraph::V(network)$name
  stable <- .strip_version(ids)
  user_edb <- is(transcript_biotypes, "EnsDb")
  table <- if (!is.null(transcript_biotypes) && !user_edb)
    .validate_biotype_table(transcript_biotypes)

  by_id <- stats::setNames(character(), character())
  vertex_biotypes <- igraph::vertex_attr(network, "transcript_biotype")
  if (!is.null(vertex_biotypes)) {
    known <- !is.na(vertex_biotypes) & nzchar(trimws(vertex_biotypes))
    vertex_table <- .validate_biotype_table(data.frame(
      isoform_id = ids[known], transcript_biotype = vertex_biotypes[known]
    ))
    by_id[vertex_table$isoform_id] <- vertex_table$transcript_biotype
  }
  if (!is.null(table))
    by_id[table$isoform_id] <- table$transcript_biotype

  filtering <- !include_nmd || !include_nsd
  # Ensembl transcript IDs are always annotated from Ensembl, so that every
  # annotated decay transcript, in the network or not, is excluded. Without
  # it, a session needs no exclusion or a biotype for every vertex.
  edb <- if (user_edb) {
    transcript_biotypes
  } else if (any(grepl("^ENST[0-9]+", stable))) {
    tryCatch(.isonet_ensdb(), error = function(e) {
      # Labels and symbol lookups later in the session do not retry either
      .annotation_cache$symbols_failed <- TRUE
      if (!filtering) {
        message("The Ensembl 107 annotation could not be loaded; transcript ",
                "biotypes and gene symbols are not annotated.")
      } else if (all(stable %in% names(by_id))) {
        warning(conditionMessage(e), "\nEvery vertex has a biotype, so the ",
                "session is created, but decay transcripts outside the ",
                "network cannot be blocked during identifier resolution.",
                call. = FALSE)
      } else {
        stop(conditionMessage(e), "\nTo create the session without ",
             "AnnotationHub, give every vertex a biotype (a ",
             "`transcript_biotype` vertex attribute or a ",
             "`transcript_biotypes` table), pass the annotation as an EnsDb ",
             "through `transcript_biotypes`, or set include_nmd = TRUE and ",
             "include_nsd = TRUE.", call. = FALSE)
      }
      NULL
    })
  }
  if (!is.null(edb)) {
    ann <- if (user_edb) isonet_transcript_biotypes(edb) else
      isonet_transcript_biotypes()
    # Every annotated transcript, so that seeds outside the network that are
    # decay transcripts are still reported as excluded
    base <- stats::setNames(ann$transcript_biotype,
                            .strip_version(ann$isoform_id))
    by_id <- c(by_id, base[!names(base) %in% names(by_id)])
    if (is.null(igraph::vertex_attr(network, "gene_name"))) {
      genes <- if (user_edb) isonet_gene_names(edb) else isonet_gene_names()
      vertex_genes <- .strip_gene_version(igraph::V(network)$gene)
      network <- igraph::set_vertex_attr(
        network, "gene_name",
        value = genes$gene_name[match(vertex_genes,
                                      .strip_gene_version(genes$gene_id))]
      )
    }
  }

  biotype <- unname(by_id[stable])
  unknown <- is.na(biotype)
  if (any(unknown) && filtering)
    stop(
      sprintf("%d network vertices have unknown transcript biotypes. Supply a ",
              sum(unknown)),
      "`transcript_biotype` vertex attribute or `transcript_biotypes` table; ",
      "gene biotypes are not sufficient. To explicitly use an unfiltered ",
      "network, set include_nmd = TRUE and include_nsd = TRUE.",
      call. = FALSE
    )
  excluded_types <- c(if (!include_nmd) "nonsense_mediated_decay",
                      if (!include_nsd) "non_stop_decay")
  excluded_ids <- names(by_id)[by_id %in% excluded_types]
  excluded <- stable %in% excluded_ids
  audit <- data.frame(isoform_id = ids, transcript_biotype = biotype,
                       excluded = excluded, stringsAsFactors = FALSE)
  network <- igraph::set_vertex_attr(network, "transcript_biotype", value = biotype)
  before_edges <- igraph::ecount(network)
  if (any(excluded)) network <- igraph::delete_vertices(network, which(excluded))
  network <- igraph::set_graph_attr(network, "excluded_transcript_ids", excluded_ids)
  annotation <- c(
    if (user_edb) "user-supplied EnsDb" else if (!is.null(edb))
      "Ensembl 107 (GENCODE v41), AnnotationHub AH104864",
    if (!is.null(vertex_biotypes)) "transcript_biotype vertex attribute",
    if (!is.null(table)) "transcript_biotypes table"
  )
  annotation <- if (length(annotation)) paste(annotation, collapse = "; ") else
    "none"
  list(network = network, audit = list(
    include_nmd = include_nmd, include_nsd = include_nsd,
    annotation = annotation,
    vertices_before = length(ids), vertices_after = igraph::vcount(network),
    edges_before = before_edges, edges_after = igraph::ecount(network),
    excluded_nmd = sum(excluded & biotype == "nonsense_mediated_decay", na.rm = TRUE),
    excluded_nsd = sum(excluded & biotype == "non_stop_decay", na.rm = TRUE),
    unknown_vertices = sum(unknown), vertex_audit = audit
  ))
}

.excluded_decay_ids <- function(network) {
  if (is.null(network)) return(character())
  ids <- igraph::graph_attr(network, "excluded_transcript_ids")
  if (is.null(ids)) character() else as.character(ids)
}
