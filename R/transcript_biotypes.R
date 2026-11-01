#' Transcript biotypes used for protein-network eligibility
#'
#' Returns pinned GENCODE v41 transcript annotation. These transcript biotypes
#' are not gene-level coding labels or evidence of stable protein production.
#' @return A data frame with `isoform_id` and `transcript_biotype` columns.
#' @examples
#' biotypes <- isonet_transcript_biotypes()
#' head(biotypes)
#' head(sort(table(biotypes$transcript_biotype), decreasing = TRUE))
#'
#' @export
isonet_transcript_biotypes <- function() {
  if (is.null(.biotype_cache$table)) {
    path <- system.file("extdata", "transcript_biotypes_gencode_v41.csv.gz",
                        package = "IsoNetAnalyzeR")
    if (!nzchar(path))
      stop("Packaged transcript-biotype annotation is missing.", call. = FALSE)
    .biotype_cache$table <- utils::read.csv(path, stringsAsFactors = FALSE)
  }
  .biotype_cache$table
}

.biotype_cache <- new.env(parent = emptyenv())

#' Gene symbols of the bundled network
#'
#' Returns the GENCODE v41 gene name, usually the HGNC symbol, of every gene
#' in the bundled network. Plots use it to label nodes, and it can annotate
#' any result table keyed by Ensembl gene ID.
#' @return A data frame with `gene_id` (Ensembl, without version) and
#'   `gene_name` columns.
#' @examples
#' gene_names <- isonet_gene_names()
#' gene_names[gene_names$gene_id == "ENSG00000070371", ]
#'
#' @export
isonet_gene_names <- function() {
  if (is.null(.gene_name_cache$table)) {
    path <- system.file("extdata", "gene_names_gencode_v41.csv.gz",
                        package = "IsoNetAnalyzeR")
    if (!nzchar(path))
      stop("Packaged gene-name annotation is missing.", call. = FALSE)
    .gene_name_cache$table <- utils::read.csv(path, stringsAsFactors = FALSE)
  }
  .gene_name_cache$table
}

.gene_name_cache <- new.env(parent = emptyenv())

# Gene symbols for Ensembl gene IDs (versioned or not); NA when unknown.
.gene_symbols <- function(gene_ids) {
  ann <- isonet_gene_names()
  ann$gene_name[match(sub("\\..*$", "", gene_ids), ann$gene_id)]
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
  if (!is.data.frame(x) || !all(required %in% names(x)))
    stop("`transcript_biotypes` must have isoform_id and transcript_biotype columns.", call. = FALSE)
  x <- x[required]
  x$isoform_id <- sub("\\.[0-9]+$", "", trimws(as.character(x$isoform_id)))
  x$transcript_biotype <- tolower(trimws(as.character(x$transcript_biotype)))
  if (anyNA(x$isoform_id) || any(!nzchar(x$isoform_id)) ||
      anyNA(x$transcript_biotype) || any(!nzchar(x$transcript_biotype)))
    stop("Transcript-biotype annotations must have non-missing IDs and biotypes.", call. = FALSE)
  x <- unique(x)
  if (anyDuplicated(x$isoform_id))
    stop("Conflicting transcript biotypes for the same normalized ID.", call. = FALSE)
  x
}

.filter_network_biotypes <- function(network, include_nmd, include_nsd,
                                     transcript_biotypes = NULL) {
  for (arg in c("include_nmd", "include_nsd")) {
    value <- get(arg)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop(sprintf("`%s` must be TRUE or FALSE.", arg), call. = FALSE)
  }
  ann <- isonet_transcript_biotypes()
  by_id <- stats::setNames(ann$transcript_biotype, ann$isoform_id)
  ids <- igraph::V(network)$name
  stable <- sub("\\.[0-9]+$", "", ids)
  vertex_biotypes <- igraph::vertex_attr(network, "transcript_biotype")
  if (!is.null(vertex_biotypes)) {
    known <- !is.na(vertex_biotypes) & nzchar(trimws(vertex_biotypes))
    override <- .validate_biotype_table(data.frame(
      isoform_id = ids[known], transcript_biotype = vertex_biotypes[known]
    ))
    by_id[override$isoform_id] <- override$transcript_biotype
  }
  if (!is.null(transcript_biotypes)) {
    override <- .validate_biotype_table(transcript_biotypes)
    by_id[override$isoform_id] <- override$transcript_biotype
  }
  biotype <- unname(by_id[stable])
  unknown <- is.na(biotype)
  if (any(unknown) && (!include_nmd || !include_nsd))
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
  list(network = network, audit = list(
    include_nmd = include_nmd, include_nsd = include_nsd,
    annotation = "GENCODE v41; optional vertex/table overrides",
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
