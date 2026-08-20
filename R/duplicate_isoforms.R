## Protein-identical isoform identifier resolution and count aggregation.

.dup_map_cache <- new.env(parent = emptyenv())

.normalize_transcript_ids <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

.validate_duplicate_map <- function(map) {
  required <- c(
    "kept_ENSG", "redundant_ENSG", "kept_ENSP", "kept_ENST",
    "redundant_ENSP", "redundant_ENST", "network_applicable"
  )
  if (!is.data.frame(map) || !all(required %in% names(map))) {
    stop(
      "`duplicate_map` must contain: ", paste(required, collapse = ", "),
      ".", call. = FALSE
    )
  }
  map <- map[, required, drop = FALSE]
  for (column in setdiff(required, "network_applicable")) {
    map[[column]] <- as.character(map[[column]])
  }
  map$network_applicable <- as.logical(map$network_applicable)
  if (anyNA(map) || any(!nzchar(as.matrix(map[setdiff(
      required, "network_applicable"
    )])))) {
    stop("`duplicate_map` contains missing or blank values.", call. = FALSE)
  }
  if (anyDuplicated(map$redundant_ENST)) {
    stop("`duplicate_map` has ambiguous redundant ENSTs.", call. = FALSE)
  }
  map
}

.read_duplicate_map <- function(refresh = FALSE) {
  if (!refresh && exists("map", envir = .dup_map_cache, inherits = FALSE)) {
    return(get("map", envir = .dup_map_cache, inherits = FALSE))
  }
  map_path <- system.file(
    "extdata", "protein_identical_isoform_map.csv",
    package = "IsoNetAnalyzeR"
  )
  provenance_path <- system.file(
    "extdata", "protein_identical_isoform_map_provenance.csv",
    package = "IsoNetAnalyzeR"
  )
  if (!nzchar(map_path) || !file.exists(map_path)) {
    stop(
      "The packaged protein-identical isoform map could not be found. ",
      "Reinstall IsoNetAnalyzeR from a complete source package.",
      call. = FALSE
    )
  }
  map <- utils::read.csv(
    map_path, colClasses = c(rep("character", 6L), "logical"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  map <- .validate_duplicate_map(map)
  provenance <- if (nzchar(provenance_path) && file.exists(provenance_path)) {
    x <- utils::read.csv(
      provenance_path, colClasses = "character", stringsAsFactors = FALSE
    )
    stats::setNames(as.list(x$value), x$key)
  } else {
    list()
  }
  attr(map, "provenance") <- provenance
  assign("map", map, envir = .dup_map_cache)
  map
}

#' Protein-identical isoform identifier map
#'
#' Returns the transcriptome-wide annotation used by IsoNetAnalyzeR to map
#' redundant Ensembl transcript identifiers to the representative transcript
#' retained for the same protein product. Rows outside the embedded network are
#' retained because they are still needed for upstream count aggregation.
#'
#' @param refresh Logical; reread the packaged CSV instead of using the
#'   package-local cache.
#'
#' @return A data frame with kept/redundant ENSG, ENSP, and ENST identifiers
#'   and a \code{network_applicable} flag. Build provenance is attached as the
#'   \code{provenance} attribute.
#' @export
isonet_duplicate_map <- function(refresh = FALSE) {
  if (!is.logical(refresh) || length(refresh) != 1L || is.na(refresh)) {
    stop("`refresh` must be a single non-NA logical value.", call. = FALSE)
  }
  .read_duplicate_map(refresh = refresh)
}

.network_resolution_context <- function(network) {
  if (is.null(network)) return(NULL)
  if (!inherits(network, "igraph")) {
    stop("`network` must be an igraph object or NULL.", call. = FALSE)
  }
  vertices <- igraph::V(network)$name
  genes <- igraph::V(network)$gene
  list(
    vertices = vertices,
    gene = stats::setNames(as.character(genes), vertices)
  )
}

.resolve_duplicate_vector <- function(ids, genes = NULL, network = NULL,
                                      duplicate_map = NULL) {
  if (!is.character(ids) || anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop("Transcript IDs must be non-missing, non-blank character values.",
         call. = FALSE)
  }
  ids <- trimws(ids)
  if (!is.null(genes)) {
    genes <- trimws(as.character(genes))
    if (length(genes) != length(ids)) {
      stop("Gene and transcript vectors must have equal length.", call. = FALSE)
    }
  }
  map <- .validate_duplicate_map(
    if (is.null(duplicate_map)) isonet_duplicate_map() else duplicate_map
  )
  context <- .network_resolution_context(network)
  normalized <- .normalize_transcript_ids(ids)
  red_index <- match(normalized, map$redundant_ENST)
  kept_index <- match(normalized, map$kept_ENST)
  resolved <- ids
  status <- rep("unchanged_unmapped", length(ids))
  kept_gene <- rep(NA_character_, length(ids))

  for (i in seq_along(ids)) {
    original <- ids[[i]]
    base <- normalized[[i]]
    red <- red_index[[i]]
    kept <- kept_index[[i]]

    if (!is.null(context) && original %in% context$vertices) {
      status[[i]] <- "direct"
      next
    }
    if (!is.null(context) && base %in% context$vertices) {
      resolved[[i]] <- base
      status[[i]] <- if (identical(original, base)) "direct" else
        "version_normalized"
      next
    }

    if (!is.na(red)) {
      target <- map$kept_ENST[[red]]
      target_gene <- map$kept_ENSG[[red]]
      kept_gene[[i]] <- target_gene
      if (!is.null(genes) && !is.na(genes[[i]]) && nzchar(genes[[i]]) &&
          !identical(genes[[i]], map$redundant_ENSG[[red]])) {
        stop(
          "Submitted redundant transcript '", original, "' belongs to gene '",
          genes[[i]], "', but the duplicate map assigns it to '",
          map$redundant_ENSG[[red]], "'.", call. = FALSE
        )
      }
      if (is.null(context)) {
        resolved[[i]] <- target
        status[[i]] <- "duplicate_remapped"
      } else if (!(base %in% context$vertices) &&
                 target %in% context$vertices &&
                 identical(unname(context$gene[[target]]), target_gene)) {
        resolved[[i]] <- target
        status[[i]] <- "duplicate_remapped"
      } else if (!(target %in% context$vertices)) {
        status[[i]] <- "map_target_unavailable"
      } else {
        status[[i]] <- "map_network_incompatible"
      }
      next
    }

    if (!is.na(kept)) {
      target_gene <- map$kept_ENSG[[kept]]
      kept_gene[[i]] <- target_gene
      if (!is.null(genes) && !is.na(genes[[i]]) && nzchar(genes[[i]]) &&
          !identical(genes[[i]], target_gene)) {
        stop(
          "Submitted representative transcript '", original,
          "' disagrees with duplicate-map gene '", target_gene, "'.",
          call. = FALSE
        )
      }
      if (is.null(context)) {
        resolved[[i]] <- base
        status[[i]] <- if (identical(original, base)) "direct" else
          "version_normalized"
      } else {
        status[[i]] <- "map_target_unavailable"
      }
    }
  }

  data.frame(
    row_index = seq_along(ids),
    submitted_id = ids,
    normalized_id = normalized,
    resolved_id = resolved,
    submitted_gene = if (is.null(genes)) rep(NA_character_, length(ids)) else genes,
    kept_gene = kept_gene,
    resolution = status,
    stringsAsFactors = FALSE
  )
}

#' Resolve protein-identical isoform identifiers
#'
#' Remaps redundant Ensembl transcript identifiers to their kept protein
#' representative. With a network supplied, a mapping is applied only when it
#' is compatible with that graph. Identifiers absent from the map are retained.
#'
#' @param x Character vector or data frame containing transcript identifiers.
#' @param network Optional igraph used to validate whether each mapping is
#'   applicable.
#' @param isoform_col Name of the transcript-ID column when \code{x} is a data
#'   frame.
#' @param gene_col Name of the gene-ID column when \code{x} is a data frame.
#' @param duplicate_map Optional map with the schema returned by
#'   \code{isonet_duplicate_map()}; primarily useful for validated custom maps
#'   and tests.
#'
#' @return The same broad type as \code{x}, with resolved identifiers and a
#'   \code{dedup_audit} attribute.
#' @export
resolve_duplicate_isoforms <- function(
    x, network = NULL, isoform_col = "isoform_id", gene_col = "gene_id",
    duplicate_map = NULL) {
  if (is.character(x)) {
    audit <- .resolve_duplicate_vector(
      x, network = network, duplicate_map = duplicate_map
    )
    out <- audit$resolved_id
    names(out) <- names(x)
    attr(out, "dedup_audit") <- audit
    return(out)
  }
  if (!is.data.frame(x)) {
    stop("`x` must be a character vector or data frame.", call. = FALSE)
  }
  if (!all(c(isoform_col, gene_col) %in% names(x))) {
    stop("The requested isoform and gene columns are not both present.",
         call. = FALSE)
  }
  audit <- .resolve_duplicate_vector(
    as.character(x[[isoform_col]]), as.character(x[[gene_col]]),
    network = network, duplicate_map = duplicate_map
  )
  out <- x
  out[[isoform_col]] <- audit$resolved_id
  attr(out, "dedup_audit") <- audit
  out
}

.warn_duplicate_resolution <- function(audit) {
  if (is.null(audit) || !nrow(audit)) return(invisible(NULL))
  incompatible <- audit$resolution %in% c(
    "map_target_unavailable", "map_network_incompatible"
  )
  if (any(incompatible)) {
    warning(
      sprintf(
        paste0(
          "%d duplicate-map candidate(s) were incompatible with the active ",
          "network and were left unchanged; inspect the resolution audit."
        ),
        sum(incompatible)
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.feature_data_from_counts <- function(counts, feature_data, map) {
  rn <- rownames(counts)
  if (is.null(rn) || anyNA(rn) || any(!nzchar(rn)) || anyDuplicated(rn)) {
    stop("Counts must have unique, non-blank row names.", call. = FALSE)
  }
  if (!is.null(feature_data)) {
    feature_data <- as.data.frame(feature_data, stringsAsFactors = FALSE)
    if (nrow(feature_data) != nrow(counts) ||
        !all(c("gene_id", "isoform_id") %in% names(feature_data))) {
      stop(
        "`feature_data` must align to count rows and contain gene_id and ",
        "isoform_id.", call. = FALSE
      )
    }
    rownames(feature_data) <- rn
    return(feature_data)
  }
  composite <- grepl(":", rn, fixed = TRUE)
  if (all(composite)) {
    return(data.frame(
      gene_id = sub(":.*$", "", rn),
      isoform_id = sub("^[^:]*:", "", rn),
      stringsAsFactors = FALSE, row.names = rn
    ))
  }
  if (any(composite)) {
    stop("Count row names must use one consistent identifier style.",
         call. = FALSE)
  }
  normalized <- .normalize_transcript_ids(rn)
  red <- match(normalized, map$redundant_ENST)
  kept <- match(normalized, map$kept_ENST)
  gene <- rep(NA_character_, length(rn))
  gene[!is.na(red)] <- map$redundant_ENSG[red[!is.na(red)]]
  gene[is.na(red) & !is.na(kept)] <- map$kept_ENSG[kept[is.na(red) & !is.na(kept)]]
  data.frame(
    gene_id = gene, isoform_id = rn,
    stringsAsFactors = FALSE, row.names = rn
  )
}

.aggregate_duplicate_matrix <- function(counts, feature_data, map) {
  original_class <- if (is.data.frame(counts)) "data.frame" else "matrix"
  matrix_counts <- as.matrix(counts)
  if (!is.numeric(matrix_counts) || anyNA(matrix_counts) ||
      any(!is.finite(matrix_counts)) || any(matrix_counts < 0) ||
      any(matrix_counts != floor(matrix_counts))) {
    stop("Counts must be finite nonnegative integers without missing values.",
         call. = FALSE)
  }
  input_integer <- is.integer(matrix_counts)
  feature_data <- .feature_data_from_counts(matrix_counts, feature_data, map)
  ids <- as.character(feature_data$isoform_id)
  genes <- as.character(feature_data$gene_id)
  normalized <- .normalize_transcript_ids(ids)
  red <- match(normalized, map$redundant_ENST)
  kept <- match(normalized, map$kept_ENST)
  mapped <- !is.na(red) | !is.na(kept)
  expected_gene <- rep(NA_character_, length(ids))
  expected_gene[!is.na(red)] <- map$redundant_ENSG[red[!is.na(red)]]
  expected_gene[is.na(red) & !is.na(kept)] <-
    map$kept_ENSG[kept[is.na(red) & !is.na(kept)]]
  conflict <- mapped & !is.na(genes) & nzchar(genes) & genes != expected_gene
  if (any(conflict)) {
    stop(
      "Count feature gene assignment disagrees with the duplicate map for: ",
      paste(unique(ids[conflict]), collapse = ", "), ".", call. = FALSE
    )
  }
  resolved_id <- ids
  resolved_id[!is.na(red)] <- map$kept_ENST[red[!is.na(red)]]
  resolved_id[is.na(red) & !is.na(kept)] <- normalized[is.na(red) & !is.na(kept)]
  resolved_gene <- genes
  resolved_gene[mapped] <- expected_gene[mapped]
  safe_gene <- ifelse(is.na(resolved_gene) | !nzchar(resolved_gene),
                      paste0(".unknown.", seq_along(ids)), resolved_gene)
  group_key <- paste(safe_gene, resolved_id, sep = "\r")
  group_levels <- unique(group_key)
  group_factor <- factor(group_key, levels = group_levels)
  aggregated <- rowsum(
    matrix_counts * 1, group = group_factor, reorder = FALSE, na.rm = FALSE
  )
  rownames(aggregated) <- NULL
  groups <- split(seq_along(group_key), group_factor)
  representative <- vapply(groups, function(index) {
    target <- resolved_id[index[[1L]]]
    existing <- index[normalized[index] == target & is.na(red[index])]
    if (length(existing)) existing[[1L]] else index[[1L]]
  }, integer(1L))
  out_features <- feature_data[representative, , drop = FALSE]
  out_features$gene_id <- resolved_gene[representative]
  out_features$isoform_id <- resolved_id[representative]
  composite_style <- all(grepl(":", rownames(matrix_counts), fixed = TRUE))
  output_names <- if (composite_style) {
    paste(out_features$gene_id, out_features$isoform_id, sep = ":")
  } else {
    out_features$isoform_id
  }
  if (anyNA(output_names) || anyDuplicated(output_names)) {
    stop("Aggregation produced invalid or duplicate output feature names.",
         call. = FALSE)
  }
  rownames(aggregated) <- output_names
  rownames(out_features) <- output_names
  if (input_integer) {
    if (any(aggregated > .Machine$integer.max)) {
      stop("Aggregated counts exceed R's integer range.", call. = FALSE)
    }
    storage.mode(aggregated) <- "integer"
  }
  output_object <- if (identical(original_class, "data.frame")) {
    as.data.frame(aggregated, check.names = FALSE)
  } else aggregated
  resolution <- rep("unchanged", length(ids))
  resolution[!is.na(kept) & is.na(red) & ids != normalized] <-
    "version_normalized"
  resolution[!is.na(red)] <- "duplicate_aggregated"
  output_index <- match(group_key, group_levels)
  audit <- data.frame(
    input_row = seq_along(ids),
    input_feature = rownames(matrix_counts),
    submitted_gene = genes,
    submitted_id = ids,
    normalized_id = normalized,
    resolved_gene = resolved_gene,
    resolved_id = resolved_id,
    resolution = resolution,
    output_feature = output_names[output_index],
    stringsAsFactors = FALSE
  )
  summary <- data.frame(
    n_input_rows = nrow(matrix_counts),
    n_output_rows = nrow(aggregated),
    n_rows_removed = nrow(matrix_counts) - nrow(aggregated),
    n_redundant_rows_present = sum(!is.na(red)),
    n_representative_groups_affected = length(unique(group_key[!is.na(red)])),
    stringsAsFactors = FALSE
  )
  list(
    object = output_object,
    feature_data = out_features,
    audit = audit,
    summary = summary
  )
}

#' Aggregate counts for protein-identical isoforms
#'
#' Sums raw count rows into the kept transcript representative for each
#' protein-identical group. This transformation must be performed before
#' differential expression or transcript-usage model fitting.
#'
#' @param x A numeric matrix/data frame of raw counts, or a single-assay
#'   \code{SummarizedExperiment}.
#' @param feature_data Optional data frame aligned to matrix rows, containing
#'   at least \code{gene_id} and \code{isoform_id}. When omitted, identifiers
#'   are parsed from ENST or \code{gene:ENST} row names.
#' @param assay Assay name used for a \code{SummarizedExperiment}; defaults to
#'   \code{"counts"}.
#' @param duplicate_map Optional validated custom map.
#'
#' @return An object of class \code{DuplicateCountAggregation}: a list with
#'   \code{object}, \code{feature_data}, \code{audit}, and \code{summary}.
#' @export
aggregate_duplicate_counts <- function(
    x, feature_data = NULL, assay = "counts", duplicate_map = NULL) {
  map <- .validate_duplicate_map(
    if (is.null(duplicate_map)) isonet_duplicate_map() else duplicate_map
  )
  is_se <- inherits(x, "SummarizedExperiment")
  if (is_se) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
        !requireNamespace("S4Vectors", quietly = TRUE)) {
      stop(
        "SummarizedExperiment and S4Vectors are required for this input.",
        call. = FALSE
      )
    }
    assay_names <- SummarizedExperiment::assayNames(x)
    if (length(assay_names) != 1L) {
      stop(
        "SummarizedExperiment input must contain exactly one raw-count assay. ",
        "For multi-assay objects, extract the count matrix and rowData first.",
        call. = FALSE
      )
    }
    if (!is.character(assay) || length(assay) != 1L ||
        !assay %in% assay_names) {
      stop("`assay` is not present in the SummarizedExperiment.", call. = FALSE)
    }
    se_feature_data <- as.data.frame(SummarizedExperiment::rowData(x))
    result <- .aggregate_duplicate_matrix(
      SummarizedExperiment::assay(x, assay), se_feature_data, map
    )
    rebuilt <- SummarizedExperiment::SummarizedExperiment(
      assays = stats::setNames(list(result$object), assay),
      rowData = S4Vectors::DataFrame(
        result$feature_data, row.names = rownames(result$feature_data)
      ),
      colData = SummarizedExperiment::colData(x),
      metadata = S4Vectors::metadata(x)
    )
    S4Vectors::metadata(rebuilt)$isonet_duplicate_aggregation <- list(
      summary = result$summary,
      audit = result$audit,
      map_provenance = attr(map, "provenance")
    )
    result$object <- rebuilt
  } else {
    if (!is.matrix(x) && !is.data.frame(x)) {
      stop("`x` must be a matrix, data frame, or SummarizedExperiment.",
           call. = FALSE)
    }
    result <- .aggregate_duplicate_matrix(x, feature_data, map)
  }
  result$map_provenance <- attr(map, "provenance")
  class(result) <- c("DuplicateCountAggregation", "list")
  result
}

#' @export
print.DuplicateCountAggregation <- function(x, ...) {
  cat("Protein-identical count aggregation\n")
  cat("  Input rows :", x$summary$n_input_rows, "\n")
  cat("  Output rows:", x$summary$n_output_rows, "\n")
  cat("  Removed    :", x$summary$n_rows_removed, "\n")
  invisible(x)
}
