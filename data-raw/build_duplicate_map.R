#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_full[grepl("^--file=", args_full)])
source_root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(file_arg[[1L]]), ".."), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required to build the duplicate map.",
         call. = FALSE)
  }
}
require_package("digest")
require_package("igraph")
require_package("SummarizedExperiment")

sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
unique_count_by <- function(value, group) {
  vapply(split(value, group), function(x) length(unique(x)), integer(1L))
}

raw_path <- file.path(
  source_root, "data-raw",
  "duplicates_long_kept_redundant_enst_ensp_ensg.csv"
)
se_path <- normalizePath(
  file.path(source_root, "..", "figure3_work", "data", "interim",
            "GSE141432_transcript_se.rds"),
  mustWork = TRUE
)
sysdata_path <- file.path(source_root, "R", "sysdata.rda")
out_dir <- file.path(source_root, "inst", "extdata")
out_path <- file.path(out_dir, "protein_identical_isoform_map.csv")
provenance_path <- file.path(
  out_dir, "protein_identical_isoform_map_provenance.csv"
)

expected_source_sha256 <-
  "cfe6fb57ed06c7231b9c77d7246b148be3cfdb1ac6aad4f9895cc422efaa36f5"
expected_se_sha256 <-
  "8ab52523b1564b061c9424e5189dcf6014f0f4e0161d1eb1e6f87c855dbe7b35"
expected_sysdata_sha256 <-
  "3130209f090b9b83bbf2ebf328f29a0ad9a8fcd89922a277d0966325bc163ebb"

assert(file.exists(raw_path), paste("Missing source map:", raw_path))
assert(identical(sha256(raw_path), expected_source_sha256),
       "Source duplicate-map SHA-256 does not match the reviewed attachment.")
assert(identical(sha256(se_path), expected_se_sha256),
       "GSE141432 SE SHA-256 does not match the pinned annotation evidence.")
assert(identical(sha256(sysdata_path), expected_sysdata_sha256),
       "R/sysdata.rda changed; refusing to build against an unreviewed network.")

map <- utils::read.csv(
  raw_path, colClasses = "character", check.names = FALSE,
  stringsAsFactors = FALSE
)
required <- c(
  "kept_ENSG", "kept_ENSP", "kept_ENST",
  "redundant_ENSP", "redundant_ENST"
)
assert(identical(names(map), required),
       "Duplicate-map source schema differs from the five reviewed columns.")
assert(nrow(map) == 16630L, "Duplicate-map source must contain 16,630 rows.")
assert(!anyNA(map) && all(nzchar(trimws(as.matrix(map)))),
       "Duplicate-map source contains missing or blank values.")
assert(!anyDuplicated(map), "Duplicate-map source contains duplicate rows.")
assert(all(grepl("^ENSG[0-9]+$", map$kept_ENSG)),
       "Invalid kept ENSG identifier.")
assert(all(grepl("^ENSP[0-9]+$", c(map$kept_ENSP, map$redundant_ENSP))),
       "Invalid ENSP identifier.")
assert(all(grepl("^ENST[0-9]+$", c(map$kept_ENST, map$redundant_ENST))),
       "Invalid ENST identifier.")
assert(!length(intersect(map$kept_ENST, map$redundant_ENST)),
       "An ENST occurs as both kept and redundant (transitive chain).")
assert(all(unique_count_by(map$kept_ENST, map$redundant_ENST) == 1L),
       "A redundant ENST maps to multiple kept ENSTs.")
assert(all(unique_count_by(map$kept_ENST, map$kept_ENSP) == 1L) &&
         all(unique_count_by(map$kept_ENSG, map$kept_ENSP) == 1L) &&
         all(unique_count_by(map$kept_ENSP, map$kept_ENST) == 1L),
       "A kept protein/transcript group is internally inconsistent.")
assert(all(unique_count_by(map$redundant_ENST, map$redundant_ENSP) == 1L),
       "A redundant protein identifier refers to multiple transcripts.")

network_env <- new.env(parent = emptyenv())
load(sysdata_path, envir = network_env)
network <- network_env$.network_data
vertices <- igraph::V(network)$name
vertex_gene <- stats::setNames(as.character(igraph::V(network)$gene), vertices)
assert(!any(map$redundant_ENST %in% vertices),
       "A purported redundant ENST is present in the shipped network.")
network_applicable <- map$kept_ENST %in% vertices
assert(sum(network_applicable) == 10592L,
       "Unexpected number of network-applicable duplicate-map rows.")
assert(sum(!network_applicable) == 6038L,
       "Unexpected number of transcriptome-only duplicate-map rows.")
assert(all(
  map$kept_ENSG[network_applicable] ==
    unname(vertex_gene[map$kept_ENST[network_applicable]])
), "A network-applicable kept transcript has a conflicting gene annotation.")

se <- readRDS(se_path)
row_data <- as.data.frame(SummarizedExperiment::rowData(se))
assert(all(c("gene_id", "isoform_id") %in% names(row_data)),
       "Pinned SE lacks gene_id/isoform_id row metadata.")
se_gene <- stats::setNames(
  as.character(row_data$gene_id), as.character(row_data$isoform_id)
)
assert(all(map$kept_ENST %in% names(se_gene)) &&
         all(map$redundant_ENST %in% names(se_gene)),
       "Pinned SE does not contain every duplicate-map transcript.")
assert(all(unname(se_gene[map$kept_ENST]) == map$kept_ENSG) &&
         all(unname(se_gene[map$redundant_ENST]) == map$kept_ENSG),
       "Pinned SE contradicts the single-gene duplicate-map assumption.")

shipped <- data.frame(
  kept_ENSG = map$kept_ENSG,
  redundant_ENSG = map$kept_ENSG,
  kept_ENSP = map$kept_ENSP,
  kept_ENST = map$kept_ENST,
  redundant_ENSP = map$redundant_ENSP,
  redundant_ENST = map$redundant_ENST,
  network_applicable = network_applicable,
  stringsAsFactors = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(shipped, out_path, row.names = FALSE, quote = TRUE,
                 na = "")

provenance <- data.frame(
  key = c(
    "source_filename", "source_sha256", "artifact_sha256", "build_date",
    "row_count", "unique_kept_enst", "unique_redundant_enst",
    "network_applicable_rows", "transcriptome_only_rows",
    "redundant_ensg_derivation", "annotation_evidence_path",
    "annotation_evidence_sha256", "network_sysdata_sha256"
  ),
  value = c(
    basename(raw_path), expected_source_sha256, sha256(out_path),
    as.character(Sys.Date()), nrow(shipped), length(unique(shipped$kept_ENST)),
    length(unique(shipped$redundant_ENST)), sum(network_applicable),
    sum(!network_applicable),
    "derived as kept_ENSG; zero mismatches in pinned GSE141432 SE",
    "figure3_work/data/interim/GSE141432_transcript_se.rds",
    expected_se_sha256, expected_sysdata_sha256
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(provenance, provenance_path, row.names = FALSE, quote = TRUE)

message(
  "Built ", out_path, ": ", nrow(shipped), " rows; ",
  sum(network_applicable), " network-applicable; artifact SHA-256 ",
  sha256(out_path)
)
