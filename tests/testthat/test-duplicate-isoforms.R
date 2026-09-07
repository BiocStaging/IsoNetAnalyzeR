make_duplicate_fixture <- function() {
  map <- isonet_duplicate_map()
  row <- map[which(map$network_applicable)[1L], , drop = FALSE]
  vertices <- data.frame(
    name = c(row$kept_ENST, "same_gene_down", "other_up", "other_down"),
    gene = c(row$kept_ENSG, row$kept_ENSG, "G_OTHER", "G_OTHER"),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    from = c(row$kept_ENST, "same_gene_down", "other_up"),
    to = c("other_up", "other_up", "other_down"),
    interact = 1L,
    pi = c(0.9, 0.8, 0.7),
    pi_score = c(0.9, 0.8, 0.7),
    stringsAsFactors = FALSE
  )
  vertices$transcript_biotype <- "protein_coding"
  network <- igraph::graph_from_data_frame(
    edges, directed = FALSE, vertices = vertices
  )
  obj <- IsoNetAnalyzeR(
    "duplicate-fixture", network,
    metadata = list(network_cutoff = 0, pi_threshold = 0.5)
  )
  list(map = map, row = row, network = network, obj = obj)
}

test_that("packaged duplicate map has reviewed dimensions and provenance", {
  map <- isonet_duplicate_map(refresh = TRUE)
  expect_equal(nrow(map), 16630L)
  expect_equal(length(unique(map$kept_ENST)), 9760L)
  expect_equal(length(unique(map$redundant_ENST)), 16630L)
  expect_equal(sum(map$network_applicable), 10592L)
  expect_equal(sum(!map$network_applicable), 6038L)
  expect_false(any(map$redundant_ENST %in% map$kept_ENST))
  provenance <- attr(map, "provenance")
  expect_identical(
    provenance$source_sha256,
    "cfe6fb57ed06c7231b9c77d7246b148be3cfdb1ac6aad4f9895cc422efaa36f5"
  )
})

test_that("resolver handles duplicate, versioned, and incompatible IDs", {
  fixture <- make_duplicate_fixture()
  row <- fixture$row
  ids <- c(
    row$redundant_ENST,
    paste0(row$kept_ENST, ".7"),
    "not_in_map"
  )
  resolved <- resolve_duplicate_isoforms(ids, fixture$network)
  expect_identical(
    as.vector(resolved), c(row$kept_ENST, row$kept_ENST, "not_in_map")
  )
  expect_identical(
    attr(resolved, "dedup_audit")$resolution,
    c("duplicate_remapped", "version_normalized", "unchanged_unmapped")
  )

  incompatible <- igraph::add_vertices(
    fixture$network, 1L,
    attr = list(name = row$redundant_ENST, gene = row$kept_ENSG)
  )
  existing <- resolve_duplicate_isoforms(row$redundant_ENST, incompatible)
  expect_identical(as.vector(existing), row$redundant_ENST)
  expect_identical(attr(existing, "dedup_audit")$resolution, "direct")
})

test_that("direct pulldown and lookup resolve duplicate seeds and retain audits", {
  fixture <- make_duplicate_fixture()
  row <- fixture$row
  graph <- virtual_pulldown(
    fixture$obj, row$redundant_ENST, order = 1L, cutoff = 0
  )
  expect_true(row$kept_ENST %in% igraph::V(graph)$name)
  audit <- igraph::graph_attr(graph, "duplicate_resolution_audit")
  expect_identical(audit$resolution, "duplicate_remapped")

  obj <- lookup_virtual_pulldown(
    fixture$obj, "duplicate", row$redundant_ENST,
    order = 1L, cutoff = 0
  )
  lookup <- get_lookup(obj, "duplicate")
  expect_identical(
    lookup@params$duplicate_resolution$resolution,
    "duplicate_remapped"
  )
  expect_error(
    lookup_virtual_pulldown(
      fixture$obj, "missing", "not_in_map_or_network",
      order = 1L, cutoff = 0
    ),
    "not found in the network"
  )
})

test_that("differential resolution preserves submitted and resolved counts", {
  fixture <- make_duplicate_fixture()
  row <- fixture$row
  input <- data.frame(
    gene_id = c(row$kept_ENSG, row$kept_ENSG, row$kept_ENSG),
    isoform_id = c(row$kept_ENST, row$redundant_ENST, "same_gene_down"),
    direction = c("up", "up", "down"),
    extra = c("kept", "redundant", "down"),
    stringsAsFactors = FALSE
  )
  obj <- suppressMessages(run_differential_pulldown(
    fixture$obj, "same", input, order = 1L, cutoff = 0,
    confirm = FALSE, verbose = FALSE
  ))
  analysis <- get_analysis(obj, "same")
  expect_equal(nrow(analysis@input_table), 2L)
  expect_true("extra" %in% names(analysis@input_table))
  overview <- analysis@cross_gene_summary
  expect_equal(overview$n_up_input, 2L)
  expect_equal(overview$n_down_input, 1L)
  expect_equal(overview$n_up_resolved, 1L)
  expect_equal(overview$n_down_resolved, 1L)
  stored <- obj@metadata$analysis_duplicate_resolution$same
  expect_equal(nrow(stored$submitted_table), 3L)
  expect_equal(sum(stored$audit$action == "collapsed_same_direction"), 1L)
  coverage <- obj@metadata$analysis_network_coverage$same
  expect_equal(coverage$overall$n_input, 3L)
  expect_equal(coverage$overall$n_resolved_representatives, 2L)
  expect_equal(coverage$overall$n_remapped_from_duplicate, 1L)

  component_gene <- igraph::make_empty_graph(1L)
  component_gene <- igraph::set_vertex_attr(
    component_gene, "name", value = row$kept_ENSG
  )
  score <- IsoNetAnalyzeR:::.compute_dominance_scores(
    component_gene, analysis@input_table
  )
  expect_equal(unname(score), 0)
})

test_that("opposite protein-identical directions drop the whole gene", {
  fixture <- make_duplicate_fixture()
  row <- fixture$row
  input <- data.frame(
    gene_id = c(row$kept_ENSG, row$kept_ENSG, "G_OTHER", "G_OTHER"),
    isoform_id = c(
      row$kept_ENST, row$redundant_ENST, "other_up", "other_down"
    ),
    direction = c("up", "down", "up", "down"),
    stringsAsFactors = FALSE
  )
  expect_warning(
    obj <- suppressMessages(run_differential_pulldown(
      fixture$obj, "opposite", input, order = 1L, cutoff = 0,
      confirm = FALSE, verbose = FALSE
    )),
    "opposite directions"
  )
  analysis <- get_analysis(obj, "opposite")
  expect_identical(unique(analysis@input_table$gene_id), "G_OTHER")
  expect_equal(
    obj@metadata$analysis_duplicate_resolution$opposite$summary$
      n_opposite_direction_genes_dropped,
    1L
  )

  only_collision <- input[seq_len(2L), , drop = FALSE]
  expect_error(
    suppressWarnings(IsoNetAnalyzeR:::.validate_diff_table(
      only_collision, fixture$network, resolve_duplicates = TRUE
    )),
    "No analyzable rows remain"
  )
})

test_that("resolve_duplicates FALSE retains 0.4.0 identifier behavior", {
  fixture <- make_duplicate_fixture()
  row <- fixture$row
  input <- data.frame(
    gene_id = row$kept_ENSG,
    isoform_id = row$redundant_ENST,
    direction = "up",
    stringsAsFactors = FALSE
  )
  cleaned <- IsoNetAnalyzeR:::.validate_diff_table(
    input, fixture$network, resolve_duplicates = FALSE
  )
  expect_identical(cleaned, input)
  expect_null(attr(cleaned, "dedup_resolution"))
})

test_that("base count aggregation preserves arithmetic and metadata", {
  fixture <- make_duplicate_fixture()
  row <- fixture$row
  ids <- c(row$kept_ENST, row$redundant_ENST, "unrelated")
  genes <- c(row$kept_ENSG, row$kept_ENSG, "G_OTHER")
  rn <- paste(genes, ids, sep = ":")
  counts <- matrix(
    c(1L, 2L, 7L, 3L, 4L, 8L), nrow = 3L,
    dimnames = list(rn, c("sample_a", "sample_b"))
  )
  feature_data <- data.frame(
    gene_id = genes, isoform_id = ids,
    symbol = c("kept metadata", "redundant metadata", "other"),
    row.names = rn, stringsAsFactors = FALSE
  )
  result <- aggregate_duplicate_counts(counts, feature_data)
  expect_s3_class(result, "DuplicateCountAggregation")
  expect_equal(nrow(result$object), 2L)
  expect_equal(colSums(result$object), colSums(counts))
  expect_equal(result$summary$n_rows_removed, 1L)
  target <- paste(row$kept_ENSG, row$kept_ENST, sep = ":")
  expect_equal(unname(result$object[target, ]), c(3L, 7L))
  expect_identical(result$feature_data[target, "symbol"], "kept metadata")
})

test_that("aggregation creates an absent kept row and supports one-assay SE", {
  map <- data.frame(
    kept_ENSG = "G", redundant_ENSG = "G",
    kept_ENSP = "P_KEEP", kept_ENST = "tx_keep",
    redundant_ENSP = "P_RED", redundant_ENST = "tx_red",
    network_applicable = FALSE, stringsAsFactors = FALSE
  )
  counts <- matrix(
    c(2L, 5L), nrow = 1L,
    dimnames = list("G:tx_red", c("a", "b"))
  )
  feature_data <- data.frame(
    gene_id = "G", isoform_id = "tx_red", note = "source",
    row.names = "G:tx_red", stringsAsFactors = FALSE
  )
  result <- aggregate_duplicate_counts(
    counts, feature_data, duplicate_map = map
  )
  expect_identical(rownames(result$object), "G:tx_keep")
  expect_identical(result$feature_data$isoform_id, "tx_keep")

  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(feature_data)
  )
  se_result <- aggregate_duplicate_counts(se, duplicate_map = map)
  expect_s4_class(se_result$object, "SummarizedExperiment")
  expect_identical(rownames(se_result$object), "G:tx_keep")
  expect_equal(
    S4Vectors::metadata(se_result$object)$isonet_duplicate_aggregation$
      summary$n_rows_removed,
    0L
  )
})
