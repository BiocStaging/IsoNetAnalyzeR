test_that("DataFrame input gives the same pulldowns as data.frame input", {
  fx <- make_multi_gene_fixture()
  obj <- suppressMessages(run_differential_pulldown(
    fx$obj, "from_DataFrame", S4Vectors::DataFrame(fx$seeds), order = 1,
    confirm = FALSE, verbose = FALSE
  ))
  expect_equal(get_cross_gene_summary(obj, "from_DataFrame"),
               get_cross_gene_summary(obj, "a"))
  expect_equal(get_isoform_stats(obj, "from_DataFrame", "GB"),
               get_isoform_stats(obj, "a", "GB"))
})

test_that("result tables are DataFrames", {
  fx <- make_multi_gene_fixture()
  obj <- fx$obj
  expect_s4_class(get_cross_gene_summary(obj, "a"), "DataFrame")
  expect_s4_class(get_isoform_stats(obj, "a", "GA"), "DataFrame")

  obj <- suppressMessages(run_component_analysis(obj, "a", label = "comp"))
  expect_s4_class(get_component_overview(obj, "a", "comp"), "DataFrame")
  expect_s4_class(get_component(obj, "a", "comp", 1), "DataFrame")

  withr::local_seed(1)
  obj <- suppressWarnings(suppressMessages(run_fgsea(
    obj, "a", label = "nodes", pathways = fx$pathways, minSize = 2
  )))
  expect_s4_class(get_gsea_table(obj, "a", "nodes"), "DataFrame")
  expect_s4_class(get_gsea_stats(obj, "a", "nodes"), "DataFrame")

  obj <- suppressMessages(run_component_ora(
    obj, "a", "comp", label = "ora", pathways = fx$pathways, minSize = 2
  ))
  expect_s4_class(get_ora_table(obj, "a", "comp", "ora"), "DataFrame")

  obj <- suppressMessages(lookup_virtual_pulldown(
    obj, "lk", seed = "GA_up", order = 1
  ))
  expect_s4_class(get_lookup_stats(obj, "lk"), "DataFrame")

  diag <- get_network_diagnostics(suppressMessages(
    run_network_diagnostics(obj)
  ))
  expect_s4_class(diag@network_stats, "DataFrame")
  expect_s4_class(diag@isoform_degree_table, "DataFrame")
})

test_that("duplicate helpers keep DataFrame input as DataFrame", {
  map <- isonet_duplicate_map()
  x <- S4Vectors::DataFrame(gene_id = map$kept_ENSG[1],
                            isoform_id = map$redundant_ENST[1])
  out <- resolve_duplicate_isoforms(x)
  expect_s4_class(out, "DataFrame")
  expect_identical(out$isoform_id, map$kept_ENST[1])
  expect_identical(S4Vectors::metadata(out)$dedup_audit$resolution,
                   "duplicate_remapped")

  # The redundant transcript's counts are added to its representative's
  counts <- S4Vectors::DataFrame(
    s1 = c(10L, 4L), s2 = c(2L, 1L),
    row.names = c(map$kept_ENST[1], map$redundant_ENST[1])
  )
  features <- S4Vectors::DataFrame(
    gene_id = c(map$kept_ENSG[1], map$redundant_ENSG[1]),
    isoform_id = rownames(counts)
  )
  res <- aggregate_duplicate_counts(counts, features)
  expect_s4_class(res$object, "DataFrame")
  expect_identical(rownames(res$object), map$kept_ENST[1])
  expect_identical(res$object$s1, 14L)
  expect_identical(res$object$s2, 3L)
})

test_that("transcript_biotypes may be a DataFrame", {
  fx <- make_isolated_seed_fixture()
  g <- get_network(fx$obj)
  annotation <- S4Vectors::DataFrame(
    isoform_id = c("tx_up", "tx_down", "tx_neighbour", "tx_isolated_up",
                   "tx_isolated_down", "tx_other"),
    transcript_biotype = c(rep("protein_coding", 5),
                           "nonsense_mediated_decay")
  )
  obj <- suppressMessages(IsoNetAnalyzeR("custom", network = g,
                                         metadata = list(network_cutoff = 0),
                                         transcript_biotypes = annotation))
  expect_false("tx_other" %in% igraph::V(get_network(obj))$name)
})
