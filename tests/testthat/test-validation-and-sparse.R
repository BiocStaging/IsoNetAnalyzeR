# Expected values follow from the design of the hand-built fixtures in
# helper-fixtures.R.
test_that("differential pulldown rejects malformed public inputs", {
  fixture <- make_isolated_seed_fixture()

  expect_error(
    run_differential_pulldown(
      fixture$obj, "", fixture$input, order = 1,
      confirm = FALSE, verbose = FALSE
    ),
    "label"
  )
  expect_error(
    run_differential_pulldown(
      fixture$obj, "bad_order", fixture$input, order = 0,
      confirm = FALSE, verbose = FALSE
    ),
    "positive integer"
  )
  expect_error(
    run_differential_pulldown(
      fixture$obj, "bad_cutoff", fixture$input, order = 1, cutoff = 2,
      confirm = FALSE, verbose = FALSE
    ),
    "\\[0, 1\\]"
  )
  expect_error(
    run_differential_pulldown(
      fixture$obj, "missing_column",
      fixture$input[, c("gene_id", "isoform_id")],
      order = 1, confirm = FALSE, verbose = FALSE
    ),
    "missing required column"
  )

  bad_direction <- fixture$input
  bad_direction$direction[[1]] <- "sideways"
  expect_error(
    run_differential_pulldown(
      fixture$obj, "bad_direction", bad_direction,
      order = 1, confirm = FALSE, verbose = FALSE
    ),
    "Unrecognised direction"
  )
})

test_that("decay defaults filter seeds and neighbours with independent opt-ins", {
  fixture <- make_isolated_seed_fixture()
  graph <- fixture$obj@network
  igraph::V(graph)$transcript_biotype <- c(
    "protein_coding", "protein_coding", "non_stop_decay",
    "nonsense_mediated_decay", "protein_coding", "protein_coding"
  )
  make <- function(nmd = FALSE, nsd = FALSE) IsoNetAnalyzeR(
    "decay", graph, metadata = list(network_cutoff = 0),
    include_nmd = nmd, include_nsd = nsd
  )
  # tx_neighbour is NSD and tx_isolated_up NMD; both edges touch tx_neighbour
  obj <- make()
  expect_false(any(c("tx_neighbour", "tx_isolated_up") %in% igraph::V(obj@network)$name))
  expect_equal(igraph::ecount(obj@network), 0L)
  expect_true("tx_isolated_up" %in% igraph::V(make(TRUE, FALSE)@network)$name)
  expect_false("tx_neighbour" %in% igraph::V(make(TRUE, FALSE)@network)$name)
  expect_true("tx_neighbour" %in% igraph::V(make(FALSE, TRUE)@network)$name)
  expect_setequal(igraph::V(make(TRUE, TRUE)@network)$name, igraph::V(graph)$name)
  input <- fixture$input[fixture$input$gene_id %in% c("G1", "G3"), ]
  result <- run_differential_pulldown(obj, "filtered", input, order = 1,
                                    confirm = FALSE, verbose = FALSE)
  expect_false("tx_isolated_up" %in% get_analysis(result, "filtered")@input_table$isoform_id)
  expect_equal(result@metadata$analysis_biotype_exclusions$filtered$isoform_id,
               "tx_isolated_up")
  expect_error(run_differential_pulldown(
    obj, "all_excluded", input[input$isoform_id == "tx_isolated_up", ],
    order = 1, confirm = FALSE, verbose = FALSE
  ), "No seeds remain")
  expect_error(make(NA, FALSE), "include_nmd")
  expect_error(make(FALSE, "yes"), "include_nsd")
})

test_that("custom networks require transcript annotation unless explicitly unfiltered", {
  graph <- make_isolated_seed_fixture()$obj@network
  graph <- igraph::delete_vertex_attr(graph, "transcript_biotype")
  expect_error(IsoNetAnalyzeR("unknown", graph), "unknown transcript biotypes")
  expect_s4_class(IsoNetAnalyzeR("unfiltered", graph,
                                include_nmd = TRUE, include_nsd = TRUE), "IsoNetAnalyzeR")
  ann <- data.frame(isoform_id = igraph::V(graph)$name,
                    transcript_biotype = "protein_coding")
  expect_s4_class(IsoNetAnalyzeR("annotated", graph, transcript_biotypes = ann),
                  "IsoNetAnalyzeR")
  conflict <- rbind(ann, transform(ann[1, ], transcript_biotype = "non_stop_decay"))
  expect_error(IsoNetAnalyzeR("conflicting", graph, transcript_biotypes = conflict),
               "Conflicting transcript biotypes")
})

test_that("decay-excluded aliases cannot bypass eligibility by duplicate remapping", {
  fixture <- make_isolated_seed_fixture()
  ann <- data.frame(isoform_id = "excluded_alias",
                    transcript_biotype = "nonsense_mediated_decay")
  obj <- IsoNetAnalyzeR("alias", fixture$obj@network, transcript_biotypes = ann)
  map <- data.frame(kept_ENST = "tx_up", redundant_ENST = "excluded_alias",
                    kept_ENSG = "G1", redundant_ENSG = "G1",
                    kept_ENSP = "p_up", redundant_ENSP = "p_alias",
                    network_applicable = TRUE)
  resolved <- resolve_duplicate_isoforms("excluded_alias.7", obj@network,
                                         duplicate_map = map)
  expect_equal(as.character(resolved), "excluded_alias.7")
  expect_equal(attr(resolved, "dedup_audit")$resolution,
               "excluded_transcript_biotype")
  expect_error(virtual_pulldown(obj, "excluded_alias.7", order = 1, cutoff = 0),
               "None of the seed vertices")
})

test_that("retained support requires a significant mapped seed and both directions", {
  fixture <- make_isolated_seed_fixture()
  input <- fixture$input
  # G1 keeps a significant mapped seed (q 0.01) in both directions; G3's
  # mapped seeds are not significant; G_missing has no mapped seed. Seeds are
  # significant when q < seed_fdr, so q = 0.05 retains no gene.
  input$q <- c(.01, NA, .2, .3, .001, .001)
  # A significant but unmapped G3 transcript cannot support its mapped pair.
  input <- rbind(input, data.frame(gene_id = "G3", isoform_id = "missing_driver",
                                 direction = "up", q = .001))
  obj <- run_differential_pulldown(fixture$obj, "supported", input,
    order = 1, cutoff = 0, confirm = FALSE, verbose = FALSE, seed_padj_col = "q")
  expect_named(get_analysis(obj, "supported")@gene_results, "G1")
  expect_setequal(get_analysis(obj, "supported")@input_table$isoform_id,
                  c("tx_up", "tx_down"))
  audit <- obj@metadata$analysis_seed_support$supported$gene_audit
  expect_equal(audit$reason[audit$gene_id == "G3"], "no_significant_retained_seed")
  expect_equal(audit$reason[audit$gene_id == "G_missing"], "not_bidirectional_after_mapping")
  expect_error(run_differential_pulldown(fixture$obj, "no_column", fixture$input,
    order = 1, confirm = FALSE, seed_padj_col = "q"), "must name")
  input$q <- .05
  expect_error(run_differential_pulldown(fixture$obj, "boundary", input,
    order = 1, confirm = FALSE, seed_padj_col = "q"), "No genes retain")
  input$q <- -1
  expect_error(run_differential_pulldown(fixture$obj, "bad_p", input,
    order = 1, confirm = FALSE, seed_padj_col = "q"), "numeric in")
  expect_error(run_differential_pulldown(fixture$obj, "bad_alpha", input,
    order = 1, confirm = FALSE, seed_padj_col = "q", seed_fdr = 0), "seed_fdr")
})

test_that("a decay-excluded significant seed cannot support a coding pair", {
  fixture <- make_isolated_seed_fixture()
  graph <- igraph::add_vertices(fixture$obj@network, 1,
    attr = list(name = "decay_driver", gene = "G3", transcript_biotype = "nonsense_mediated_decay"))
  obj <- IsoNetAnalyzeR("decay_support", graph, metadata = list(network_cutoff = 0))
  input <- fixture$input[1:4, ]
  input$q <- c(.01, .2, .2, .2)
  input <- rbind(input, data.frame(gene_id = "G3", isoform_id = "decay_driver",
                                 direction = "up", q = .001))
  obj <- run_differential_pulldown(obj, "checked", input, order = 1, cutoff = 0,
    confirm = FALSE, verbose = FALSE, seed_padj_col = "q")
  expect_named(get_analysis(obj, "checked")@gene_results, "G1")
  expect_equal(obj@metadata$analysis_biotype_exclusions$checked$isoform_id, "decay_driver")
})

test_that("lookup and component entry points reject invalid state", {
  fixture <- make_isolated_seed_fixture()

  expect_error(
    lookup_virtual_pulldown(
      fixture$obj, "missing", "tx_not_in_network",
      order = 1, cutoff = 0
    ),
    "not found in the network"
  )
  expect_error(
    lookup_virtual_pulldown(
      fixture$obj, "bad_order", "tx_up",
      order = 0, cutoff = 0
    ),
    "positive integer"
  )
  expect_error(
    run_component_analysis(fixture$obj, "missing", "components"),
    "No analysis labelled"
  )
  expect_error(
    run_network_diagnostics(fixture$obj, sticky_cutoff = 0),
    "positive numeric"
  )
})

test_that("fully unmapped and single-direction inputs remain typed and reportable", {
  fixture <- make_isolated_seed_fixture()
  unmapped <- data.frame(
    gene_id = c("G_missing", "G_missing"),
    isoform_id = c("tx_absent_up", "tx_absent_down"),
    direction = c("up", "down"),
    stringsAsFactors = FALSE
  )
  suppressWarnings(
    obj <- run_differential_pulldown(
      fixture$obj, "unmapped", unmapped, order = 1, cutoff = 0,
      confirm = FALSE, verbose = FALSE
    )
  )
  result <- get_gene_result(obj, "unmapped", "G_missing")
  expect_null(result@up_graph)
  expect_null(result@down_graph)
  expect_named(
    result@up_stats,
    c(
      "isoform", "is_seed", "degree", "clustering", "distance_from_seed",
      "mean_pi", "pi_entropy"
    )
  )
  expect_equal(nrow(result@up_stats), 0L)

  obj <- run_component_analysis(
    obj, "unmapped", "empty_components",
    direction = "both", level = "gene"
  )
  overview <- get_component_overview(obj, "unmapped", "empty_components")
  expect_s4_class(overview, "DataFrame")
  expect_equal(nrow(overview), 0L)
  expect_true(all(c(
    "component_id", "n_nodes", "n_edges", "n_genes_up", "n_genes_down",
    "mixed"
  ) %in% names(overview)))

  one_direction <- data.frame(
    gene_id = "G1",
    isoform_id = "tx_up",
    direction = "up",
    stringsAsFactors = FALSE
  )
  obj <- run_differential_pulldown(
    fixture$obj, "up_only", one_direction, order = 1, cutoff = 0,
    confirm = FALSE, verbose = FALSE
  )
  result <- get_gene_result(obj, "up_only", "G1")
  expect_s4_class(result@up_stats, "DataFrame")
  expect_null(result@down_graph)
  expect_equal(nrow(result@down_stats), 0L)
})
