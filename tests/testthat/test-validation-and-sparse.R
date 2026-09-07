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
  expect_s3_class(overview, "data.frame")
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
  expect_s3_class(result@up_stats, "data.frame")
  expect_null(result@down_graph)
  expect_equal(nrow(result@down_stats), 0L)
})
