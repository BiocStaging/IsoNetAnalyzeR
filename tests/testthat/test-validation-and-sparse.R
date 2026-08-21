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
