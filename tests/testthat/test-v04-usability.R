test_that("coverage has three states, warns diagnostically, and persists", {
  fixture <- make_isolated_seed_fixture()
  obj <- suppressWarnings(run_differential_pulldown(
    fixture$obj, "coverage", fixture$input, order = 1, cutoff = 0,
    confirm = FALSE, verbose = FALSE
  ))

  coverage <- obj@metadata$analysis_network_coverage$coverage
  expect_named(coverage, c("overall", "by_direction", "genes",
                           "warning_threshold"))
  expect_equal(coverage$overall$n_unmapped, 2L)
  expect_equal(coverage$overall$n_mapped_no_predicted_edge, 2L)
  expect_equal(coverage$overall$n_mapped_with_predicted_edge, 2L)
  expect_type(coverage$overall$n_input, "integer")
  expect_type(coverage$overall$pct_mapped, "double")
  expect_setequal(coverage$by_direction$direction, c("up", "down"))

  low <- data.frame(
    gene_id = paste0("G", seq_len(11)),
    isoform_id = c("tx_up", paste0("absent_", seq_len(10))),
    direction = rep("up", 11),
    stringsAsFactors = FALSE
  )
  coverage_warnings <- character(0)
  withCallingHandlers(
    IsoNetAnalyzeR:::.report_network_coverage(
      low,
      igraph::V(fixture$obj@network)$name,
      unique(low$gene_id),
      confirm = FALSE,
      interacting_vertices = c("tx_up", "tx_down", "tx_neighbour")
    ),
    warning = function(w) {
      coverage_warnings <<- c(coverage_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("annotation-release|identifier-type", coverage_warnings)))
  expect_length(coverage_warnings, 2L)
})

test_that("input cleaning handles missing values, duplicates, and conflicts", {
  fixture <- make_isolated_seed_fixture()

  duplicated_input <- rbind(fixture$input[1, ], fixture$input[1, ])
  expect_warning(
    IsoNetAnalyzeR:::.validate_diff_table(duplicated_input, fixture$obj@network),
    "duplicate"
  )
  cleaned <- suppressWarnings(IsoNetAnalyzeR:::.validate_diff_table(
    duplicated_input, fixture$obj@network
  ))
  expect_equal(nrow(cleaned), 1L)

  missing_id <- rbind(
    fixture$input[1, ],
    data.frame(gene_id = NA_character_, isoform_id = "x", direction = "up")
  )
  expect_warning(
    IsoNetAnalyzeR:::.validate_diff_table(missing_id, fixture$obj@network),
    "missing or empty"
  )
  cleaned_missing <- suppressWarnings(IsoNetAnalyzeR:::.validate_diff_table(
    missing_id, fixture$obj@network
  ))
  expect_equal(nrow(cleaned_missing), 1L)

  missing_direction <- fixture$input[1, ]
  missing_direction$direction <- NA_character_
  expect_error(
    IsoNetAnalyzeR:::.validate_diff_table(
      missing_direction, fixture$obj@network
    ),
    "<missing>"
  )

  direction_conflict <- rbind(
    fixture$input[1, ],
    transform(fixture$input[1, ], direction = "down")
  )
  expect_error(
    IsoNetAnalyzeR:::.validate_diff_table(
      direction_conflict, fixture$obj@network
    ),
    "one direction"
  )

  gene_conflict <- rbind(
    fixture$input[1, ],
    transform(fixture$input[1, ], gene_id = "G_other")
  )
  expect_error(
    IsoNetAnalyzeR:::.validate_diff_table(gene_conflict, fixture$obj@network),
    "one input gene"
  )

  mapping_conflict <- transform(fixture$input[1, ], gene_id = "wrong_gene")
  expect_error(
    IsoNetAnalyzeR:::.validate_diff_table(
      mapping_conflict, fixture$obj@network
    ),
    "disagrees with the network"
  )
})

test_that("integer and cutoff arguments reject fractional and non-finite values", {
  fixture <- make_isolated_seed_fixture()
  for (bad_order in list(1.5, NA_real_, NaN, Inf)) {
    expect_error(
      virtual_pulldown(fixture$obj, "tx_up", bad_order, cutoff = 0),
      "finite positive integer"
    )
  }
  for (bad_cutoff in list(NA_real_, NaN, Inf, -Inf, 1.1)) {
    expect_error(
      virtual_pulldown(fixture$obj, "tx_up", 1, cutoff = bad_cutoff),
      "finite numeric value"
    )
  }

  expect_error(
    IsoNetAnalyzeR(
      "bad-cutoff", network = fixture$obj@network,
      metadata = list(network_cutoff = Inf)
    ),
    "metadata\\$network_cutoff"
  )
})

test_that("lookup reports partially missing isoform seeds", {
  fixture <- make_isolated_seed_fixture()
  expect_error(
    lookup_virtual_pulldown(
      fixture$obj, "partial", c("tx_up", "not_present"),
      order = 1, cutoff = 0
    ),
    "not_present"
  )
})

test_that("custom network contract is enforced", {
  fixture <- make_isolated_seed_fixture()
  directed <- igraph::as_directed(fixture$obj@network)
  expect_error(IsoNetAnalyzeR("directed", directed), "undirected")

  missing_gene <- igraph::delete_vertex_attr(fixture$obj@network, "gene")
  expect_error(IsoNetAnalyzeR("missing-gene", missing_gene), "gene")

  missing_interact <- igraph::delete_edge_attr(fixture$obj@network, "interact")
  expect_error(IsoNetAnalyzeR("missing-interact", missing_interact), "interact")

  nonbinary <- igraph::set_edge_attr(
    fixture$obj@network, "interact", value = c(1, 0.5)
  )
  expect_error(IsoNetAnalyzeR("nonbinary", nonbinary), "binary")

  blank_name <- igraph::set_vertex_attr(
    fixture$obj@network, "name", value = c("", igraph::V(fixture$obj@network)$name[-1])
  )
  expect_error(IsoNetAnalyzeR("blank-name", blank_name), "vertex names")

  duplicate_name <- igraph::set_vertex_attr(
    fixture$obj@network, "name",
    value = c("duplicate", "duplicate", igraph::V(fixture$obj@network)$name[-c(1, 2)])
  )
  expect_error(IsoNetAnalyzeR("duplicate-name", duplicate_name), "unique")
})

test_that("all-isolated custom networks work throughout", {
  g <- igraph::make_empty_graph(n = 2L, directed = FALSE)
  g <- igraph::set_vertex_attr(g, "name", value = c("iso_a", "iso_b"))
  g <- igraph::set_vertex_attr(g, "gene", value = c("GA", "GB"))
  g <- igraph::set_edge_attr(g, "interact", value = numeric(0))

  obj <- IsoNetAnalyzeR(
    "all-isolated", g, metadata = list(network_cutoff = 0)
  )
  pulled <- virtual_pulldown(obj, "iso_a", order = 1)
  expect_equal(igraph::vcount(pulled), 1L)
  expect_equal(igraph::ecount(pulled), 0L)

  obj <- run_network_diagnostics(obj)
  diagnostics <- get_network_diagnostics(obj)
  expect_equal(diagnostics@network_stats$full_predicted_edges, 0L)
  expect_equal(diagnostics@network_stats$isolated_isoforms, 2L)
  expect_equal(nrow(diagnostics@component_table), 0L)
  expect_equal(diagnostics@isoform_degree_table$degree, c(0L, 0L))
  expect_true(is.na(diagnostics@sticky_isoform_threshold))
})

test_that("an analysis stops clearly when every gene genuinely fails", {
  fixture <- make_isolated_seed_fixture()
  testthat::local_mocked_bindings(
    .process_one_gene = function(gid, ...) list(error = paste("forced", gid)),
    .package = "IsoNetAnalyzeR"
  )
  expect_error(
    suppressWarnings(run_differential_pulldown(
      fixture$obj, "all-fail", fixture$input[1:2, ], order = 1, cutoff = 0,
      confirm = FALSE, verbose = FALSE
    )),
    "Every gene failed.*G1",
    fixed = FALSE
  )
})

test_that("plots can be assigned and saved without a default device file", {
  skip_if_not_installed("ggraph")
  fixture <- make_isolated_seed_fixture()
  obj <- suppressWarnings(run_differential_pulldown(
    fixture$obj, "plot", fixture$input[3, ], order = 1, cutoff = 0,
    confirm = FALSE, verbose = FALSE
  ))
  result <- get_gene_result(obj, "plot", "G3")
  default_plot <- file.path(getwd(), "Rplots.pdf")
  unlink(default_plot)

  visible_plot <- withVisible(suppressMessages(
    plot_network(result, direction = "up")
  ))
  expect_true(visible_plot$visible)
  p <- visible_plot$value
  expect_s3_class(p, "ggplot")
  output <- tempfile(fileext = ".pdf")
  expect_no_error(ggplot2::ggsave(output, plot = p, width = 4, height = 3))
  expect_true(file.exists(output))
  expect_false(file.exists(default_plot))
})
