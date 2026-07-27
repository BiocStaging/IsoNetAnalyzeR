make_isolated_seed_fixture <- function() {
  vertices <- data.frame(
    name = c("tx_up", "tx_down", "tx_neighbour", "tx_isolated", "tx_other"),
    gene = c("G1", "G1", "G2", "G3", "G4"),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    from = c("tx_up", "tx_down"),
    to = c("tx_neighbour", "tx_neighbour"),
    interact = c(1L, 1L),
    pi = c(0.9, 0.8),
    pi_score = c(0.9, 0.8),
    stringsAsFactors = FALSE
  )
  network <- igraph::graph_from_data_frame(
    edges, directed = FALSE, vertices = vertices
  )
  obj <- IsoNetAnalyzeR(
    "isolated-seed-test",
    network = network,
    metadata = list(network_cutoff = 0, pi_threshold = 0.5)
  )
  input <- data.frame(
    gene_id = c("G1", "G1", "G3", "G3", "G_missing", "G_missing"),
    isoform_id = c(
      "tx_up", "tx_down", "tx_isolated", "tx_isolated",
      "tx_absent_up", "tx_absent_down"
    ),
    direction = c("up", "down", "up", "down", "up", "down"),
    stringsAsFactors = FALSE
  )
  list(obj = obj, input = input)
}

test_that("isolated and unmapped seeds are distinguished end to end", {
  fixture <- make_isolated_seed_fixture()
  suppressWarnings(
    obj <- run_differential_pulldown(
      fixture$obj, "test", fixture$input, order = 1, cutoff = 0,
      confirm = FALSE, verbose = FALSE
    )
  )

  isolated <- get_gene_result(obj, "test", "G3")
  expect_equal(igraph::vcount(isolated@up_graph), 1L)
  expect_equal(igraph::ecount(isolated@up_graph), 0L)
  expect_equal(igraph::V(isolated@up_graph)$name, "tx_isolated")
  expect_true(isolated@up_stats$is_seed)
  expect_equal(isolated@graph_summary$density, c(0, 0))
  expect_equal(isolated@graph_summary$clustering, c(0, 0))
  expect_true(all(is.na(isolated@graph_summary$starness)))

  missing <- get_gene_result(obj, "test", "G_missing")
  expect_null(missing@up_graph)
  expect_null(missing@down_graph)
  expect_equal(missing@gene_summary$n_up_found, 0L)
  expect_equal(missing@gene_summary$n_down_found, 0L)

  stats <- get_isoform_stats(obj, "test", "G3")
  expect_true(all(stats$is_seed))

  obj <- run_component_analysis(
    obj, "test", "iso_components", direction = "both", level = "isoform"
  )
  iso_overview <- get_component_overview(obj, "test", "iso_components")
  expect_equal(nrow(iso_overview), 2L)
  expect_true(any(iso_overview$n_nodes == 1L & iso_overview$n_edges == 0L))

  obj <- run_component_analysis(
    obj, "test", "gene_components", direction = "both", level = "gene"
  )
  gene_overview <- get_component_overview(obj, "test", "gene_components")
  expect_equal(nrow(gene_overview), 2L)
  singleton <- gene_overview[gene_overview$n_nodes == 1L, , drop = FALSE]
  expect_equal(singleton$n_seeds_up, 1L)
  expect_equal(singleton$n_seeds_down, 1L)
  expect_equal(singleton$n_genes_up, 1L)
  expect_equal(singleton$n_genes_down, 1L)
  expect_true(singleton$mixed)

  obj <- run_component_analysis(
    obj, "test", "gene_edges_only", direction = "both", level = "gene",
    min_edges = 1
  )
  filtered <- get_component_overview(obj, "test", "gene_edges_only")
  expect_equal(nrow(filtered), 1L)
  expect_true(all(filtered$n_edges >= 1L))
})

test_that("isolated seed plots and serialization remain valid", {
  skip_if_not_installed("ggraph")
  default_plot_paths <- c(
    file.path(getwd(), "Rplots.pdf"),
    file.path(getwd(), "tests", "testthat", "Rplots.pdf")
  )

  with_temp_plot_device({
    fixture <- make_isolated_seed_fixture()
    suppressWarnings(
      obj <- run_differential_pulldown(
        fixture$obj, "test", fixture$input, order = 1, cutoff = 0,
        confirm = FALSE, verbose = FALSE
      )
    )
    obj <- run_component_analysis(
      obj, "test", "gene_components", direction = "both", level = "gene"
    )

    isolated <- get_gene_result(obj, "test", "G3")
    p_isolated_up <- expect_no_warning(
      suppressMessages(plot_network(isolated, direction = "up"))
    )
    expect_s3_class(p_isolated_up, "ggplot")
    p_isolated_overlay <- expect_no_warning(
      suppressMessages(plot_network(isolated, direction = "overlay"))
    )
    expect_s3_class(p_isolated_overlay, "ggplot")

    singleton_id <- get_component_overview(
      obj, "test", "gene_components"
    )$component_id[
      get_component_overview(obj, "test", "gene_components")$n_nodes == 1L
    ]
    p_component <- expect_no_warning(
      suppressMessages(plot_component(
        obj, "test", "gene_components", singleton_id,
        colour_by = "direction"
      ))
    )
    expect_s3_class(p_component, "ggplot")

    obj <- lookup_virtual_pulldown(
      obj, "isolated_lookup", seed = "tx_isolated", order = 1, cutoff = 0
    )
    lookup_stats <- get_lookup_stats(obj, "isolated_lookup")
    expect_equal(nrow(lookup_stats), 1L)
    expect_true(lookup_stats$is_seed)
    p_lookup <- expect_no_warning(
      suppressMessages(plot_network(get_lookup(obj, "isolated_lookup")))
    )
    expect_s3_class(p_lookup, "ggplot")

    path <- tempfile(fileext = ".rds")
    saveRDS(obj, path)
    restored <- readRDS(path)
    expect_equal(
      get_component_overview(restored, "test", "gene_components"),
      get_component_overview(obj, "test", "gene_components")
    )
    p_restored_lookup <- expect_no_warning(
      suppressMessages(plot_network(get_lookup(restored, "isolated_lookup")))
    )
    expect_s3_class(p_restored_lookup, "ggplot")
  })

  expect_false(any(file.exists(default_plot_paths)))
})

test_that("plot graph filtering preserves isolated seeds without warnings", {
  fixture <- make_isolated_seed_fixture()
  filtered <- expect_no_warning(
    IsoNetAnalyzeR:::.build_plot_graph(
      fixture$obj@network,
      level = "isoform",
      all_edges = FALSE,
      threshold = 0.85,
      seed_vertices = "tx_isolated"
    )
  )

  expect_equal(igraph::ecount(filtered), 1L)
  expect_setequal(
    igraph::V(filtered)$name,
    c("tx_up", "tx_neighbour", "tx_isolated")
  )
  expect_equal(igraph::E(filtered)$pi_score, 0.9)
})

test_that("explicit seed flags survive disconnected multi-seed distances", {
  g <- igraph::make_empty_graph(n = 2, directed = FALSE)
  igraph::V(g)$name <- c("seed_a", "seed_b")
  igraph::V(g)$gene <- c("GA", "GB")
  igraph::V(g)$order <- c(0, 0)
  igraph::V(g)$frac_internal <- c(NA_real_, NA_real_)

  stats <- analyze_subgraph(g)
  expect_true(all(stats$is_seed))
  expect_true(all(is.infinite(stats$distance_from_seed)))

  legacy <- stats[, setdiff(names(stats), "is_seed"), drop = FALSE]
  expect_false("is_seed" %in% names(legacy))

  fixture <- make_isolated_seed_fixture()
  suppressWarnings(
    obj <- run_differential_pulldown(
      fixture$obj, "legacy", fixture$input, order = 1, cutoff = 0,
      confirm = FALSE, verbose = FALSE
    )
  )
  legacy_result <- obj@analyses[["legacy"]]@gene_results[["G3"]]
  legacy_result@up_stats <- legacy_result@up_stats[
    , setdiff(names(legacy_result@up_stats), "is_seed"), drop = FALSE
  ]
  obj@analyses[["legacy"]]@gene_results[["G3"]] <- legacy_result
  legacy_stats <- get_isoform_stats(obj, "legacy", "G3")
  expect_true(all(legacy_stats$is_seed))
})
