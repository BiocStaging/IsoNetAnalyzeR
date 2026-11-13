# Expected values follow from the design of the hand-built fixtures in
# helper-fixtures.R.
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
  expect_equal(igraph::V(isolated@up_graph)$name, "tx_isolated_up")
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
  # Isoform level: G1's two seeds and their shared neighbour, and one
  # component for each isolated G3 seed. Gene level: G1-G2, and G3 alone
  # with both seeds.
  iso_overview <- get_component_overview(obj, "test", "iso_components")
  expect_equal(nrow(iso_overview), 3L)
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
    expect_true(ggplot2::is_ggplot(p_isolated_up))
    p_isolated_overlay <- expect_no_warning(
      suppressMessages(plot_network(isolated, direction = "overlay"))
    )
    expect_true(ggplot2::is_ggplot(p_isolated_overlay))

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
    expect_true(ggplot2::is_ggplot(p_component))

    obj <- run_component_analysis(
      obj, "test", "gene_edges_only",
      direction = "both", level = "gene", min_edges = 1
    )
    p_filtered_component <- suppressWarnings(
      suppressMessages(plot_component(
        obj, "test", "gene_edges_only", component_id = 1,
        colour_by = "direction"
      ))
    )
    expect_true(ggplot2::is_ggplot(p_filtered_component))
    expect_setequal(as.character(p_filtered_component$data$name), c("G1", "G2"))

    obj <- lookup_virtual_pulldown(
      obj, "isolated_lookup", seed = "tx_isolated_up", order = 1, cutoff = 0
    )
    lookup_stats <- get_lookup_stats(obj, "isolated_lookup")
    expect_equal(nrow(lookup_stats), 1L)
    expect_true(lookup_stats$is_seed)
    p_lookup <- expect_no_warning(
      suppressMessages(plot_network(get_lookup(obj, "isolated_lookup")))
    )
    expect_true(ggplot2::is_ggplot(p_lookup))

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
    expect_true(ggplot2::is_ggplot(p_restored_lookup))
  })

  expect_false(any(file.exists(default_plot_paths)))
})

test_that("overlay plots colour focal seeds independently from neighbourhood status", {
  skip_if_not_installed("ggraph")

  with_temp_plot_device({
    fixture <- make_isolated_seed_fixture()
    obj <- suppressWarnings(run_differential_pulldown(
      fixture$obj, "seed-colour", fixture$input[1:2, ],
      order = 1, cutoff = 0, confirm = FALSE, verbose = FALSE
    ))
    result <- get_gene_result(obj, "seed-colour", "G1")
    p <- expect_no_warning(
      suppressMessages(plot_network(result, direction = "overlay"))
    )

    seed_rows <- as.character(p$data$seed_direction) != "non_seed"
    expect_true(any(seed_rows))
    expect_true(all(as.character(p$data$node_fill[seed_rows]) == "focal_seed"))
    expect_true(all(as.character(p$data$node_fill[!seed_rows]) == "shared"))
    # Okabe-Ito reddish purple for focal seeds, blue for shared neighbours
    built_nodes <- ggplot2::ggplot_build(p)$data[[2L]]
    expect_equal(built_nodes$fill[seed_rows], rep("#CC79A7", sum(seed_rows)))
    expect_equal(built_nodes$fill[!seed_rows], "#0072B2")
  })
})

test_that("plot graphs keep every edge and isolated seeds without warnings", {
  fixture <- make_isolated_seed_fixture()
  prepared <- expect_no_warning(
    IsoNetAnalyzeR:::.build_plot_graph(
      fixture$obj@network,
      level = "isoform",
      seed_vertices = "tx_isolated_up"
    )
  )

  # Both fixture edges are drawn; the isolated seed stays and the isolated
  # non-seeds tx_isolated_down and tx_other are dropped
  expect_equal(igraph::ecount(prepared), 2L)
  expect_setequal(
    igraph::V(prepared)$name,
    c("tx_up", "tx_down", "tx_neighbour", "tx_isolated_up")
  )
})

test_that("all_edges is deprecated and has no effect", {
  fixture <- make_isolated_seed_fixture()
  obj <- suppressMessages(lookup_virtual_pulldown(
    fixture$obj, "lk", "tx_up", order = 1, cutoff = 0))
  lookup <- get_lookup(obj, "lk")
  default <- expect_no_warning(suppressMessages(plot_network(lookup)))
  expect_warning(named <- suppressMessages(plot_network(lookup,
                                                        all_edges = TRUE)),
                 "deprecated")
  # Third position is still all_edges, not group_by_gene
  expect_warning(positional <- suppressMessages(plot_network(lookup, "isoform",
                                                             TRUE)),
                 "deprecated")
  expect_identical(named$data$name, default$data$name)
  expect_identical(positional$data$name, default$data$name)
})

test_that("patchwork loads with the package, for combined plots", {
  # Its `+`, `|` and `/` methods for ggplots exist only once it is loaded
  expect_true("patchwork" %in% names(getNamespaceImports("IsoNetAnalyzeR")))
})

test_that("distance_from_seed is the distance to the nearest seed", {
  # seed_a has one neighbour, n1; seed_b is isolated
  g <- igraph::make_graph(~ seed_a - n1, seed_b)
  igraph::V(g)$gene <- c("GA", "GN", "GB")
  igraph::V(g)$order <- c(0, 1, 0)

  stats <- analyze_subgraph(g)
  expect_identical(stats$is_seed, c(TRUE, FALSE, TRUE))
  # Each seed is 0 from itself; n1 is 1 hop from seed_a, although seed_b
  # cannot reach it
  expect_equal(stats$distance_from_seed, c(0, 1, 0))
})

test_that("summarize_subgraphs gives one row per graph, named or not", {
  fixture <- make_isolated_seed_fixture()
  up <- virtual_pulldown(fixture$obj, "tx_up", order = 1, cutoff = 0)
  down <- virtual_pulldown(fixture$obj, "tx_down", order = 1, cutoff = 0)
  # Rows follow the list, tx_up's graph then tx_down's, whatever the names
  for (graphs in list(list(up, down), list(a = up, b = down),
                      list(a = up, a = down), list(a = up, down)))
    expect_identical(summarize_subgraphs(graphs)$seed, c("tx_up", "tx_down"))
})

test_that("results saved without is_seed take seeds from the seed lists", {
  fixture <- make_isolated_seed_fixture()
  obj <- suppressWarnings(run_differential_pulldown(
    fixture$obj, "legacy", fixture$input[1:2, ], order = 1, cutoff = 0,
    confirm = FALSE, verbose = FALSE
  ))
  gr <- obj@analyses$legacy@gene_results$G1
  gr@up_stats <- gr@up_stats[, names(gr@up_stats) != "is_seed"]
  gr@down_stats <- gr@down_stats[, names(gr@down_stats) != "is_seed"]
  obj@analyses$legacy@gene_results$G1 <- gr
  stats <- get_isoform_stats(obj, "legacy", "G1")
  # G1's seeds are tx_up (up) and tx_down (down); tx_neighbour is in both
  # neighbourhoods and is no seed
  expect_identical(nrow(stats), 4L)
  expect_setequal(paste(stats$direction, stats$isoform)[stats$is_seed],
                  c("up tx_up", "down tx_down"))
})
