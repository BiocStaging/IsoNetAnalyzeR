# Expected values are known properties of Zachary's karate club graph
# (igraph::make_graph("Zachary")): 34 vertices, 78 edges, diameter 5, and
# two hubs of degree 17 and 16. Pairs of vertices form 17 "genes".
make_karate_session <- function() {
  g <- igraph::make_graph("Zachary")
  n <- igraph::vcount(g)
  igraph::V(g)$name <- paste0("tx", seq_len(n))
  igraph::V(g)$gene <- paste0("gene", ceiling(seq_len(n) / 2))
  igraph::V(g)$transcript_biotype <- "protein_coding"
  igraph::E(g)$interact <- 1L
  suppressMessages(IsoNetAnalyzeR("karate", network = g,
                                  metadata = list(network_cutoff = 0.25)))
}

test_that("network diagnostics summarise a known graph", {
  toy <- suppressMessages(run_network_diagnostics(make_karate_session()))
  d <- get_network_diagnostics(toy)
  expect_s4_class(d, "NetworkDiagnostics")

  s <- d@network_stats
  expect_equal(s$pred_isoforms, 34)
  expect_equal(s$pred_genes, 17)
  expect_equal(s$pred_edges, 78)
  expect_equal(s$pred_density, 78 / choose(34, 2))
  expect_equal(s$pred_diameter, 5)
  expect_equal(s$pred_n_components, 1)
  expect_equal(s$isolated_isoforms, 0)

  # The two hubs of the karate club network have degrees 17 and 16. The
  # default rule keeps degrees above the 0.95 quantile, 12 + 0.35 * (16 - 12)
  # = 13.4, which only they exceed.
  expect_identical(d@sticky_method, "quantile")
  expect_setequal(d@sticky_isoform_table$isoform_id, c("tx1", "tx34"))

  out <- cli::cli_fmt(show(d))
  expect_true(any(grepl("Sticky isoforms *: 2", out)))
  expect_true(any(grepl("All connected *: 17 genes", out)))
})

test_that("the standard-deviation rule replaces earlier diagnostics", {
  toy <- suppressMessages(run_network_diagnostics(make_karate_session()))
  expect_warning(
    toy <- suppressMessages(run_network_diagnostics(
      toy, sticky_method = "sd", sticky_cutoff = 2
    )),
    "Overwriting"
  )
  d <- get_network_diagnostics(toy)
  deg <- igraph::degree(get_network(toy))
  expect_identical(d@sticky_method, "sd")
  expect_equal(d@sticky_isoform_threshold, mean(deg) + 2 * stats::sd(deg))
  expect_true(all(d@sticky_isoform_table$degree > d@sticky_isoform_threshold))
})

test_that("diagnostics count isolated isoforms and partly connected genes", {
  fx <- make_multi_gene_fixture()
  obj <- suppressMessages(run_network_diagnostics(fx$obj))
  d <- get_network_diagnostics(obj)
  # GE's down isoform has no predicted interaction, so GE is the one gene
  # with both connected and isolated isoforms
  expect_equal(d@network_stats$isolated_isoforms, 1)
  out <- cli::cli_fmt(show(d))
  expect_true(any(grepl("Partial *: 1 gene$", out)))
})

test_that("diagnostics plots draw and require a diagnostics run", {
  toy <- make_karate_session()
  expect_error(get_network_diagnostics(toy), "run_network_diagnostics")
  expect_error(plot_network_diagnostics(toy), "run_network_diagnostics")

  toy <- suppressMessages(run_network_diagnostics(toy))
  p <- plot_network_diagnostics(toy, title = "Karate club")
  expect_true(inherits(p, "patchwork"))
  with_temp_plot_device(expect_no_error(print(p)))
})
