test_that("bundled network keeps every vertex and only predicted interactions", {
  ns <- asNamespace("IsoNetAnalyzeR")
  g <- get(".network_data", envir = ns)
  expect_equal(igraph::vcount(g), 58164L)
  expect_equal(igraph::ecount(g), 158458L)
  expect_true(all(igraph::E(g)$interact == 1L))
  expect_setequal(igraph::edge_attr_names(g), c("pi", "interact"))
  expect_setequal(igraph::vertex_attr_names(g), c("name", "gene"))
  expect_equal(unname(get(".network_cutoff", envir = ns)), 0.2608696,
               tolerance = 1e-6)
  expect_null(get(".pi_threshold", envir = ns))
})

test_that("default session applies eligibility to the reduced network", {
  obj <- suppressMessages(IsoNetAnalyzeR("bundled"))
  audit <- get_metadata(obj, "network_biotype_filter")
  expect_equal(audit$vertices_before, 58164L)
  expect_equal(audit$vertices_after, 46433L)
  expect_equal(audit$excluded_nmd, 11662L)
  expect_equal(audit$excluded_nsd, 69L)
  expect_equal(igraph::vcount(get_network(obj)), 46433L)
})

test_that("CLTCL1 example reproduces the published neighbourhood sizes", {
  seeds <- utils::read.csv(system.file(
    "extdata", "gse141432_t2d_vs_lean_seeds.csv", package = "IsoNetAnalyzeR"
  ))
  expect_equal(nrow(seeds), 137L)
  expect_equal(length(unique(seeds$gene_id)), 65L)
  obj <- suppressMessages(IsoNetAnalyzeR("cltcl1"))
  obj <- suppressMessages(run_differential_pulldown(
    obj, "t2d", seeds[seeds$gene_name == "CLTCL1", ], order = 1,
    cutoff = 0.2609, confirm = FALSE, verbose = FALSE,
    seed_padj_col = "within_gene_padj", seed_fdr = 0.05
  ))
  gr <- get_gene_result(obj, "t2d", "ENSG00000070371")
  expect_equal(igraph::ecount(get_pulldown_graph(gr, "up")), 8L)
  expect_equal(igraph::ecount(get_pulldown_graph(gr, "down")), 33L)
})

test_that("metadata, network and pulldown-graph accessors", {
  fixture <- make_isolated_seed_fixture()
  obj <- set_metadata(fixture$obj, "note", "pilot")
  expect_identical(get_metadata(obj, "note"), "pilot")
  expect_type(get_metadata(obj), "list")
  expect_true("network_biotype_filter" %in% names(get_metadata(obj)))
  expect_error(get_metadata(obj, "absent"), "No metadata entry 'absent'")
  expect_error(get_metadata(obj, c("a", "b")), "single non-empty")
  expect_error(get_metadata(list(), "note"), "IsoNetAnalyzeR object")

  expect_s3_class(get_network(obj), "igraph")
  expect_equal(igraph::vcount(get_network(obj)), 6L)
  expect_error(get_network(list()), "IsoNetAnalyzeR object")

  suppressWarnings(
    obj <- run_differential_pulldown(
      obj, "test", fixture$input, order = 1, cutoff = 0,
      confirm = FALSE, verbose = FALSE
    )
  )
  expect_true(
    "test" %in% names(get_metadata(obj, "analysis_network_coverage"))
  )
  g1 <- get_gene_result(obj, "test", "G1")
  expect_identical(get_pulldown_graph(g1, "up"), g1@up_graph)
  expect_identical(get_pulldown_graph(g1, "down"), g1@down_graph)
  expect_identical(get_pulldown_graph(g1), g1@up_graph)
  expect_error(get_pulldown_graph(g1, "sideways"))
  missing <- get_gene_result(obj, "test", "G_missing")
  expect_null(get_pulldown_graph(missing, "down"))
})

test_that("biotype table is cached without global assignment", {
  first <- isonet_transcript_biotypes()
  expect_identical(isonet_transcript_biotypes(), first)
  expect_equal(nrow(first), 250161L)
})

test_that("gene symbols cover the bundled network", {
  gene_names <- isonet_gene_names()
  expect_named(gene_names, c("gene_id", "gene_name"))
  expect_false(anyDuplicated(gene_names$gene_id) > 0L)
  g <- get(".network_data", envir = asNamespace("IsoNetAnalyzeR"))
  expect_true(all(unique(igraph::V(g)$gene) %in% gene_names$gene_id))
  expect_identical(
    gene_names$gene_name[gene_names$gene_id == "ENSG00000070371"], "CLTCL1"
  )
  expect_identical(.gene_label(c("ENSG00000070371.9", "gene_X")),
                   c("CLTCL1", "gene_X"))
  expect_identical(.gene_title("ENSG00000070371"),
                   "CLTCL1 (ENSG00000070371)")
})

test_that("node labels use symbols, identifiers or nothing", {
  g <- igraph::make_graph(~ tx1 - tx2 - tx3)
  igraph::V(g)$gene <- c("ENSG00000070371", "ENSG00000070371", "gene_X")
  expect_identical(igraph::V(.add_node_labels(g, "symbol"))$node_label,
                   c("CLTCL1", "CLTCL1", "tx3"))
  expect_identical(igraph::V(.add_node_labels(g, "id"))$node_label,
                   c("tx1", "tx2", "tx3"))
  igraph::V(g)$gene_name <- c("A", "B", "")
  expect_identical(igraph::V(.add_node_labels(g, "symbol"))$node_label,
                   c("A", "B", "tx3"))
  expect_null(.label_geom(2.5, "none"))

  seeds <- utils::read.csv(system.file(
    "extdata", "gse141432_t2d_vs_lean_seeds.csv", package = "IsoNetAnalyzeR"
  ))
  obj <- suppressMessages(IsoNetAnalyzeR("labels"))
  obj <- suppressMessages(run_differential_pulldown(
    obj, "t2d", seeds[seeds$gene_name == "CLTCL1", ], order = 1,
    confirm = FALSE, verbose = FALSE
  ))
  gr <- get_gene_result(obj, "t2d", "ENSG00000070371")
  with_labels <- plot_network(gr, direction = "overlay")
  without <- plot_network(gr, direction = "overlay", node_labels = "none")
  expect_equal(length(with_labels$layers), length(without$layers) + 1L)
  expect_true("CLTCL1" %in% with_labels$data$node_label)
  expect_match(with_labels$labels$title, "CLTCL1 (ENSG00000070371)",
               fixed = TRUE)
  expect_error(plot_network(gr, node_labels = "gene"))
})

test_that("printed summaries keep their column alignment", {
  x <- 5
  expect_identical(cli::cli_fmt(.cli_aligned("Name       : {x}")),
                   "Name       : 5")
  expect_identical(cli::cli_fmt(.cli_aligned("Up : {x}
                                             and more")),
                   "Up : 5 and more")
  obj <- suppressMessages(IsoNetAnalyzeR("printing"))
  out <- cli::cli_fmt(show(obj))
  expect_true(any(grepl("Network     : 46,433 transcripts, 146,267 edges",
                        out, fixed = TRUE)))
})

test_that("result hints name the caller's session variable", {
  expect_identical(.obj_name(quote(my_session)), "my_session")
  expect_identical(.obj_name(quote(list()[[1]])), "obj")
  fixture <- make_isolated_seed_fixture()
  my_session <- fixture$obj
  expect_message(
    suppressWarnings(run_differential_pulldown(
      my_session, "test", fixture$input, order = 1, cutoff = 0,
      confirm = FALSE, verbose = FALSE
    )),
    'get_analysis(my_session, "test")', fixed = TRUE
  )
})
