# Expected values were fixed outside the package: network counts from the
# 0.6.5 network file, biotype counts from the Ensembl 107 annotation, the
# example seed table, and the CLTCL1 neighbourhood sizes of the study.
# inst/validation/test_ground_truth.R recomputes all of them without package
# code, including the CLTCL1 and APAF1 neighbourhoods with a plain igraph
# implementation of the documented pulldown.
test_that("bundled network keeps every vertex and only predicted interactions", {
  # 58,164 transcripts; 158,458 interact == 1 rows of the 0.6.5 network, less
  # 6 repeated pairs and 5 self-pairs (inst/scripts/reduce_network.R)
  ns <- asNamespace("IsoNetAnalyzeR")
  g <- get(".network_data", envir = ns)
  expect_equal(igraph::vcount(g), 58164L)
  expect_equal(igraph::ecount(g), 158447L)
  expect_true(igraph::is_simple(g))
  expect_true(all(igraph::E(g)$interact == 1L))
  expect_setequal(igraph::edge_attr_names(g), c("pi", "interact"))
  expect_setequal(igraph::vertex_attr_names(g), c("name", "gene"))
  expect_equal(unname(get(".network_cutoff", envir = ns)), 0.2608696,
               tolerance = 1e-6)
  expect_null(get(".pi_threshold", envir = ns))
})

test_that("default session applies eligibility to the reduced network", {
  skip_if_no_ensembl()
  obj <- suppressMessages(IsoNetAnalyzeR("bundled"))
  audit <- get_metadata(obj, "network_biotype_filter")
  expect_equal(audit$vertices_before, 58164L)
  expect_equal(audit$vertices_after, 46433L)
  expect_equal(audit$excluded_nmd, 11662L)
  expect_equal(audit$excluded_nsd, 69L)
  expect_equal(igraph::vcount(get_network(obj)), 46433L)
})

test_that("CLTCL1 example reproduces the published neighbourhood sizes", {
  skip_if_no_ensembl()
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

test_that("repeated source pairs are merged and self-pairs removed", {
  g <- get(".network_data", envir = asNamespace("IsoNetAnalyzeR"))
  # PRPF18-PRPF8 was listed twice in the source table, with pi 0.5911 and
  # 0.5594
  e <- igraph::get_edge_ids(g, c("ENST00000378572", "ENST00000304992"))
  expect_equal(igraph::E(g)$pi[e], mean(c(0.5911, 0.5594)))
  expect_false(any(igraph::which_loop(g)))
})

test_that("APAF1 neighbourhoods match the simplified network", {
  skip_if_no_ensembl()
  # In the 0.99.0 network the DIABLO transcript ENST00000464942 entered the
  # up-labelled APAF1 neighbourhood only because its repeated edges and its
  # self-pair inflated the degrees behind frac_internal. Expected sizes are
  # recomputed without package code in inst/validation/test_ground_truth.R.
  seeds <- utils::read.csv(system.file(
    "extdata", "gse141432_t2d_vs_lean_seeds.csv", package = "IsoNetAnalyzeR"
  ))
  obj <- suppressMessages(IsoNetAnalyzeR("apaf1"))
  obj <- suppressMessages(run_differential_pulldown(
    obj, "t2d", seeds[seeds$gene_name == "APAF1", ], order = 1,
    confirm = FALSE, verbose = FALSE,
    seed_padj_col = "within_gene_padj", seed_fdr = 0.05
  ))
  gr <- get_gene_result(obj, "t2d", "ENSG00000120868")
  up <- get_pulldown_graph(gr, "up")
  down <- get_pulldown_graph(gr, "down")
  expect_false("ENST00000464942" %in% igraph::V(up)$name)
  expect_equal(c(igraph::vcount(up), igraph::ecount(up)), c(3, 2))
  expect_equal(c(igraph::vcount(down), igraph::ecount(down)), c(9, 11))
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
  # The isolated-seed fixture has six transcripts
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

test_that("transcript biotypes come from Ensembl 107 and are cached", {
  skip_if_no_ensembl()
  first <- isonet_transcript_biotypes()
  expect_identical(isonet_transcript_biotypes(), first)
  expect_s4_class(first, "DataFrame")
  expect_named(first, c("isoform_id", "transcript_biotype"))
  # Every network transcript is annotated; the counts are those of the
  # GENCODE v41 table the network was built with (inst/PROVENANCE.md)
  g <- get(".network_data", envir = asNamespace("IsoNetAnalyzeR"))
  biotype <- first$transcript_biotype[match(igraph::V(g)$name,
                                            first$isoform_id)]
  expect_false(anyNA(biotype))
  expect_equal(as.vector(table(biotype)[c("protein_coding",
                                          "nonsense_mediated_decay",
                                          "non_stop_decay")]),
               c(46433L, 11662L, 69L))
})

test_that("gene symbols cover the bundled network", {
  skip_if_no_ensembl()
  gene_names <- isonet_gene_names()
  expect_s4_class(gene_names, "DataFrame")
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

test_that("networks with Ensembl transcript IDs always use Ensembl 107", {
  skip_if_no_ensembl()
  ann <- isonet_transcript_biotypes()
  decay <- ann$isoform_id[ann$transcript_biotype %in%
                            c("nonsense_mediated_decay", "non_stop_decay")]
  # ENST00000477548 is an NMD transcript outside the network whose
  # protein-identical alias would map to a coding representative
  g <- get(".network_data", envir = asNamespace("IsoNetAnalyzeR"))
  expect_true("ENST00000477548" %in% decay)
  expect_false("ENST00000477548" %in% igraph::V(g)$name)

  # A table annotating every vertex still excludes every annotated decay
  # transcript of Ensembl 107, as the default session does
  complete <- ann[ann$isoform_id %in% igraph::V(g)$name, ]
  obj <- suppressMessages(IsoNetAnalyzeR("table", g,
                                         transcript_biotypes = complete))
  expect_setequal(
    igraph::graph_attr(get_network(obj), "excluded_transcript_ids"), decay
  )
  resolved <- resolve_duplicate_isoforms("ENST00000477548", get_network(obj))
  expect_identical(attr(resolved, "dedup_audit")$resolution,
                   "excluded_transcript_biotype")

  # With both exclusions off, every bundled transcript is still annotated
  both <- suppressMessages(IsoNetAnalyzeR("both", include_nmd = TRUE,
                                          include_nsd = TRUE))
  audit <- get_metadata(both, "network_biotype_filter")
  expect_equal(audit$unknown_vertices, 0L)
  expect_equal(audit$vertices_after, 58164L)
  expect_match(audit$annotation, "Ensembl 107")
})

test_that("GENCODE PAR_Y gene IDs get their symbol", {
  skip_if_no_ensembl()
  # ENSG00000002586 is CD99; GENCODE marks its Y copy with _PAR_Y
  expect_identical(.gene_label("ENSG00000002586.20_PAR_Y"), "CD99")
})

test_that("a failed annotation load is retried only once it has loaded", {
  # The symbol table loads only once the EnsDb is in the session cache
  calls <- 0L
  local_mocked_bindings(isonet_gene_names = function(edb = NULL) {
    calls <<- calls + 1L
    if (is.null(.annotation_cache$ensdb)) stop("offline")
    S4Vectors::DataFrame(gene_id = "ENSG00000070371", gene_name = "CLTCL1")
  })
  cached <- mget(c("ensdb", "gene_names", "symbols_failed"),
                 envir = .annotation_cache, ifnotfound = list(NULL))
  withr::defer(list2env(cached, envir = .annotation_cache))
  .annotation_cache$ensdb <- NULL
  .annotation_cache$gene_names <- NULL
  .annotation_cache$symbols_failed <- NULL
  # Two label lookups, one load attempt; the IDs are kept
  expect_identical(.gene_label("ENSG00000070371"), "ENSG00000070371")
  expect_identical(.gene_label("ENSG00000070371"), "ENSG00000070371")
  expect_identical(calls, 1L)
  # Once the EnsDb is loaded by another route (a session), symbols return
  .annotation_cache$ensdb <- "loaded"
  expect_identical(.gene_label("ENSG00000070371"), "CLTCL1")
})

test_that("node labels use symbols, identifiers or nothing", {
  skip_if_no_ensembl()
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
  skip_if_no_ensembl()
  obj <- suppressMessages(IsoNetAnalyzeR("printing"))
  out <- cli::cli_fmt(show(obj))
  expect_true(any(grepl("Network     : 46,433 transcripts, 146,256 edges",
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
