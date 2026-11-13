# Gene results are built by hand so that the expected component counts can
# be read off the graphs below.
make_gene_result <- function(gene_id, up_isoforms, down_isoforms) {
  empty_graph <- igraph::make_empty_graph(directed = FALSE)
  methods::new(
    "GeneDiffPulldown",
    gene_id = gene_id,
    up_isoforms = up_isoforms,
    down_isoforms = down_isoforms,
    up_graph = empty_graph,
    down_graph = empty_graph,
    up_stats = S4Vectors::DataFrame(),
    down_stats = S4Vectors::DataFrame(),
    graph_summary = S4Vectors::DataFrame(),
    gene_summary = list()
  )
}

test_that("gene-level component summaries count ENSG seed genes", {
  graph <- igraph::make_graph(c("ENSG000001", "ENSG000002"), directed = FALSE)
  igraph::V(graph)$source_gene <- c("ENSG000001", "multi-gene")
  igraph::V(graph)$seed_direction <- c("both", NA_character_)

  gene_results <- list(
    ENSG000001 = make_gene_result(
      "ENSG000001",
      up_isoforms = "ENST_UP_1",
      down_isoforms = "ENST_DOWN_1"
    )
  )

  component_table <- IsoNetAnalyzeR:::.build_component_table(
    list(graph),
    gene_results,
    direction = "both",
    level = "gene"
  )

  expect_equal(component_table$n_seeds_up, 1L)
  expect_equal(component_table$n_seeds_down, 1L)
  expect_equal(component_table$n_genes_up, 1L)
  expect_equal(component_table$n_genes_down, 1L)
  expect_true(component_table$mixed)
})

test_that("isoform-level component summaries keep ENST matching", {
  graph <- igraph::make_graph(c("ENST_UP_1", "ENST_DOWN_1"), directed = FALSE)
  igraph::V(graph)$source_gene <- c("ENSG000001", "ENSG000001")
  igraph::V(graph)$seed_direction <- c("up", "down")

  gene_results <- list(
    ENSG000001 = make_gene_result(
      "ENSG000001",
      up_isoforms = "ENST_UP_1",
      down_isoforms = "ENST_DOWN_1"
    )
  )

  component_table <- IsoNetAnalyzeR:::.build_component_table(
    list(graph),
    gene_results,
    direction = "both",
    level = "isoform"
  )

  expect_equal(component_table$n_genes_up, 1L)
  expect_equal(component_table$n_genes_down, 1L)
  expect_true(component_table$mixed)
})

test_that("merged graphs record which seed genes reach each node", {
  fx <- make_multi_gene_fixture()
  results <- get_analysis(fx$obj, "a")@gene_results
  source_of <- function(level, nodes) {
    g <- .build_cross_gene_graph(results, "both", level)$graph
    igraph::V(g)$source_gene[match(nodes, igraph::V(g)$name)]
  }
  # N4 is reached only from GB's seeds and N7 only from GC's; N1 from GA's
  # seeds and GF's down seed; GA's seeds only from GA
  expect_identical(source_of("isoform", c("tx_N4", "tx_N7", "tx_N1", "GA_up")),
                   c("GB", "GC", "multi-gene", "GA"))
  expect_identical(source_of("gene", c("N4", "N7", "N1", "GA")),
                   c("GB", "GC", "multi-gene", "GA"))
})
