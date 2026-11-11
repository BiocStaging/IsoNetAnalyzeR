# Expected values follow from the multi-gene fixture in helper-fixtures.R.
# With cutoff 0 and no edges between neighbours, a seed's neighbourhood is the
# seed and its neighbours: up / down nodes are GA 4 / 2, GB 2 / 4, GC 3 / 2,
# GD 2 / 2, GE 3 / 1 and GF 2 / 3.
test_that("fixture neighbourhood sizes match the design", {
  fx <- make_multi_gene_fixture()
  cgs <- get_cross_gene_summary(fx$obj, "a")
  delta <- setNames(cgs$up_nodes - cgs$down_nodes, cgs$gene_id)
  expect_equal(delta[c("GA", "GB", "GC", "GD", "GE", "GF")],
               c(GA = 2, GB = -2, GC = 1, GD = 0, GE = 2, GF = -1))
})

test_that("run_fgsea ranks genes by a preset and stores the result", {
  fx <- make_multi_gene_fixture()
  withr::local_seed(1)
  obj <- suppressWarnings(suppressMessages(run_fgsea(
    fx$obj, "a", label = "nodes", pathways = fx$pathways, minSize = 2
  )))

  res <- get_gsea_result(obj, "a", "nodes")
  expect_s4_class(res, "FgseaResult")
  expect_identical(res@rank_by, "delta_nodes")
  expect_equal(res@params$n_genes_ranked, 6)

  # delta_nodes, highest first: GA 2 and GE 2 (tied), GC 1, GD 0, GF -1, GB -2
  stats <- get_gsea_stats(obj, "a", "nodes")
  expect_identical(names(stats), c("gene_id", "stat", "rank"))
  expect_equal(stats$stat, c(2, 2, 1, 0, -1, -2))
  expect_setequal(stats$gene_id[1:2], c("GA", "GE"))
  expect_identical(stats$gene_id[3:6], c("GC", "GD", "GF", "GB"))
  expect_equal(stats$rank, 1:6)

  tbl <- get_gsea_table(obj, "a", "nodes")
  # neighbour_set has no ranked genes, so fgsea drops it. up_set holds the
  # three highest-ranked genes and down_set the two lowest, so their
  # enrichment scores are positive and negative.
  expect_setequal(tbl$pathway, c("up_set", "down_set", "mixed_set"))
  expect_gt(tbl$NES[tbl$pathway == "up_set"], 0)
  expect_lt(tbl$NES[tbl$pathway == "down_set"], 0)
  expect_false(is.unsorted(tbl$padj))
  expect_true(all(tbl$padj >= 0 & tbl$padj <= 1))

  expect_identical(nrow(get_gsea_table(obj, "a", "nodes", padj_cutoff = 1)),
                   nrow(tbl))
  # With nperm = 1000 no p-value can fall below 1 / 1001
  expect_message(none <- get_gsea_table(obj, "a", "nodes", padj_cutoff = 1e-6),
                 "No pathways")
  expect_identical(nrow(none), 0L)

  out <- cli::cli_fmt(show(res))
  expect_true(any(grepl("Rank by *: delta_nodes", out)))
})

test_that("run_fgsea accepts other presets and custom statistics", {
  fx <- make_multi_gene_fixture()
  withr::local_seed(1)
  obj <- suppressWarnings(suppressMessages(run_fgsea(
    fx$obj, "a", label = "edges", rank_by = "delta_edges",
    pathways = fx$pathways, minSize = 2
  )))
  expect_identical(get_gsea_result(obj, "a", "edges")@rank_by, "delta_edges")

  # GD's NA is dropped; the others rank GA 3, GE 2, GC 1, GF -1, GB -3
  custom <- c(GA = 3, GB = -3, GC = 1, GD = NA, GE = 2, GF = -1)
  msgs <- capture_messages(obj <- suppressWarnings(run_fgsea(
    obj, "a", label = "custom", stats = custom, pathways = fx$pathways,
    minSize = 2
  )))
  expect_true(any(grepl("Dropping 1 gene", msgs)))
  res <- get_gsea_result(obj, "a", "custom")
  expect_identical(res@rank_by, "custom")
  expect_equal(res@params$n_genes_dropped_na, 1)
  expect_identical(get_gsea_stats(obj, "a", "custom")$gene_id,
                   c("GA", "GE", "GC", "GF", "GB"))
  expect_setequal(names(get_analysis(obj, "a")@gsea_results),
                  c("edges", "custom"))
})

test_that("run_fgsea rejects invalid arguments", {
  fx <- make_multi_gene_fixture()
  obj <- fx$obj
  pw <- fx$pathways
  expect_error(run_fgsea(obj, "a", label = "", pathways = pw), "label")
  expect_error(run_fgsea(obj, "a", label = "x", pathways = pw,
                         minSize = 5, maxSize = 2), "maxSize")
  expect_error(run_fgsea(obj, "a", label = "x", stats = c(1, 2),
                         pathways = pw), "named numeric")
  expect_error(suppressMessages(run_fgsea(
    obj, "a", label = "x", stats = c(GA = NA_real_), pathways = pw
  )), "non-NA")
  expect_error(suppressMessages(run_fgsea(
    obj, "a", label = "x", pathways = list(c("GA", "GB"))
  )), "named list")
  expect_error(run_fgsea(obj, "a", label = "x", rank_by = "nonsense",
                         pathways = pw))
  expect_error(run_fgsea(list(), "a", label = "x", pathways = pw),
               "IsoNetAnalyzeR")
  expect_error(get_gsea_result(obj, "a", "missing"))
})

test_that("component ORA reproduces hypergeometric p-values", {
  fx <- make_multi_gene_fixture()
  obj <- suppressMessages(run_component_analysis(fx$obj, "a", label = "comp"))
  obj <- suppressMessages(run_component_ora(
    obj, "a", "comp", label = "ora", pathways = fx$pathways, minSize = 2
  ))
  ora <- get_component_ora_result(obj, "a", "comp", "ora")
  expect_s4_class(ora, "ComponentOraResult")
  # Universe: the 14 genes of the network, GA-GF and N1-N8. Components: the
  # GA, GB and GD seeds, GE's up seed, GF's down seed and N1-N6 (11 genes);
  # GC's seeds, GF's up seed, N7 and N8 (4 genes); GE's isolated down seed.
  expect_equal(ora@params$universe_size, 14)
  expect_equal(ora@params$n_components_tested, 3)
  rt <- ora@result_table

  # Component 3 is GE's isolated down isoform: one gene drawn from 14,
  # tested against the three up_set genes.
  ge <- rt[rt$component_id == 3 & rt$pathway == "up_set", ]
  # p-values are stored to six decimals
  expect_equal(ge$pval, stats::phyper(0, 3, 11, 1, lower.tail = FALSE),
               tolerance = 1e-5)
  expect_equal(ge$overlap, 1)

  # Component 1 holds 11 genes, including all four mixed_set genes.
  mixed <- rt[rt$component_id == 1 & rt$pathway == "mixed_set", ]
  expect_equal(mixed$pval, stats::phyper(3, 4, 10, 11, lower.tail = FALSE),
               tolerance = 1e-5)
  expect_setequal(mixed$overlap_genes[[1]], c("GA", "GB", "GD", "N4"))
  expect_equal(unique(rt$n_genes_tested[rt$component_id == 1]), 11)

  # One top set per component; only up_set overlaps component 3. All four
  # sets have at least minSize = 2 genes in the universe.
  top <- get_ora_table(obj, "a", "comp", "ora", top_n = 1)
  expect_identical(nrow(top), 3L)
  expect_identical(top$pathway[top$component_id == 3], "up_set")
  first <- get_ora_table(obj, "a", "comp", "ora", component_id = 1)
  expect_true(all(first$component_id == 1))
  expect_identical(nrow(first), 4L)

  out <- cli::cli_fmt(show(ora))
  expect_true(any(grepl("Universe size *: 14 genes", out)))
})

test_that("component ORA supports seed genes and a custom universe", {
  fx <- make_multi_gene_fixture()
  obj <- suppressMessages(run_component_analysis(fx$obj, "a", label = "comp"))
  obj <- suppressMessages(run_component_ora(
    obj, "a", "comp", label = "seeds", pathways = fx$pathways,
    nodes = "seeds", minSize = 2
  ))
  rt <- get_component_ora_result(obj, "a", "comp", "seeds")@result_table
  # Component 1 contains the seeds of GA, GB and GD (both directions), GE's
  # up seed and GF's down seed: five genes
  expect_equal(unique(rt$n_genes_tested[rt$component_id == 1]), 5)

  genes <- c("GA", "GB", "GC", "GD", "GE", "GF")
  obj <- suppressMessages(run_component_ora(
    obj, "a", "comp", label = "custom", pathways = fx$pathways,
    universe = genes, minSize = 2
  ))
  ora <- get_component_ora_result(obj, "a", "comp", "custom")
  expect_equal(ora@params$universe_size, 6)
})

test_that("component ORA rejects invalid arguments", {
  fx <- make_multi_gene_fixture()
  obj <- suppressMessages(run_component_analysis(fx$obj, "a", label = "comp"))
  pw <- fx$pathways
  expect_error(run_component_ora(obj, "a", "comp", label = "", pathways = pw),
               "label")
  expect_error(run_component_ora(obj, "a", "comp", label = "x", pathways = pw,
                                 minSize = 5, maxSize = 2), "maxSize")
  expect_error(suppressMessages(run_component_ora(
    obj, "a", "comp", label = "x", pathways = pw, universe = character()
  )), "universe")
  expect_error(suppressMessages(run_component_ora(
    obj, "a", "comp", label = "x", pathways = list(c("GA", "GB"))
  )), "named list")
  expect_error(run_component_ora(obj, "a", "missing", label = "x",
                                 pathways = pw))
  expect_error(get_component_ora_result(obj, "a", "comp", "missing"))
})

test_that("gene sets can be a GeneSetCollection of Ensembl IDs", {
  skip_if_not_installed("GSEABase")
  fx <- make_multi_gene_fixture()
  gsc <- GSEABase::GeneSetCollection(lapply(names(fx$pathways), function(n)
    GSEABase::GeneSet(fx$pathways[[n]], setName = n,
                      geneIdType = GSEABase::ENSEMBLIdentifier())))
  expect_identical(.as_pathway_list(gsc, genes = character()), fx$pathways)

  run <- function(pathways) {
    withr::local_seed(1)
    obj <- suppressWarnings(suppressMessages(run_fgsea(
      fx$obj, "a", label = "nodes", pathways = pathways, minSize = 2
    )))
    get_gsea_table(obj, "a", "nodes")
  }
  expect_equal(run(gsc), run(fx$pathways))

  obj <- suppressMessages(run_component_analysis(fx$obj, "a", label = "comp"))
  ora <- function(pathways) {
    get_ora_table(suppressMessages(run_component_ora(
      obj, "a", "comp", label = "ora", pathways = pathways, minSize = 2
    )), "a", "comp", "ora")
  }
  expect_equal(ora(gsc), ora(fx$pathways))
})

test_that("gene symbols in a GeneSetCollection map to network genes", {
  skip_if_not_installed("GSEABase")
  sym <- GSEABase::SymbolIdentifier()
  gsc <- GSEABase::GeneSetCollection(list(
    GSEABase::GeneSet(c("CLTCL1", "APAF1", "NOT_A_GENE"), setName = "s",
                      geneIdType = sym),
    GSEABase::GeneSet("NOT_A_GENE", setName = "empty", geneIdType = sym)
  ))
  # CLTCL1 is ENSG00000070371 and APAF1 ENSG00000120868; the unknown symbol
  # and the set left empty are dropped
  out <- .as_pathway_list(gsc, genes = c("ENSG00000070371", "ENSG00000120868"))
  expect_identical(out, list(s = c("ENSG00000070371", "ENSG00000120868")))
  # Symbols map only to the genes given, here the network's
  expect_identical(.as_pathway_list(gsc, genes = "ENSG00000070371"),
                   list(s = "ENSG00000070371"))

  entrez <- GSEABase::GeneSetCollection(list(GSEABase::GeneSet(
    "1213", setName = "e", geneIdType = GSEABase::EntrezIdentifier()
  )))
  expect_error(.as_pathway_list(entrez, genes = "x"), "symbols or Ensembl")
  expect_error(.as_pathway_list(c(a = "x"), genes = "x"), "named list")
})

test_that("GO gene sets are fetched from MSigDB with the msigdb package", {
  skip_if_not_installed("msigdb")
  skip_if_not_installed("GSEABase")
  genes <- c("ENSG00000070371", "ENSG00000120868", "ENSG00000141968")
  go <- suppressMessages(.fetch_go_pathways("GO:MF", genes))
  expect_true(length(go) > 0L)
  expect_true(all(startsWith(names(go), "GOMF_")))
  expect_true(all(unlist(go) %in% genes))
})
