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

  out <- cli::cli_fmt(show(res))
  expect_true(any(grepl("Rank by *: delta_nodes", out)))
})

test_that("run_fgsea uses fgseaMultilevel unless nperm is given", {
  fx <- make_multi_gene_fixture()
  withr::local_seed(1)
  warnings <- character()
  obj <- withCallingHandlers(
    suppressMessages(run_fgsea(fx$obj, "a", label = "multilevel",
                               pathways = fx$pathways, minSize = 2)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("fgseaSimple", warnings)))
  params <- get_gsea_result(obj, "a", "multilevel")@params
  expect_null(params$nperm)
  expect_identical(params$method, "fgseaMultilevel")

  obj <- suppressWarnings(suppressMessages(run_fgsea(
    obj, "a", label = "simple", pathways = fx$pathways, minSize = 2,
    nperm = 1000
  )))
  res <- get_gsea_result(obj, "a", "simple")
  expect_identical(res@params$method, "fgseaSimple")
  # With nperm = 1000 no p-value can fall below 1 / 1001
  expect_message(none <- get_gsea_table(obj, "a", "simple", padj_cutoff = 1e-6),
                 "No pathways")
  expect_identical(nrow(none), 0L)
  expect_true(any(grepl("fgseaSimple, nperm = 1000", cli::cli_fmt(show(res)))))
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
  expect_equal(ge$pval, stats::phyper(0, 3, 11, 1, lower.tail = FALSE))
  expect_equal(ge$overlap, 1)

  # Component 1 holds 11 genes, including all four mixed_set genes.
  mixed <- rt[rt$component_id == 1 & rt$pathway == "mixed_set", ]
  expect_equal(mixed$pval, stats::phyper(3, 4, 10, 11, lower.tail = FALSE))
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
  genes <- unique(igraph::V(get_network(fx$obj))$gene)
  expect_identical(.as_pathway_list(gsc, genes = genes), fx$pathways)

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
  skip_if_no_ensembl()
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
  # Versioned network gene IDs match and are returned as the network has them
  expect_identical(
    .as_pathway_list(gsc, genes = c("ENSG00000070371.16",
                                    "ENSG00000120868.13")),
    list(s = c("ENSG00000070371.16", "ENSG00000120868.13"))
  )

  entrez <- GSEABase::GeneSetCollection(list(GSEABase::GeneSet(
    "1213", setName = "e", geneIdType = GSEABase::EntrezIdentifier()
  )))
  expect_error(.as_pathway_list(entrez, genes = "x"), "symbols or Ensembl")
  expect_error(.as_pathway_list(c(a = "x"), genes = "x"), "named list")
})

test_that("collections without an identifier type are read as IDs or symbols", {
  skip_if_not_installed("GSEABase")
  expect_error(.as_pathway_list(GSEABase::GeneSetCollection(list()),
                                genes = character()), "no gene sets")
  # GeneSet() and GSEABase::getGmt() give NullIdentifier gene sets. IDs that
  # all start with ENSG are Ensembl gene IDs
  cltcl1_apaf1 <- c("ENSG00000070371", "ENSG00000120868")
  ensembl <- GSEABase::GeneSetCollection(list(
    GSEABase::GeneSet(cltcl1_apaf1, setName = "e")))
  expect_identical(.as_pathway_list(ensembl, genes = cltcl1_apaf1),
                   list(e = cltcl1_apaf1))
  # A collection with no network gene is an error, not zero gene sets
  expect_error(.as_pathway_list(ensembl, genes = "ENSG00000141968"),
               "No gene of the gene sets is a network gene")
  # Other IDs are gene symbols: CLTCL1 and APAF1
  skip_if_no_ensembl()
  symbols <- GSEABase::GeneSetCollection(list(
    GSEABase::GeneSet(c("CLTCL1", "APAF1"), setName = "s")))
  expect_identical(.as_pathway_list(symbols, genes = cltcl1_apaf1),
                   list(s = cltcl1_apaf1))
  expect_error(.as_pathway_list(symbols, genes = "ENSG00000141968"),
               "No gene of the gene sets is a network gene")
})

test_that("the MSigDB GO gene sets are loaded and converted once", {
  skip_if_not_installed("msigdb")
  skip_if_not_installed("GSEABase")
  cached <- mget(c("msigdb_go", "go_pathways"), envir = .annotation_cache,
                 ifnotfound = list(NULL))
  withr::defer(list2env(cached, envir = .annotation_cache))
  # A stand-in for the loaded collection: one GO:MF set holding CLTCL1,
  # whose gene is ENSG00000070371
  local_mocked_bindings(isonet_gene_names = function(edb = NULL)
    S4Vectors::DataFrame(gene_id = "ENSG00000070371", gene_name = "CLTCL1"))
  .annotation_cache$go_pathways <- NULL
  .annotation_cache$msigdb_go <- list(
    sets = GSEABase::GeneSetCollection(list(GSEABase::GeneSet(
      "CLTCL1", setName = "GOMF_CACHED",
      geneIdType = GSEABase::SymbolIdentifier()))),
    subcategory = "GO:MF")
  expect_identical(
    suppressMessages(.fetch_go_pathways("GO:MF", "ENSG00000070371")),
    list(GOMF_CACHED = "ENSG00000070371")
  )
  # For the same genes the converted sets are reused, not the collection
  .annotation_cache$msigdb_go$sets <- GSEABase::GeneSetCollection(list())
  expect_identical(
    suppressMessages(.fetch_go_pathways("GO:MF", "ENSG00000070371")),
    list(GOMF_CACHED = "ENSG00000070371")
  )
})

test_that("GO gene sets without a network gene give a GO error", {
  skip_if_not_installed("msigdb")
  skip_if_not_installed("GSEABase")
  cached <- mget(c("msigdb_go", "go_pathways"), envir = .annotation_cache,
                 ifnotfound = list(NULL))
  withr::defer(list2env(cached, envir = .annotation_cache))
  local_mocked_bindings(isonet_gene_names = function(edb = NULL)
    S4Vectors::DataFrame(gene_id = "ENSG00000070371", gene_name = "CLTCL1"))
  # The only GO:MF set holds a symbol that is not a network gene
  .annotation_cache$go_pathways <- NULL
  .annotation_cache$msigdb_go <- list(
    sets = GSEABase::GeneSetCollection(list(GSEABase::GeneSet(
      "NOT_A_GENE", setName = "GOMF_OTHER",
      geneIdType = GSEABase::SymbolIdentifier()))),
    subcategory = "GO:MF")
  expect_error(
    suppressMessages(.fetch_go_pathways("GO:MF", "ENSG00000070371")),
    "No GO:MF gene set of MSigDB 2023.1 contains a network gene",
    fixed = TRUE)
})

test_that("named gene lists are matched to versioned network gene IDs", {
  versioned <- c("ENSG00000070371.16", "ENSG00000120868.13")
  sets <- list(s = c("ENSG00000070371", "ENSG00000120868", "ENSG_OTHER"))
  # Matched without versions and returned as the network has them; IDs
  # that are not network genes are kept as given
  expect_identical(.as_pathway_list(sets, genes = versioned),
                   list(s = c(versioned, "ENSG_OTHER")))
  # A list with no matching gene is an error that does not advise symbols
  msg <- tryCatch(.as_pathway_list(list(s = "ENSG00000141968"),
                                   genes = versioned),
                  error = conditionMessage)
  expect_match(msg, "No gene of the named list is a gene of the analysis")
  expect_match(msg, "For gene symbols, pass a GeneSetCollection")
})

test_that("gene-set IDs that are target genes stay themselves", {
  # Clone-based gene names are no Ensembl IDs: both stay as given
  clones <- c("AC000061.1", "AC000061.2")
  expect_identical(.as_pathway_list(list(s = clones), genes = clones),
                   list(s = clones))
  # A network with both copies of the PAR gene ENSG00000002586 (CD99): the
  # exact chrY copy stays, the unversioned ID matches both copies
  copies <- c("ENSG00000002586.20", "ENSG00000002586.20_PAR_Y")
  expect_identical(.as_pathway_list(list(y = copies[2]), genes = copies),
                   list(y = copies[2]))
  expect_identical(.as_pathway_list(list(u = "ENSG00000002586"),
                                    genes = copies),
                   list(u = copies))
  # The symbol CD99 maps to both copies too
  skip_if_not_installed("GSEABase")
  local_mocked_bindings(isonet_gene_names = function(edb = NULL)
    S4Vectors::DataFrame(gene_id = "ENSG00000002586", gene_name = "CD99"))
  cd99 <- GSEABase::GeneSetCollection(list(GSEABase::GeneSet(
    "CD99", setName = "s", geneIdType = GSEABase::SymbolIdentifier())))
  expect_identical(.as_pathway_list(cd99, genes = copies), list(s = copies))
})

test_that("untyped collections of other species' Ensembl IDs stay IDs", {
  skip_if_not_installed("GSEABase")
  # ENSMUSG IDs, with a blank entry, are Ensembl gene IDs, not symbols
  mouse <- c("ENSMUSG00000000001", "ENSMUSG00000000028")
  gsc <- GSEABase::GeneSetCollection(list(GSEABase::GeneSet(
    c(mouse, ""), setName = "m")))
  expect_identical(.as_pathway_list(gsc, genes = mouse),
                   list(m = c(mouse, "")))
})

test_that("run_fgsea matches gene sets to the ranked genes' IDs", {
  fx <- make_multi_gene_fixture()
  withr::local_seed(1)
  # Statistics keyed by symbols, and a named list of the same symbols
  stats <- c(SYM_A = 3, SYM_B = 2, SYM_C = 1, SYM_D = -1, SYM_E = -2,
             SYM_F = -3)
  sets <- list(top = c("SYM_A", "SYM_B", "SYM_C"),
               bottom = c("SYM_E", "SYM_F"))
  obj <- suppressWarnings(suppressMessages(run_fgsea(
    fx$obj, "a", label = "symbols", stats = stats, pathways = sets,
    minSize = 2)))
  tbl <- get_gsea_table(obj, "a", "symbols")
  # Every member is ranked: sizes 3 and 2; "top" holds the three highest
  expect_equal(tbl$size[match(c("top", "bottom"), tbl$pathway)], c(3, 2))
  expect_gt(tbl$NES[tbl$pathway == "top"], 0)
})

test_that("gene sets match GENCODE PAR_Y network gene IDs", {
  skip_if_not_installed("GSEABase")
  # ENSG00000002586 is CD99; GENCODE marks its Y copy with _PAR_Y
  par_y <- "ENSG00000002586.20_PAR_Y"
  gsc <- GSEABase::GeneSetCollection(list(GSEABase::GeneSet(
    "ENSG00000002586", setName = "e",
    geneIdType = GSEABase::ENSEMBLIdentifier())))
  expect_identical(.as_pathway_list(gsc, genes = par_y), list(e = par_y))
})

test_that("Ensembl gene sets match versioned network gene IDs", {
  skip_if_not_installed("GSEABase")
  versioned <- c("ENSG00000070371.16", "ENSG00000120868.13")
  for (type in c(GSEABase::ENSEMBLIdentifier(), GSEABase::NullIdentifier())) {
    gsc <- GSEABase::GeneSetCollection(list(GSEABase::GeneSet(
      c("ENSG00000070371", "ENSG00000120868"), setName = "e",
      geneIdType = type)))
    # Matched without the versions and returned as the network has them
    expect_identical(.as_pathway_list(gsc, genes = versioned),
                     list(e = versioned))
  }
})

test_that("component ORA resolves its analysis once", {
  fx <- make_multi_gene_fixture()
  obj <- suppressMessages(run_component_analysis(fx$obj, "a", label = "comp"))
  obj@analyses$a@cross_gene_summary <- obj@analyses$a@cross_gene_summary[0, ]
  # One warning for the empty cross-gene summary
  warnings <- character()
  withCallingHandlers(
    suppressMessages(run_component_ora(obj, "a", "comp", label = "ora",
                                       pathways = fx$pathways, minSize = 2)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_identical(sum(grepl("empty cross-gene summary", warnings)), 1L)
})

test_that("GO gene sets are fetched from MSigDB with the msigdb package", {
  skip_if_no_ensembl()
  skip_if_no_msigdb()
  genes <- c("ENSG00000070371", "ENSG00000120868", "ENSG00000141968")
  go <- suppressMessages(.fetch_go_pathways("GO:MF", genes))
  expect_true(length(go) > 0L)
  expect_true(all(startsWith(names(go), "GOMF_")))
  expect_true(all(unlist(go) %in% genes))
})

test_that("factor-valued gene sets keep their genes", {
  fx <- make_multi_gene_fixture()
  genes <- c("GA", "GB", "GC", "GD", "GE", "GF", paste0("N", 1:8))
  as_factors <- lapply(fx$pathways, factor)
  # A factor's labels are its genes; its level numbers are not
  expect_identical(.network_gene_sets(as_factors, genes, symbols = FALSE),
                   .network_gene_sets(fx$pathways, genes, symbols = FALSE))
  expect_identical(.network_gene_sets(as_factors, genes, symbols = FALSE),
                   fx$pathways)
})
