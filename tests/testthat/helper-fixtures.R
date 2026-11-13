# Tests that need the Ensembl 107 annotation (AnnotationHub AH104864) or the
# MSigDB gene sets (ExperimentHub) are skipped when these cannot be loaded,
# for example offline without a hub cache. Each hub is tried once per run.
hub_loaded <- new.env(parent = emptyenv())

skip_if_no_ensembl <- function() {
  if (is.null(hub_loaded$ensembl))
    hub_loaded$ensembl <- tryCatch({
      suppressMessages(.isonet_ensdb())
      TRUE
    }, error = function(e) FALSE)
  if (!hub_loaded$ensembl)
    skip("The Ensembl 107 annotation (AnnotationHub) cannot be loaded.")
}

skip_if_no_msigdb <- function() {
  skip_if_not_installed("msigdb")
  skip_if_not_installed("GSEABase")
  if (is.null(hub_loaded$msigdb))
    hub_loaded$msigdb <- tryCatch({
      suppressMessages(.msigdb_go_sets())
      TRUE
    }, error = function(e) FALSE)
  if (!hub_loaded$msigdb)
    skip("The MSigDB gene sets (ExperimentHub) cannot be loaded.")
}

# Until the calling test ends: AnnotationHub in local mode with an empty
# cache and no annotation loaded in the session, as offline without a cache
local_offline_annotation <- function(env = parent.frame()) {
  options <- list(CACHE = AnnotationHub::getAnnotationHubOption("CACHE"),
                  LOCAL = AnnotationHub::getAnnotationHubOption("LOCAL"))
  cached <- mget(c("ensdb", "biotypes", "gene_names", "symbols_failed"),
                 envir = .annotation_cache, ifnotfound = list(NULL))
  withr::defer({
    for (o in names(options))
      AnnotationHub::setAnnotationHubOption(o, options[[o]])
    list2env(cached, envir = .annotation_cache)
  }, envir = env)
  AnnotationHub::setAnnotationHubOption(
    "CACHE", withr::local_tempdir(.local_envir = env))
  AnnotationHub::setAnnotationHubOption("LOCAL", TRUE)
  list2env(list(ensdb = NULL, biotypes = NULL, gene_names = NULL),
           envir = .annotation_cache)
}

# Gene G1's up and down seeds share one neighbour (tx_neighbour), so each
# neighbourhood has 2 nodes and 1 edge. G3's seeds are in the network without
# edges (one-node neighbourhoods); G_missing's seeds are not in the network
# (no neighbourhood); tx_other is an unconnected non-seed transcript.
make_isolated_seed_fixture <- function() {
  vertices <- data.frame(
    name = c("tx_up", "tx_down", "tx_neighbour", "tx_isolated_up",
             "tx_isolated_down", "tx_other"),
    gene = c("G1", "G1", "G2", "G3", "G3", "G4"),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    from = c("tx_up", "tx_down"),
    to = c("tx_neighbour", "tx_neighbour"),
    interact = c(1L, 1L),
    pi = c(0.9, 0.8),
    stringsAsFactors = FALSE
  )
  vertices$transcript_biotype <- "protein_coding"
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
      "tx_up", "tx_down", "tx_isolated_up", "tx_isolated_down",
      "tx_absent_up", "tx_absent_down"
    ),
    direction = c("up", "down", "up", "down", "up", "down"),
    stringsAsFactors = FALSE
  )
  list(obj = obj, input = input)
}

# Six seed genes with designed neighbourhood sizes. Up-minus-down node counts
# are GA +2, GB -2, GC +1, GD 0, GE +2 and GF -1. GE's down isoform has no
# predicted interaction. Shared neighbours join genes into three components.
make_multi_gene_fixture <- function() {
  genes <- c("GA", "GB", "GC", "GD", "GE", "GF")
  neighbours <- paste0("N", 1:8)
  vertices <- data.frame(
    name = c(paste0(genes, "_up"), paste0(genes, "_down"),
             paste0("tx_", neighbours)),
    gene = c(genes, genes, neighbours),
    transcript_biotype = "protein_coding",
    stringsAsFactors = FALSE
  )
  pairs <- rbind(
    c("GA_up", "tx_N1"), c("GA_up", "tx_N2"), c("GA_up", "tx_N3"),
    c("GA_down", "tx_N1"),
    c("GB_up", "tx_N4"),
    c("GB_down", "tx_N4"), c("GB_down", "tx_N5"), c("GB_down", "tx_N6"),
    c("GC_up", "tx_N7"), c("GC_up", "tx_N8"), c("GC_down", "tx_N7"),
    c("GD_up", "tx_N2"), c("GD_down", "tx_N5"),
    c("GE_up", "tx_N3"), c("GE_up", "tx_N6"),
    c("GF_up", "tx_N8"), c("GF_down", "tx_N1"), c("GF_down", "tx_N2")
  )
  edges <- data.frame(from = pairs[, 1], to = pairs[, 2], interact = 1L,
                      pi = 0.8, stringsAsFactors = FALSE)
  network <- igraph::graph_from_data_frame(edges, directed = FALSE,
                                           vertices = vertices)
  obj <- suppressMessages(IsoNetAnalyzeR(
    "multi-gene", network = network, metadata = list(network_cutoff = 0)
  ))
  seeds <- data.frame(
    gene_id = rep(genes, each = 2),
    isoform_id = as.vector(rbind(paste0(genes, "_up"),
                                 paste0(genes, "_down"))),
    direction = rep(c("up", "down"), length(genes)),
    stringsAsFactors = FALSE
  )
  obj <- suppressMessages(run_differential_pulldown(
    obj, "a", seeds, order = 1, confirm = FALSE, verbose = FALSE
  ))
  pathways <- list(
    up_set = c("GA", "GC", "GE"),
    down_set = c("GB", "GF"),
    neighbour_set = c("N1", "N2", "N3"),
    mixed_set = c("GA", "GB", "GD", "N4")
  )
  list(obj = obj, pathways = pathways, seeds = seeds)
}
