# Analysis-agnostic IsoNetAnalyzeR workflow using synthetic data.
# Seed derivation is intentionally outside this script and package.

library(IsoNetAnalyzeR)

vertices <- data.frame(
  name = c("tx_A_up", "tx_A_down", "tx_B", "tx_C_isolated"),
  gene = c("gene_A", "gene_A", "gene_B", "gene_C")
)

edges <- data.frame(
  from = c("tx_A_up", "tx_A_down"),
  to = c("tx_B", "tx_B"),
  interact = c(1L, 1L),
  pi = c(0.90, 0.75),
  pi_score = c(0.90, 0.75)
)

vertices$transcript_biotype <- "protein_coding"
network <- igraph::graph_from_data_frame(
  edges, directed = FALSE, vertices = vertices
)

obj <- IsoNetAnalyzeR(
  "synthetic-workflow",
  network = network,
  metadata = list(network_cutoff = 0)
)

# The caller supplies this validated three-column contract. Versioned ENSTs and
# graph-compatible protein-identical duplicates resolve automatically; use
# resolve_duplicates = FALSE only to reproduce the exact-ID v0.4.0 behaviour.
seed_table <- data.frame(
  gene_id = c("gene_A", "gene_A", "gene_C"),
  isoform_id = c("tx_A_up", "tx_A_down", "tx_C_isolated"),
  direction = c("up", "down", "up")
)

obj <- run_differential_pulldown(
  obj, "synthetic", seed_table,
  order = 1, cutoff = 0, confirm = FALSE, verbose = FALSE
)

coverage <- obj@metadata$analysis_network_coverage$synthetic
coverage$overall
coverage$by_direction
obj@metadata$analysis_duplicate_resolution$synthetic$summary

summary <- get_cross_gene_summary(obj, "synthetic")
summary

gene_A <- get_gene_result(obj, "synthetic", "gene_A")
gene_A@graph_summary

p <- plot_network(gene_A, layout = "circle")
# ggplot2::ggsave("gene_A.pdf", p, width = 7, height = 5)

obj <- run_component_analysis(
  obj, "synthetic", "components", direction = "both", level = "gene"
)
get_component_overview(obj, "synthetic", "components")

obj <- lookup_virtual_pulldown(
  obj, "isolated_lookup", "tx_C_isolated", order = 1, cutoff = 0
)
get_lookup_stats(obj, "isolated_lookup")

# Count aggregation belongs before differential fitting. The packaged map is
# transcriptome-wide, including mappings outside the embedded network.
# aggregated <- aggregate_duplicate_counts(count_matrix, feature_metadata)
# count_matrix <- aggregated$object
# feature_metadata <- aggregated$feature_data

# Reproducible enrichment requires a frozen, versioned pathway collection.
# pathways <- readRDS("frozen_pathways.rds")
# set.seed(2026)
# obj <- run_fgsea(
#   obj, "synthetic", "frozen_gsea", pathways = pathways,
#   rank_by = "delta_nodes", minSize = 10, maxSize = 500, nperm = 10000
# )

# Save the complete session for reproducibility.
# saveRDS(obj, "synthetic_isonet_session.rds")
