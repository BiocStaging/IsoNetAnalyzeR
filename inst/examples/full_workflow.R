# ==============================================================================
#  IsoNetAnalyzeR — Full workflow example
#
#  This script walks through a complete analysis from raw DEXSeq output
#  to gene set enrichment, using a type 2 diabetes vs healthy comparison
#  as the example dataset.
#
#  Sections:
#    0. Setup
#    1. Input preparation
#    2. Session initialisation
#    3. Network diagnostics
#    4. Differential pulldown
#    5. Cross-gene exploration
#    6. Single-gene deep dive
#    7. Component analysis
#    8. ORA per component
#    9. GSEA on topology deltas
#   10. Ad hoc lookup
#   11. Session management
# ==============================================================================


# ==============================================================================
#  0. Setup
# ==============================================================================

library("IsoNetAnalyzeR")
library("future")   # parallel execution
library("dplyr")


# ==============================================================================
#  1. Input preparation
#
#  IsoNetAnalyzeR expects a long-format data frame with three columns:
#    gene_id    — Ensembl gene ID
#    isoform_id — transcript ID matching vertex names in the network
#    direction  — "up" or "down"
#
#  The code below converts a typical DEXSeq output into this format.
#  Genes that have only up- or only down-regulated isoforms in the
#  significant set are filled with the most-extreme non-significant
#  isoform from the same gene as a proxy, so every gene enters the
#  analysis with both directions represented.
# ==============================================================================

# -- Load raw DEXSeq results ---------------------------------------------------
tmp    <- readRDS("path/to/dexseq_results.RDS")
dex_df <- as.data.frame(tmp)
rm(tmp); gc()

# -- Filter to significant isoforms -------------------------------------------
sig <- dex_df[
  !is.na(dex_df$padj) &
    !is.na(dex_df$log2fold_C_B) &
    dex_df$padj < 0.05,
]

# -- Fill missing direction with best proxy isoform ---------------------------
#    Genes with only up- or only down-regulated significant isoforms are
#    supplemented with the most-extreme non-significant isoform in the
#    opposite direction so that both conditions produce a pulldown subgraph.
fill_missing_direction <- function(sig_df, all_df) {
  genes      <- unique(sig_df$groupID)
  extra_rows <- list()

  for (g in genes) {
    gene_sig <- sig_df[sig_df$groupID == g, ]
    has_up   <- any(gene_sig$log2fold_C_B > 0, na.rm = TRUE)
    has_down <- any(gene_sig$log2fold_C_B < 0, na.rm = TRUE)

    if (has_up && has_down) next   # already has both — nothing to do

    sig_isoforms <- gene_sig$featureID
    pool <- all_df[
      all_df$groupID == g &
        !all_df$featureID %in% sig_isoforms &
        !is.na(all_df$log2fold_C_B),
    ]

    if (nrow(pool) == 0) {
      warning(sprintf(
        "Gene '%s': only %s isoforms found, no proxy available.",
        g, ifelse(has_up, "up", "down")
      ))
      next
    }

    if (has_up && !has_down) {
      best <- pool[which.min(pool$log2fold_C_B), ]
      if (best$log2fold_C_B >= 0) next   # no negative proxy available
    } else {
      best <- pool[which.max(pool$log2fold_C_B), ]
      if (best$log2fold_C_B <= 0) next   # no positive proxy available
    }

    best$proxy_fill <- TRUE
    extra_rows[[length(extra_rows) + 1]] <- best
  }

  if (length(extra_rows) > 0) {
    sig_df$proxy_fill <- FALSE
    sig_df <- rbind(sig_df, do.call(rbind, extra_rows))
  } else {
    sig_df$proxy_fill <- FALSE
  }
  sig_df
}

sig_filled <- fill_missing_direction(sig_df = sig, all_df = dex_df)

# -- Convert to IsoNetAnalyzeR input format -----------------------------------
input_df <- data.frame(
  gene_id    = sig_filled$groupID,
  isoform_id = sig_filled$featureID,
  direction  = ifelse(sig_filled$log2fold_C_B > 0, "up", "down"),
  stringsAsFactors = FALSE
)

cat(sprintf(
  "Input table: %d rows | %d genes | %d unique isoforms\n",
  nrow(input_df),
  length(unique(input_df$gene_id)),
  length(unique(input_df$isoform_id))
))


# ==============================================================================
#  2. Session initialisation
# ==============================================================================

obj <- IsoNetAnalyzeR("T2D_vs_healthy")

# Attach analysis metadata
obj <- set_metadata(obj, "author",    "Your Name")
obj <- set_metadata(obj, "condition", "T2D vs healthy")
obj <- set_metadata(obj, "notes",     "Full cohort, padj < 0.05")

print_metadata(obj)


# ==============================================================================
#  3. Network diagnostics
#
#  Run before the analysis to understand the network structure. Identifies
#  sticky nodes (non-specific high-degree binders) that should be interpreted
#  with caution when they appear in pulldown results.
# ==============================================================================

obj <- run_network_diagnostics(
  obj,
  sticky_method = "quantile",
  sticky_cutoff = 0.95        # top 5% of nodes by degree are flagged sticky
)

diag <- get_network_diagnostics(obj)
diag   # show() summary: 19 global statistics

# Six-panel diagnostic figure
plot_network_diagnostics(obj)

# Inspect sticky nodes before proceeding
head(diag@sticky_isoform_table, 20)
head(diag@sticky_gene_table,    20)


# ==============================================================================
#  4. Differential pulldown
#
#  Core analysis: for each gene, runs virtual pulldowns for the up- and
#  down-regulated isoform sets separately and compares their subgraph topology.
# ==============================================================================

plan(multisession, workers = 4)   # adjust workers to your machine

obj <- run_differential_pulldown(
  obj     = obj,
  label   = "T2D_pulldown",
  input   = input_df,
  order   = 1,          # BFS hops: 1 = direct interactors only
  cutoff  = NULL,       # use stored network_cutoff from session metadata
  parallel = TRUE,
  confirm  = FALSE,     # suppress interactive prompt in scripts
  verbose  = FALSE
)

plan(sequential)


# ==============================================================================
#  5. Cross-gene exploration
#
#  The 35-column cross-gene summary has one row per gene. Scatter and boxplot
#  overviews give a quick visual read on which topology metrics differ most
#  between the up and down conditions across all genes.
# ==============================================================================

# -- Summary table -------------------------------------------------------------
cgs <- get_cross_gene_summary(obj, "T2D_pulldown")
head(cgs)

# Most rewired genes by edge gain/loss
cgs |>
  mutate(delta_edges = up_edges - down_edges) |>
  arrange(desc(abs(delta_edges))) |>
  select(gene_id, up_edges, down_edges, n_edges_gained, n_edges_lost, jaccard) |>
  head()

# Genes with very different interaction partners (low Jaccard)
low_overlap <- cgs |>
  filter(!is.na(jaccard), jaccard < 0.1)

nrow(low_overlap)

# -- Scatter overview (six panels: up vs down for key topology metrics) --------
plot_cross_gene_overview(obj, "T2D_pulldown", log_counts = TRUE)

# Custom scatter: any two numeric columns
plot_cross_gene_scatter(
  obj, "T2D_pulldown",
  x_col   = "up_edges",
  y_col   = "down_edges",
  x_axis_log = TRUE,
  y_axis_log = TRUE,
  label_n = 5,
  label   = "delta",
  title   = "Edge count: up vs down"
)

# -- Boxplot overview (topology deltas) ----------------------------------------
plot_cross_gene_boxplot_overview(obj, "T2D_pulldown", log_ratio = TRUE)

# Custom boxplot: single column pair
plot_cross_gene_boxplot(
  obj, "T2D_pulldown",
  x_col   = "up_edges",
  y_col   = "down_edges",
  log_ratio = FALSE,
  label_n = 5
)


# ==============================================================================
#  6. Single-gene deep dive
#
#  Retrieve per-isoform statistics and visualise the interaction subgraphs
#  for individual genes of interest.
# ==============================================================================

gene_of_interest <- "ENSG00000120333"

# Per-isoform statistics for both conditions
iso <- get_isoform_stats(obj, "T2D_pulldown", gene_of_interest)
head(iso)

# Seed isoforms only (supplied as input)
iso |> filter(is_seed)

# Up-condition isoforms ranked by degree
iso |>
  filter(direction == "up") |>
  arrange(desc(degree))

# Retrieve the full GeneDiffPulldown object
gr <- get_gene_result(obj, "T2D_pulldown", gene_of_interest)
gr   # show() summary

# Visualise
plot_network(gr)                                    # side-by-side (default)
plot_network(gr, direction = "overlay")             # up/down/shared membership
plot_network(gr, direction = "overlay",
             group_by_gene = TRUE)                  # gene fill + membership border
plot_network(gr, direction = "up")                  # up condition only
plot_network(gr, level = "gene")                    # gene-level contraction


# ==============================================================================
#  7. Component analysis
#
#  Merges all per-gene subgraphs into a cross-gene network and decomposes it
#  into connected components. Run at isoform level (both directions) for
#  rewiring detection, and additionally for up and down separately to
#  characterise condition-specific communities.
# ==============================================================================

# -- Isoform level: both directions (primary analysis) -------------------------
obj <- run_component_analysis(
  obj,
  analysis            = "T2D_pulldown",
  label               = "iso_comp_both",
  direction           = "both",
  level               = "isoform",
  min_component_nodes = 5,    # drop tiny components
  min_shared_nodes    = 2     # require at least 2 cross-gene shared nodes
)

# -- Isoform level: condition-specific -----------------------------------------
obj <- run_component_analysis(
  obj,
  analysis            = "T2D_pulldown",
  label               = "iso_comp_up",
  direction           = "up",
  level               = "isoform",
  min_component_nodes = 5,
  min_shared_nodes    = 2
)

obj <- run_component_analysis(
  obj,
  analysis            = "T2D_pulldown",
  label               = "iso_comp_down",
  direction           = "down",
  level               = "isoform",
  min_component_nodes = 5,
  min_shared_nodes    = 2
)

# -- Gene level: both directions -----------------------------------------------
obj <- run_component_analysis(
  obj,
  analysis            = "T2D_pulldown",
  label               = "gene_comp_both",
  direction           = "both",
  level               = "gene",
  min_component_nodes = 5,
  min_shared_nodes    = 2
)

# -- Explore component overview ------------------------------------------------
overview <- get_component_overview(obj, "T2D_pulldown", "iso_comp_both")
head(overview)

# Rewiring candidates: components with both up and down seeds
rewiring <- overview |>
  filter(mixed) |>
  select(component_id, n_nodes, n_seeds_up, n_seeds_down, topology_class)

print(rewiring)

# Hub-spoke components with shared nodes (potential master regulators)
overview |>
  filter(!is.na(topology_class),
         topology_class == "hub-spoke",
         n_shared_nodes > 0) |>
  arrange(desc(degree_centralization))

# -- Per-component node table --------------------------------------------------
# Component 1 is always the largest (ordered by node count)
comp1 <- get_component(obj, "T2D_pulldown", "iso_comp_both", component_id = 1)
head(comp1)

# Nodes bridging across gene communities (high normalised betweenness)
comp1 |> filter(!is.na(betweenness_norm), betweenness_norm > 0.5)

# -- Component visualisation ---------------------------------------------------
# Cross-gene network: all components coloured by source gene
result <- get_analysis(obj, "T2D_pulldown")
plot_network(result, direction = "both", level = "isoform")
plot_network(result, direction = "both", level = "gene")

# Focus on shared hubs (hub_degree mode)
plot_network(result, direction = "both", level = "isoform",
             colour_by        = "hub_degree",
             min_shared_nodes = 2,
             min_component_nodes = 5)

# Single component plots
plot_component(obj, "T2D_pulldown", "iso_comp_both", component_id = 1,
               colour_by = "direction")                          # seeds coloured

plot_component(obj, "T2D_pulldown", "iso_comp_both", component_id = 1,
               colour_by = "direction", colour_neighbours = TRUE) # territory map

plot_component(obj, "T2D_pulldown", "iso_comp_both", component_id = 1,
               colour_by = "hub_degree")                         # shared hub focus

plot_component(obj, "T2D_pulldown", "gene_comp_both", component_id = 1,
               colour_by = "gene")                               # gene-level view


# ==============================================================================
#  8. ORA per component
#
#  Tests GO terms against each component's gene set using the full network
#  as the background universe. Run for each component analysis label.
# ==============================================================================

# -- GO:BP on both-direction components ----------------------------------------
obj <- run_component_ora(
  obj,
  analysis           = "T2D_pulldown",
  component_analysis = "iso_comp_both",
  label              = "ora_bp",
  go_category        = "GO:BP",
  nodes              = "all"    # or "seeds" to use only seed genes as foreground
)

# -- GO:BP on condition-specific components ------------------------------------
obj <- run_component_ora(
  obj,
  analysis           = "T2D_pulldown",
  component_analysis = "iso_comp_up",
  label              = "ora_bp_up",
  go_category        = "GO:BP"
)

obj <- run_component_ora(
  obj,
  analysis           = "T2D_pulldown",
  component_analysis = "iso_comp_down",
  label              = "ora_bp_down",
  go_category        = "GO:BP"
)

# -- Retrieve results ----------------------------------------------------------
# All significant terms across all components
get_ora_table(obj, "T2D_pulldown", "iso_comp_both", "ora_bp",
              padj_cutoff = 0.05)

# Top 3 GO terms per component
get_ora_table(obj, "T2D_pulldown", "iso_comp_both", "ora_bp",
              top_n = 3)

# Top 5 significant terms for component 1 only
get_ora_table(obj, "T2D_pulldown", "iso_comp_both", "ora_bp",
              component_id = 1, padj_cutoff = 0.05, top_n = 5)

# Full S4 result object (inspect params, universe, foreground genes)
ora_obj <- get_component_ora_result(obj, "T2D_pulldown", "iso_comp_both", "ora_bp")
ora_obj
ora_obj@params
ora_obj@universe


# ==============================================================================
#  9. GSEA on topology deltas
#
#  Ranks all genes by a topology delta score and tests GO pathways.
#  Use a built-in preset or supply a custom ranking vector.
# ==============================================================================

# -- Built-in preset: delta_nodes (default) ------------------------------------
obj <- run_fgsea(
  obj,
  analysis = "T2D_pulldown",
  label    = "gsea_nodes",
  rank_by  = "delta_nodes"   # up_nodes - down_nodes
)

# -- Built-in preset: edge rewiring balance ------------------------------------
obj <- run_fgsea(
  obj,
  analysis = "T2D_pulldown",
  label    = "gsea_rewiring",
  rank_by  = "edges_gained_vs_lost"
)

# -- Custom ranking: weighted combination of node and edge delta ---------------
my_stats <- get_cross_gene_summary(obj, "T2D_pulldown") |>
  mutate(stat = (up_nodes - down_nodes) * 0.7 +
           (up_edges - down_edges) * 0.3) |>
  select(gene_id, stat) |>
  tibble::deframe()

obj <- run_fgsea(
  obj,
  analysis = "T2D_pulldown",
  label    = "gsea_weighted",
  stats    = my_stats   # rank_by is ignored when stats is supplied
)

# -- Retrieve results ----------------------------------------------------------
# Significant pathways
sig_gsea <- get_gsea_table(obj, "T2D_pulldown", "gsea_nodes",
                           padj_cutoff = 0.05) |>
  arrange(NES) |>
  select(pathway, NES, padj, size)

# Genes with the largest up-condition network expansion
ranks <- get_gsea_stats(obj, "T2D_pulldown", "gsea_nodes")
head(ranks, 10)
tail(ranks, 10)   # genes with the largest down-condition network expansion

# Full S4 result object
gsea_obj <- get_gsea_result(obj, "T2D_pulldown", "gsea_nodes")
gsea_obj
gsea_obj@rank_by
gsea_obj@params


# ==============================================================================
#  10. Ad hoc lookup
#
#  Quick single-seed pulldown for any isoform or gene of interest, without
#  needing a full differential input table.
# ==============================================================================

# Single isoform seed
obj <- lookup_virtual_pulldown(
  obj,
  label = "BRCA1_isoform",
  seed  = "ENST00000357654",
  order = 1,
  level = "isoform"
)

# Multiple seeds — all treated as order-0 anchors in one pulldown
obj <- lookup_virtual_pulldown(
  obj,
  label = "BRCA1_combined",
  seed  = c("ENST00000357654", "ENST00000471181"),
  order = 1
)

# Gene-level seed — expands from all isoforms of the gene
obj <- lookup_virtual_pulldown(
  obj,
  label = "BRCA1_gene",
  seed  = "ENSG00000012048",
  order = 1,
  level = "gene"
)

# Retrieve results
get_lookup_stats(obj, "BRCA1_isoform")
get_lookup_subgraph(obj, "BRCA1_isoform")
plot_network(get_lookup(obj, "BRCA1_isoform"))


# ==============================================================================
#  11. Session management
# ==============================================================================

# List all stored analyses and lookups
list_analyses(obj)
list_lookups(obj)

# Save the full session object for later
saveRDS(obj, "isonetanalyzer_T2D_session.rds")

# Reload
obj <- readRDS("isonetanalyzer_T2D_session.rds")
