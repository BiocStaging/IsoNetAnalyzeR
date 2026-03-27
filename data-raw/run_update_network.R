# ==============================================================================
# data-raw/run_update_network.R
#
# One-time script to recompute the interaction network cutoff and update the
# internal package data (R/sysdata.rda).
#
# WHEN TO RUN THIS:
#   - A new version of the raw interaction model output is available
#   - You want to evaluate a different pi_threshold
#   - You want to regenerate diagnostic plots
#
# HOW TO RUN:
#   Open an R session at the package root and run:
#     source("data-raw/run_update_network.R")
#   Or line-by-line interactively.
#
# AFTER RUNNING:
#   Rebuild the package so the new sysdata.rda is picked up:
#     devtools::load_all()   # for interactive development
#     R CMD INSTALL .        # for a clean install
#
# OUTPUT:
#   R/sysdata.rda is overwritten with three internal objects:
#     .network_data   — isoform-level igraph (vertices = transcripts,
#                       edges weighted by pi, interact flag set by pi_threshold)
#     .network_cutoff — scalar ratio cutoff (ROC q75) used by the virtual
#                       pulldown to define a significant neighbourhood overlap
#     .pi_threshold   — the pi_threshold value used to derive any_interact;
#                       stored so package functions can report it to the user
#
# INPUTS:
#   network_path  — CSV produced by the Isoform Analysis group.
#                   Expected columns: enst_1, gene_1, enst_2, gene_2, pi
#                   One row per isoform pair; pi in [0, 1].
#
#   pi_threshold  — Minimum pi value for an isoform pair to count as
#                   interacting (any_interact = 1 at the gene level).
#                   Default 0.5. Increasing this makes the network sparser
#                   and typically raises the ratio cutoff.
#
# REPRODUCIBILITY:
#   The KEGG gene sets are fetched live from MSigDB via msigdbr at runtime.
#   To pin to a specific MSigDB version, install a versioned msigdbr release
#   and pass a pre-fetched tibble via the kegg_genesets argument (see below).
# ==============================================================================

source("data-raw/update_network.R")

# ------------------------------------------------------------------------------
# Run with defaults
#   pi_threshold = 0.5  (recommended starting point)
#   plot         = FALSE (set TRUE to save diagnostic plots to data-raw/plots/)
# ------------------------------------------------------------------------------
result <- update_package_network(
  network_path = "data-raw/results_PHYSICAL_Prob_Model_16_02_26.csv",
  pi_threshold = 0.5
)

# ------------------------------------------------------------------------------
# The returned result list contains:
#
#   result$best_cutoff        — recommended ratio cutoff (ROC q75); this is
#                               what gets written to .network_cutoff
#   result$gene_pairs         — gene-level edge table with pi summaries and
#                               any_interact flag
#   result$gene_degree        — per-gene binary and pi-weighted degree
#   result$pathway_results    — per-pathway ROC/PR metrics and best cutoffs
#   result$validation         — precision/recall/FPR/F1 for all 6 candidate
#                               cutoffs (ROC and PR/F1 at q50/q75/q90)
#   result$cutoff_roc_summary — AUC-weighted quantile summary (ROC-based)
#   result$cutoff_f1_summary  — AUC-weighted quantile summary (PR/F1-based)
#   result$pi_threshold       — pi_threshold used in this run
#   result$kegg_genesets      — MSigDB KEGG gene set table used (Legacy +
#                               Medicus, Homo sapiens)
#
# Quick inspection:
# ------------------------------------------------------------------------------
cat("Best cutoff (ROC q75):", result$best_cutoff, "\n")
print(result$validation[ , c("label", "cutoff", "precision", "recall", "fpr", "f1")])

# ------------------------------------------------------------------------------
# Optional: regenerate plots from this result without rerunning the analysis
# ------------------------------------------------------------------------------
# plot_cutoff_diagnostics(result, plot_dir = "data-raw/plots")

# ------------------------------------------------------------------------------
# Optional: re-run with a different pi_threshold to compare cutoffs.
# Pass the already-fetched kegg_genesets to skip the MSigDB download.
# ------------------------------------------------------------------------------
# result_07 <- update_package_network(
#   network_path  = "data-raw/results_PHYSICAL_Prob_Model_16_02_26.csv",
#   pi_threshold  = 0.7,
#   kegg_genesets = result$kegg_genesets,   # reuse cached gene sets
#   plot          = TRUE,
#   plot_dir      = "data-raw/plots/pi_07"
# )
