# ==============================================================================
# data-raw/update_network.R
#
# Top-level orchestrator: sources all data-raw helpers, runs the full network
# update pipeline, and writes new internal package data to R/sysdata.rda.
#
# WHEN TO RUN THIS:
#   - A new version of the raw interaction model output is available.
#   - You want to evaluate a different pi_threshold.
#   - You want to regenerate diagnostic plots.
#
# HOW TO RUN:
#   Open an R session at the package root and run:
#     source("data-raw/update_network.R")
#   Or step through the sections below interactively.
#
# AFTER RUNNING:
#   Rebuild the package so the new sysdata.rda is picked up:
#     devtools::load_all()   # for interactive development
#     R CMD INSTALL .        # for a clean install
#
# OUTPUT:
#   R/sysdata.rda is overwritten with three internal objects:
#     .network_data   — isoform-level igraph (transcripts as vertices,
#                       edges weighted by pi, interact flag set by pi_threshold)
#     .network_cutoff — scalar ratio cutoff (PR/F1-based q50) used by the
#                       virtual pulldown to define significant neighbourhood overlap
#     .pi_threshold   — the pi_threshold value used in this run
#
#   data-raw/logs/update_network_<timestamp>.log — full console log
#
# SOURCE FILES (all in data-raw/):
#   logger.R           — tee-logging utility
#   network_helpers.R  — pure data-manipulation helpers
#   network_analysis.R — compute_network_cutoff()
#   network_plots.R    — plot_cutoff_diagnostics(), plot_roc_curves(),
#                        plot_pr_curves()
#   kegg_plots.R       — plot_kegg_characterisation(),
#                        plot_gene_pathway_membership()
#
# INPUTS:
#   network_path  — CSV produced by the Isoform Analysis group.
#                   Expected columns: enst_1, gene_1, enst_2, gene_2, pi
#                   One row per isoform pair; pi in [0, 1].
#
#   pi_threshold  — Minimum pi value for an isoform pair to count as
#                   interacting. Default 0.5.
#
# REPRODUCIBILITY:
#   KEGG gene sets are fetched live from MSigDB via msigdbr at runtime.
#   To pin a specific MSigDB version, pass a pre-fetched tibble via the
#   kegg_genesets argument (see commented examples below).
# ==============================================================================


# ==============================================================================
# SECTION 1: Load packages and source helpers
# ==============================================================================

library("dplyr")
library("tidyr")
library("purrr")
library("msigdbr")
library("ggplot2")
library("ggrepel")
library("cowplot")
library("future")
library("future.apply")
library("tibble")
library("progressr")
library("scales")
library("stringr")
library("forcats")
library("usethis")
library("igraph")

# Source the data-raw helper modules in dependency order
source("data-raw/logger.R")
source("data-raw/network_helpers.R")
source("data-raw/network_analysis.R")
source("data-raw/network_plots.R")
source("data-raw/kegg_plots.R")

# Progress bars: with_progress() in each heavy step handles this locally, so
# a global handler is optional. Set one here if you want a specific style.
# progressr::handlers("cli")

# ==============================================================================
# update_package_network()
# ==============================================================================

#' Run the cutoff analysis and write results to internal package data
#'
#' Reads the raw isoform interaction CSV once, builds the isoform-level igraph,
#' runs \code{compute_network_cutoff()}, then writes three internal objects to
#' \code{R/sysdata.rda}.  All console output is simultaneously written to a
#' timestamped log file in \code{log_dir}.
#'
#' The parallel plan is managed inside this function via \code{n_workers}.
#' Workers are reused across repeated calls when \code{n_workers} is unchanged,
#' and respawned only when a different count is requested.
#'
#' @param network_path  Character. Path to the raw isoform interaction CSV.
#' @param pi_threshold  Numeric scalar or \code{NULL}. Set to \code{NULL} to
#'   use a pre-computed interaction column instead of a pi cutoff. Default
#'   \code{0.5}.
#' @param interact_col  Character scalar. Name of the binary (0/1) interaction
#'   column to use when \code{pi_threshold = NULL}. Ignored otherwise.
#'   E.g. \code{"conservative"} or \code{"soft"}.
#' @param kegg_genesets Optional pre-fetched MSigDB tibble. Pass a cached copy
#'   to avoid repeated downloads during iterative development.
#' @param plot          Logical. Save the full diagnostic plot bundle via
#'   \code{plot_cutoff_diagnostics()}. Default \code{FALSE}.
#' @param plot_dir      Character. Directory for diagnostic plots.
#'   Default \code{"data-raw/plots"}.
#' @param log_dir       Character. Directory for log files.
#'   Default \code{"data-raw/logs"}.
#' @param n_workers     Integer. Number of parallel workers. Defaults to all
#'   available cores minus one. If a multisession plan is already active with
#'   a different worker count, the plan is respawned with the requested count.
#'
#' @return Invisibly returns the result list from
#'   \code{compute_network_cutoff()}.
update_package_network <- function(network_path,
                                   pi_threshold     = 0.5,
                                   interact_col     = NULL,
                                   kegg_genesets    = NULL,
                                   min_neighbors    = 5L,
                                   min_pathway_pos  = 10L,
                                   min_coverage     = 0.5,
                                   plot             = FALSE,
                                   plot_dir         = "data-raw/plots",
                                   log_dir          = "data-raw/logs") {

  # ── Start logging ────────────────────────────────────────────────────────────
  log_file <- start_log(log_dir)
  on.exit(stop_log(), add = TRUE, after = FALSE)

  # ── Read input CSV ───────────────────────────────────────────────────────────
  cat(">>> Reading raw network CSV...\n")
  isoform_network_df <- read.csv(network_path, encoding = "utf-8")
  cat(sprintf("    Read %d rows from %s\n\n", nrow(isoform_network_df),
              network_path))

  # ── Build isoform-level igraph ───────────────────────────────────────────────
  cat(">>> Building isoform-level igraph...\n")
  .network_data <- .construct_igraph(isoform_network_df, pi_threshold,
                                     interact_col = interact_col)
  cat(sprintf("    igraph: %d vertices, %d edges\n\n",
              vcount(.network_data), ecount(.network_data)))

  # ── Run gene-level cutoff analysis ──────────────────────────────────────────
  result <- compute_network_cutoff(
    isoform_network_df = isoform_network_df,
    pi_threshold        = pi_threshold,
    interact_col        = interact_col,
    kegg_genesets       = kegg_genesets,
    min_neighbors       = min_neighbors,
    min_pathway_pos     = min_pathway_pos,
    min_coverage        = min_coverage,
    plot                = plot,
    plot_dir            = plot_dir,
    run_multipathway    = TRUE
  )

  .network_cutoff <- result$best_cutoff
  .pi_threshold   <- result$pi_threshold

  # ── Write internal package data ──────────────────────────────────────────────
  cat(">>> Writing internal package data to R/sysdata.rda...\n")
  usethis::use_data(
    .network_data,
    .network_cutoff,
    .pi_threshold,
    internal  = TRUE,
    overwrite = TRUE
  )
  cat(sprintf("    .network_data   : %d vertices, %d edges (isoform level)\n",
              vcount(.network_data), ecount(.network_data)))
  cat(sprintf("    .network_cutoff : %.4f (PR/F1-based q50)\n",
              .network_cutoff))
  if (is.null(pi_threshold)) {
    cat(sprintf("    .pi_threshold   : NULL (interact_col = '%s')\n", interact_col))
  } else {
    cat(sprintf("    .pi_threshold   : %.3f\n", .pi_threshold))
  }
  cat(sprintf("    Log file        : %s\n", log_file))
  cat("    Done. Rebuild the package to make the new data active.\n\n")

  invisible(result)
}


# ==============================================================================
# SECTION 2: Run the update
# ==============================================================================
result <- update_package_network(
  network_path  = "data-raw/results_PHYSICAL_700_FULL_threshold_22_06_26.csv",
  pi_threshold  = NULL,
  interact_col  = "soft",
  min_neighbors    = 5L,
  min_pathway_pos  = 10L,
  min_coverage     = 0.5,
  plot          = TRUE
)

# ==============================================================================
# SECTION 3: Post-hoc plots (can be re-run independently)
# ==============================================================================
plots_kegg <- plot_kegg_characterisation(result)
p_membership <- plot_gene_pathway_membership(result)

res <- build_kegg_category_table(result)
View(res$table)

n_other <- length(res$other_pathways)
n_total <- nrow(.build_pathway_summary(result))
cat(sprintf("%.1f%% of KEGG pathways (%d/%d) were unclassified ('Other').\n",
            100 * n_other / n_total, n_other, n_total))

# ==============================================================================
# SECTION 4: Cache result so plots can be re-styled without rerunning
# ==============================================================================
saveRDS(result, "data-raw/result_cache.rds")
