# ==============================================================================
# data-raw/network_analysis.R
#
# Main analysis function: compute_network_cutoff()
#
# Depends on: network_helpers.R  (must be sourced first)
#
# USAGE:
#   source("data-raw/network_helpers.R")
#   source("data-raw/network_analysis.R")
#
#   result <- compute_network_cutoff(
#     isoform_network_df = read.csv("data-raw/results_*.csv"),
#     pi_threshold       = NULL,
#     interact_col       = "soft",
#     min_neighbors      = 5L,    # exclude genes with < 5 interaction partners
#     min_pathway_pos    = 10L    # exclude pathways with < 10 positives after
#                                 # neighbour filtering
#   )
# ==============================================================================


#' Compute the network-wide interaction cutoff from KEGG pathway benchmarking
#'
#' Accepts an already-loaded data frame — the caller is responsible for reading
#' the CSV before calling this function.
#'
#' The analysis proceeds in nine steps:
#' \enumerate{
#'   \item Inspect the input and echo filter settings.
#'   \item Aggregate to gene-level pairs.
#'   \item Compute gene-level degree.
#'   \item Load KEGG gene sets (or accept a pre-cached copy).
#'   \item Build neighbor lookup table.
#'   \item Per-pathway ROC + PR analysis.
#'   \item AUC-weighted cutoff summaries.
#'   \item Pooled validation of the recommended cutoff candidates.
#'   \item Multi-pathway gene handling strategy comparison (optional).
#'   \item Choose and report the best cutoff.
#' }
#'
#' @param isoform_network_df Data frame loaded from the raw interaction CSV.
#'   Expected columns: \code{gene_1}, \code{gene_2}, \code{pi}, and optionally
#'   the column named by \code{interact_col}.
#' @param pi_threshold Numeric scalar or \code{NULL}. Minimum \code{pi} value
#'   to call an isoform pair interacting. Set to \code{NULL} to use a
#'   pre-computed interaction column instead (see \code{interact_col}).
#'   Default \code{0.5}.
#' @param interact_col Character scalar. Name of a binary (0/1) column in
#'   \code{isoform_network_df} that directly encodes whether each isoform pair
#'   interacts. Only used when \code{pi_threshold = NULL}. Must be present in
#'   the data frame.
#' @param kegg_genesets Optional pre-fetched MSigDB tibble (columns
#'   \code{gs_name}, \code{ensembl_gene}, \code{gs_url}). If \code{NULL}
#'   (default), fetched fresh from MSigDB. Pass a cached copy during iterative
#'   development to avoid repeated downloads.
#' @param min_neighbors Integer. Minimum number of interaction partners a gene
#'   must have for its ratio score to be included in any pathway analysis.
#'   Genes below this floor have high-variance scores (a single partner
#'   in/out of the pathway dominates the ratio) and are excluded from scoring.
#'   They remain as potential positives for computing \code{n_pos} counts.
#'   Default \code{1L} (no filter). Recommended: \code{5L}.
#' @param min_pathway_pos Integer. Minimum number of pathway-positive genes
#'   that must remain in the scored set (after neighbour filtering) for a
#'   pathway to be included in the analysis. Pathways below this floor are
#'   skipped (all-NA row in pathway_results). Default \code{1L} (original
#'   behaviour). Recommended: \code{10L}.
#' @param plot Logical. Save diagnostic plots to \code{plot_dir}. Default
#'   \code{FALSE}.
#' @param plot_dir Character. Directory for plots. Created if absent.
#' @param run_multipathway Logical. Run the multi-pathway gene handling
#'   strategy comparison. Default \code{TRUE}.
#' @param min_coverage Numeric in [0, 1]. Minimum fraction of a pathway's
#'   genes that must be present in the network (i.e. have >= 1 interaction
#'   partner) for the pathway to be included in the analysis. Pathways below
#'   this floor are dropped entirely before scoring, since their cutoff
#'   estimate would be based on a small, unrepresentative subset of the
#'   pathway. Default \code{0} (no filter). Recommended: \code{0.3}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{gene_pairs}{Gene-level edge table with pi summaries and
#'     \code{any_interact} flag.}
#'   \item{gene_degree}{Per-gene binary and pi-weighted degree.}
#'   \item{pathway_results}{Per-pathway ROC/PR metrics and best cutoffs.}
#'   \item{validation}{Precision/recall/FPR/F1/sens/spec for the PR/F1
#'     candidate cutoffs, validated on the pooled gene × pathway scores.}
#'   \item{cutoff_roc_summary}{AUC-weighted quantile summary (ROC-based).}
#'   \item{cutoff_f1_summary}{AUC-weighted quantile summary (PR/F1-based).}
#'   \item{best_cutoff}{Recommended cutoff (PR/F1-based unweighted q50).}
#'   \item{pi_threshold}{The \code{pi_threshold} value used.}
#'   \item{interact_col}{The \code{interact_col} value used.}
#'   \item{min_neighbors}{The \code{min_neighbors} value used.}
#'   \item{min_pathway_pos}{The \code{min_pathway_pos} value used.}
#'   \item{kegg_genesets}{MSigDB KEGG gene set table used.}
#'   \item{multipathway_comparison}{Named list with one entry per strategy
#'     (\code{"all"}, \code{"remove"}, \code{"random"}, \code{"best"}), or
#'     \code{NULL} when \code{run_multipathway = FALSE}.}
#' }
compute_network_cutoff <- function(isoform_network_df,
                                   pi_threshold     = 0.5,
                                   interact_col     = NULL,
                                   kegg_genesets    = NULL,
                                   min_neighbors    = 1L,
                                   min_pathway_pos  = 1L,
                                   min_coverage     = 0,
                                   plot             = TRUE,
                                   plot_dir         = "data-raw/plots",
                                   run_multipathway = TRUE) {

  min_neighbors   <- as.integer(min_neighbors)
  min_pathway_pos <- as.integer(min_pathway_pos)
  if (!is.numeric(min_coverage) || min_coverage < 0 || min_coverage > 1)
    stop("`min_coverage` must be a single numeric value in [0, 1].", call. = FALSE)

  cat("================================================================\n")
  cat("  KEGG Pathway Interaction Cutoff Analysis\n")
  cat("================================================================\n\n")

  # --------------------------------------------------------------------------
  # [1/9] Inspect input
  # --------------------------------------------------------------------------
  cat(">>> [1/9] Inspecting isoform interaction network...\n")
  cat(sprintf("    Rows (isoform pairs)    : %d\n",   nrow(isoform_network_df)))
  cat(sprintf("    Unique genes (col 1)    : %d\n",   n_distinct(isoform_network_df$gene_1)))
  cat(sprintf("    Unique genes (col 2)    : %d\n",   n_distinct(isoform_network_df$gene_2)))
  if (is.null(pi_threshold)) {
    if (is.null(interact_col))
      stop("compute_network_cutoff: 'interact_col' must be supplied when pi_threshold = NULL")
    cat(sprintf("    Interaction mode        : column '%s'\n", interact_col))
  } else {
    cat(sprintf("    pi threshold            : >= %.3f\n", pi_threshold))
  }
  cat(sprintf("    min_neighbors filter    : >= %d interaction partners\n",
              min_neighbors))
  cat(sprintf("    min_pathway_pos filter  : >= %d positives after neighbour filter\n",
              min_pathway_pos))
  cat(sprintf("    min_coverage filter     : >= %.2f fraction of pathway genes in network\n\n",
              min_coverage))

  # --------------------------------------------------------------------------
  # [2/9] Aggregate to gene level
  # --------------------------------------------------------------------------
  cat(">>> [2/9] Aggregating to gene-level pairs...\n")
  gene_pairs    <- .build_gene_pairs(isoform_network_df, pi_threshold,
                                     interact_col = interact_col)
  n_interacting <- sum(gene_pairs$any_interact == 1)
  cat(sprintf("    Total gene pairs          : %d\n", nrow(gene_pairs)))
  cat(sprintf("    Pairs with any_interact=1 : %d\n", n_interacting))
  cat(sprintf("    Interaction rate          : %.1f%%\n\n",
              100 * n_interacting / nrow(gene_pairs)))

  # --------------------------------------------------------------------------
  # [3/9] Gene-level degree
  # --------------------------------------------------------------------------
  cat(">>> [3/9] Computing gene-level degree...\n")
  degree_A <- gene_pairs %>%
    group_by(gene = gene_A) %>%
    summarise(degree_binary = sum(any_interact),
              degree_pi     = sum(mean_pi), .groups = "drop")
  degree_B <- gene_pairs %>%
    group_by(gene = gene_B) %>%
    summarise(degree_binary = sum(any_interact),
              degree_pi     = sum(mean_pi), .groups = "drop")
  gene_degree <- bind_rows(degree_A, degree_B) %>%
    group_by(gene) %>%
    summarise(degree_binary = sum(degree_binary),
              degree_pi     = sum(degree_pi), .groups = "drop")
  n_low_degree <- sum(gene_degree$degree_binary < min_neighbors)
  cat(sprintf("    Genes in network       : %d\n",     nrow(gene_degree)))
  cat(sprintf("    Median binary degree   : %.0f\n",   median(gene_degree$degree_binary)))
  cat(sprintf("    Max binary degree      : %.0f\n",   max(gene_degree$degree_binary)))
  cat(sprintf("    Genes below min_neighbors (%d) : %d (%.1f%% excluded from scoring)\n\n",
              min_neighbors, n_low_degree,
              100 * n_low_degree / nrow(gene_degree)))

  # --------------------------------------------------------------------------
  # [4/9] KEGG gene sets
  # --------------------------------------------------------------------------
  cat(">>> [4/9] Loading KEGG gene sets...\n")
  if (is.null(kegg_genesets)) {
    kegg_genesets <- bind_rows(
      msigdbr(species = "Homo sapiens", collection = "C2",
              subcollection = "CP:KEGG_LEGACY") %>%
        select(gs_name, ensembl_gene, gs_url),
      msigdbr(species = "Homo sapiens", collection = "C2",
              subcollection = "CP:KEGG_MEDICUS") %>%
        select(gs_name, ensembl_gene, gs_url)
    )
    cat("    Fetched fresh from MSigDB\n")
  } else {
    cat("    Using pre-supplied gene set table\n")
  }
  n_pathways   <- n_distinct(kegg_genesets$gs_name)
  n_kegg_genes <- n_distinct(kegg_genesets$ensembl_gene)
  n_multi      <- kegg_genesets %>%
    group_by(ensembl_gene) %>%
    summarise(n = n_distinct(gs_name), .groups = "drop") %>%
    filter(n > 1) %>% nrow()
  cat(sprintf("    Total unique pathways  : %d\n",   n_pathways))
  cat(sprintf("    Total unique genes     : %d\n",   n_kegg_genes))
  cat(sprintf("    Multi-pathway genes    : %d (%.1f%%)\n\n",
              n_multi, 100 * n_multi / n_kegg_genes))

  # --------------------------------------------------------------------------
  # [5/9] Neighbor lookup
  # --------------------------------------------------------------------------
  cat(">>> [5/9] Building neighbor lookup table...\n")
  neighbor_list     <- .build_neighbor_list(gene_pairs)
  all_network_genes <- names(neighbor_list)
  n_overlap         <- length(intersect(kegg_genesets$ensembl_gene,
                                        all_network_genes))
  cat(sprintf("    Genes with >= 1 interaction : %d\n",
              length(all_network_genes)))
  cat(sprintf("    KEGG genes in network       : %d / %d (%.1f%%)\n\n",
              n_overlap, n_kegg_genes,
              100 * n_overlap / n_kegg_genes))

  # --------------------------------------------------------------------------
  # [5b] Coverage pre-filter
  # --------------------------------------------------------------------------
  cat(">>> Applying min_coverage filter...\n")

  pathway_coverage <- kegg_genesets %>%
    group_by(gs_name) %>%
    summarise(
      n_pathway_genes = n_distinct(ensembl_gene),
      n_in_network    = length(intersect(unique(ensembl_gene), all_network_genes)),
      coverage        = n_in_network / n_pathway_genes,
      .groups = "drop"
    )

  low_cov_pathways <- pathway_coverage$gs_name[pathway_coverage$coverage < min_coverage]
  n_pathways_total <- nrow(pathway_coverage)
  n_low_coverage   <- length(low_cov_pathways)

  cat(sprintf("    Pathways below min_coverage (%.2f) : %d / %d (%.1f%%) excluded\n",
              min_coverage, n_low_coverage, n_pathways_total,
              100 * n_low_coverage / n_pathways_total))
  cat(sprintf("    Median coverage (all pathways)     : %.2f\n",
              median(pathway_coverage$coverage)))
  cat(sprintf("    Median coverage (retained)         : %.2f\n\n",
              median(pathway_coverage$coverage[pathway_coverage$coverage >= min_coverage])))

  kegg_genesets <- kegg_genesets %>% filter(!gs_name %in% low_cov_pathways)
  n_pathways    <- n_distinct(kegg_genesets$gs_name)

  # --------------------------------------------------------------------------
  # [6/9] Per-pathway ROC + PR analysis
  # --------------------------------------------------------------------------
  cat(">>> [6/9] Running per-pathway ROC + PR analysis...\n")
  cat(sprintf("    Scanning %d pathways  |  min_neighbors = %d  |  min_pathway_pos = %d\n\n",
              n_pathways, min_neighbors, min_pathway_pos))

  pathway_list <- kegg_genesets %>% group_by(gs_name) %>% group_split()

  pathway_results <- dplyr::bind_rows(lapply(
    pathway_list,
    FUN = function(pathway_df) {
      .best_cutoff_for_pathway(
        pathway_name      = unique(pathway_df$gs_name),
        pathway_genes     = pathway_df$ensembl_gene,
        neighbor_list     = neighbor_list,
        all_network_genes = all_network_genes,
        min_neighbors     = min_neighbors,
        min_pathway_pos   = min_pathway_pos
      )
    }
  ))

  n_valid   <- sum(!is.na(pathway_results$cutoff_f1))
  n_skipped <- sum( is.na(pathway_results$cutoff_f1))
  cat(sprintf("    Valid pathways : %d | Skipped (below filters): %d\n", n_valid, n_skipped))
  cat(sprintf("    Median ROC-AUC : %.3f\n",
              median(pathway_results$auc_roc, na.rm = TRUE)))
  cat(sprintf("    Median PR-AUC  : %.3f\n",
              median(pathway_results$auc_pr,  na.rm = TRUE)))
  cat(sprintf("    Median sensitivity (best Youden's J) : %.3f\n",
              median(pathway_results$sensitivity, na.rm = TRUE)))
  cat(sprintf("    Median specificity (best Youden's J) : %.3f\n\n",
              median(pathway_results$specificity, na.rm = TRUE)))

  # --------------------------------------------------------------------------
  # [7/9] Cutoff summaries — AUC-weighted and unweighted (raw count)
  # --------------------------------------------------------------------------
  cat(">>> [7/9] Computing cutoff summaries (AUC-weighted and unweighted)...\n")

  # ── AUC-weighted summaries ─────────────────────────────────────────────────
  cutoff_roc_summary <- pathway_results %>%
    filter(!is.na(cutoff_roc), auc_roc > 0.5) %>%
    mutate(w = auc_roc - 0.5) %>%
    summarise(
      n_pathways_used       = n(),
      # AUC-weighted
      q50                   = .weighted_quantile(cutoff_roc, w, 0.50),
      q75                   = .weighted_quantile(cutoff_roc, w, 0.75),
      q90                   = .weighted_quantile(cutoff_roc, w, 0.90),
      mean_auc_roc          = weighted.mean(auc_roc, w),
      mean_sensitivity      = weighted.mean(sensitivity, w, na.rm = TRUE),
      mean_specificity      = weighted.mean(specificity, w, na.rm = TRUE),
      # Unweighted (raw count)
      unweighted_q50        = quantile(cutoff_roc, 0.50),
      unweighted_q75        = quantile(cutoff_roc, 0.75),
      unweighted_q90        = quantile(cutoff_roc, 0.90),
      unweighted_mean       = mean(cutoff_roc),
      unweighted_median     = median(cutoff_roc),
      unweighted_mean_auc   = mean(auc_roc),
      unweighted_mean_sens  = mean(sensitivity, na.rm = TRUE),
      unweighted_mean_spec  = mean(specificity, na.rm = TRUE)
    )

  cutoff_f1_summary <- pathway_results %>%
    filter(!is.na(cutoff_f1), !is.na(auc_pr)) %>%
    mutate(prevalence = n_pos / (n_pos + n_neg),
           w          = pmax(auc_pr - prevalence, 0)) %>%
    filter(w > 0) %>%
    summarise(
      n_pathways_used     = n(),
      # AUC-weighted
      q50                 = .weighted_quantile(cutoff_f1, w, 0.50),
      q75                 = .weighted_quantile(cutoff_f1, w, 0.75),
      q90                 = .weighted_quantile(cutoff_f1, w, 0.90),
      mean_auc_pr         = weighted.mean(auc_pr, w),
      # Unweighted (raw count)
      unweighted_q50      = quantile(cutoff_f1, 0.50),
      unweighted_q75      = quantile(cutoff_f1, 0.75),
      unweighted_q90      = quantile(cutoff_f1, 0.90),
      unweighted_mean     = mean(cutoff_f1),
      unweighted_median   = median(cutoff_f1),
      unweighted_mean_auc = mean(auc_pr)
    )

  cat("  --- AUC-weighted ---\n")
  cat(sprintf("    Weighted mean AUC — ROC: %.3f | PR: %.3f\n",
              cutoff_roc_summary$mean_auc_roc, cutoff_f1_summary$mean_auc_pr))
  cat(sprintf("    Weighted mean sensitivity : %.3f | specificity : %.3f\n",
              cutoff_roc_summary$mean_sensitivity,
              cutoff_roc_summary$mean_specificity))
  cat(sprintf("    Cutoff q50 — ROC: %.4f | PR/F1: %.4f\n",
              cutoff_roc_summary$q50, cutoff_f1_summary$q50))
  cat(sprintf("    Cutoff q75 — ROC: %.4f | PR/F1: %.4f\n",
              cutoff_roc_summary$q75, cutoff_f1_summary$q75))
  cat(sprintf("    Cutoff q90 — ROC: %.4f | PR/F1: %.4f\n",
              cutoff_roc_summary$q90, cutoff_f1_summary$q90))
  cat("  --- Unweighted (raw count) ---\n")
  cat(sprintf("    Mean AUC — ROC: %.3f | PR: %.3f\n",
              cutoff_roc_summary$unweighted_mean_auc,
              cutoff_f1_summary$unweighted_mean_auc))
  cat(sprintf("    Mean sensitivity : %.3f | specificity : %.3f\n",
              cutoff_roc_summary$unweighted_mean_sens,
              cutoff_roc_summary$unweighted_mean_spec))
  cat(sprintf("    Cutoff median — ROC: %.4f | PR/F1: %.4f\n",
              cutoff_roc_summary$unweighted_median,
              cutoff_f1_summary$unweighted_median))
  cat(sprintf("    Cutoff mean   — ROC: %.4f | PR/F1: %.4f\n",
              cutoff_roc_summary$unweighted_mean,
              cutoff_f1_summary$unweighted_mean))
  cat(sprintf("    Cutoff q50    — ROC: %.4f | PR/F1: %.4f    <-- recommended\n",
              cutoff_roc_summary$unweighted_q50, cutoff_f1_summary$unweighted_q50))
  cat(sprintf("    Cutoff q75    — ROC: %.4f | PR/F1: %.4f\n",
              cutoff_roc_summary$unweighted_q75, cutoff_f1_summary$unweighted_q75))
  cat(sprintf("    Cutoff q90    — ROC: %.4f | PR/F1: %.4f\n\n",
              cutoff_roc_summary$unweighted_q90, cutoff_f1_summary$unweighted_q90))

  # --------------------------------------------------------------------------
  # [8/9] Pooled validation of PR/F1 candidate cutoffs
  #
  # all_scores is slimmed to (gene, ratio, in_pathway, pathway) before
  # bind_rows to avoid carrying the heavy neighbors list-column in memory.
  # --------------------------------------------------------------------------
  cat(">>> [8/9] Pooled validation of PR/F1 candidate cutoffs...\n")

  all_scores <- dplyr::bind_rows(lapply(
    pathway_list,
    FUN = function(pathway_df) {
      .compute_scores_vectorized(
        pathway_df$ensembl_gene, neighbor_list, all_network_genes,
        min_neighbors = min_neighbors
      ) %>%
        dplyr::select(gene, ratio, in_pathway) %>%
        dplyr::mutate(pathway = unique(pathway_df$gs_name))
    }
  ))

  # Validate the PR/F1 quantile cutoffs only (ROC cutoffs are noted but not
  # used as the recommended value; pooled-Youden / SensSpec methods are
  # omitted — they degenerate under heavy class imbalance).
  candidate_cutoffs <- c(
    PR_F1_q50 = unname(cutoff_f1_summary$unweighted_q50),
    PR_F1_q75 = unname(cutoff_f1_summary$unweighted_q75),
    PR_F1_q90 = unname(cutoff_f1_summary$unweighted_q90),
    ROC_q50 = unname(cutoff_roc_summary$unweighted_q50),
    ROC_q75 = unname(cutoff_roc_summary$unweighted_q75),
    ROC_q90 = unname(cutoff_roc_summary$unweighted_q90)
  )

  validation <- imap_dfr(candidate_cutoffs, function(co, nm) {
    parts  <- strsplit(nm, "_")[[1]]
    method <- paste(parts[-length(parts)], collapse = "/")
    quant  <- paste0("unweighted_", parts[length(parts)])
    .validate_cutoff(co, all_scores) %>%
      mutate(method = method, quantile = quant,
             label  = paste0(method, " | ", quant))
  }) %>%
    select(label, method, quantile, everything())

  for (q in c("unweighted_q50", "unweighted_q75", "unweighted_q90")) {

    roc_row <- validation %>% filter(method == "ROC", quantile == q)
    f1_row  <- validation %>% filter(method == "PR/F1", quantile == q)

    cat(sprintf(
      "%s | ROC:   cutoff=%.4f  sens=%.3f  spec=%.3f  F1=%.3f\n",
      q, roc_row$cutoff, roc_row$sensitivity, roc_row$specificity, roc_row$f1))
    cat(sprintf(
      "%s | PR/F1: cutoff=%.4f  prec=%.3f  recall=%.3f  F1=%.3f\n\n",
      q, f1_row$cutoff, f1_row$precision, f1_row$recall, f1_row$f1))
  }

  # --------------------------------------------------------------------------
  # [9/9] Multi-pathway gene handling strategies (optional)
  # --------------------------------------------------------------------------
  if (run_multipathway) {
    cat(">>> [9/9] Comparing multi-pathway gene handling strategies...\n")

    strategies <- c("remove", "random", "best")
    multipathway_comparison <- lapply(strategies, function(strat) {
      cat(sprintf("    Strategy: %s\n", strat))
      .run_strategy(
        kegg_genesets, neighbor_list, all_network_genes,
        strategy        = strat,
        all_scores_ref  = all_scores,
        min_neighbors   = min_neighbors,
        min_pathway_pos = min_pathway_pos
      )
    })
    names(multipathway_comparison) <- strategies

    multipathway_comparison[["all"]] <- list(
      strategy           = "all",
      cutoff_roc_summary = cutoff_roc_summary,
      cutoff_f1_summary  = cutoff_f1_summary,
      validation         = validation %>%
        filter(quantile == "unweighted_q50", method == "PR/F1") %>%
        mutate(strategy = "all"),
      n_genes_used       = n_distinct(kegg_genesets$ensembl_gene),
      n_pathways_used    = n_valid,
      pathway_results    = pathway_results
    )
  } else {
    cat(">>> [9/9] Skipping multi-pathway gene handling strategies",
        "(run_multipathway = FALSE)...\n")
    multipathway_comparison <- NULL
  }

  # Discard all_scores — not needed in the returned result list
  rm(all_scores)

  # --------------------------------------------------------------------------
  # Final: choose best cutoff and report
  # --------------------------------------------------------------------------
  best_cutoff <- cutoff_f1_summary$unweighted_q50

  cat("\n================================================================\n")
  cat("  ANALYSIS COMPLETE — Summary\n")
  cat("================================================================\n")
  cat(sprintf("  Gene pairs in network              : %d\n",  nrow(gene_pairs)))
  cat(sprintf("  Interacting pairs (any_interact=1) : %d\n",  n_interacting))
  cat(sprintf("  Genes excluded (min_neighbors < %d) : %d\n",
              min_neighbors, n_low_degree))
  cat(sprintf("  KEGG pathways scanned (total)      : %d\n",  n_pathways_total))
  cat(sprintf("  Pathways excluded (min_coverage)   : %d\n",  n_low_coverage))
  cat(sprintf("  Pathways remaining after coverage  : %d\n",  n_pathways))
  cat(sprintf("  Pathways used (>= %d pos, valid)   : %d\n",
              min_pathway_pos, n_valid))
  cat(sprintf("  Pathways skipped (neighbour/pos)   : %d\n",  n_skipped))
  cat(sprintf("  Median ROC-AUC                     : %.3f\n",
              median(pathway_results$auc_roc, na.rm = TRUE)))
  cat(sprintf("  Median PR-AUC                      : %.3f\n",
              median(pathway_results$auc_pr,  na.rm = TRUE)))
  cat("  -----------------------------------------------\n")
  cat("  Cutoffs — unweighted quantile of per-pathway PR/F1 bests:\n")
  cat(sprintf(
    "    PR/F1-based median : %.4f | q50 : %.4f | q75 : %.4f | q90 : %.4f\n",
    cutoff_f1_summary$unweighted_median,
    cutoff_f1_summary$unweighted_q50,
    cutoff_f1_summary$unweighted_q75,
    cutoff_f1_summary$unweighted_q90))
  cat(sprintf("\n  >>> Returning best cutoff (PR/F1-based unweighted q50): %.4f\n",
              best_cutoff))
  cat("================================================================\n\n")

  result <- list(
    gene_pairs              = gene_pairs,
    gene_degree              = gene_degree,
    pathway_coverage         = pathway_coverage,
    pathway_results          = pathway_results,
    validation               = validation,
    cutoff_roc_summary       = cutoff_roc_summary,
    cutoff_f1_summary        = cutoff_f1_summary,
    best_cutoff              = best_cutoff,
    pi_threshold             = pi_threshold,
    interact_col             = interact_col,
    min_neighbors            = min_neighbors,
    min_pathway_pos          = min_pathway_pos,
    min_coverage             = min_coverage,
    kegg_genesets            = kegg_genesets,
    multipathway_comparison  = multipathway_comparison
  )

  if (plot) plot_cutoff_diagnostics(result, plot_dir = plot_dir)

  result
}
