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
#     pi_threshold       = 0.5
#   )
# ==============================================================================


#' Compute the network-wide interaction cutoff from KEGG pathway benchmarking
#'
#' Accepts an already-loaded data frame — the caller is responsible for reading
#' the CSV and setting the parallel plan before calling this function. This
#' avoids redundant file I/O and double RAM usage when called from
#' \code{update_package_network()}.
#'
#' The analysis proceeds in ten steps:
#' \enumerate{
#'   \item Inspect the input.
#'   \item Aggregate to gene-level pairs.
#'   \item Compute gene-level degree.
#'   \item Load KEGG gene sets (or accept a pre-cached copy).
#'   \item Build neighbor lookup table.
#'   \item Per-pathway ROC + PR analysis (parallelised).
#'   \item AUC-weighted cutoff summaries.
#'   \item Pooled validation + direct sens/spec cutoff derivation.
#'   \item Multi-pathway gene handling strategy comparison.
#'   \item Choose and report the best cutoff.
#' }
#'
#' @param isoform_network_df Data frame loaded from the raw interaction CSV.
#'   Expected columns: \code{gene_1}, \code{gene_2}, \code{pi}.
#' @param pi_threshold Numeric scalar. Minimum \code{pi} value to call an
#'   isoform pair interacting. Default \code{0.5}.
#' @param kegg_genesets Optional pre-fetched MSigDB tibble (columns
#'   \code{gs_name}, \code{ensembl_gene}, \code{gs_url}). If \code{NULL}
#'   (default), fetched fresh from MSigDB. Pass a cached copy during iterative
#'   development to avoid repeated downloads.
#' @param plot Logical. Save diagnostic plots to \code{plot_dir}. Default
#'   \code{FALSE}.
#' @param plot_dir Character. Directory for plots. Created if absent.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{gene_pairs}{Gene-level edge table with pi summaries and
#'     \code{any_interact} flag.}
#'   \item{gene_degree}{Per-gene binary and pi-weighted degree.}
#'   \item{pathway_results}{Per-pathway ROC/PR metrics and best cutoffs.}
#'   \item{validation}{Precision/recall/FPR/F1/sens/spec for all candidate
#'     cutoffs.}
#'   \item{cutoff_roc_summary}{AUC-weighted quantile summary (ROC-based).}
#'   \item{cutoff_f1_summary}{AUC-weighted quantile summary (PR/F1-based).}
#'   \item{sens_spec_cutoffs}{Three-row tibble of direct sens/spec-based cutoffs.}
#'   \item{best_cutoff}{Recommended cutoff (pooled Youden's J).}
#'   \item{pi_threshold}{The \code{pi_threshold} value used.}
#'   \item{kegg_genesets}{MSigDB KEGG gene set table used.}
#'   \item{multipathway_comparison}{Named list with one entry per strategy
#'     (\code{"all"}, \code{"remove"}, \code{"random"}, \code{"best"}).}
#' }
compute_network_cutoff <- function(isoform_network_df,
                                   pi_threshold  = 0.5,
                                   kegg_genesets = NULL,
                                   plot          = FALSE,
                                   plot_dir      = "data-raw/plots") {

  cat("================================================================\n")
  cat("  KEGG Pathway Interaction Cutoff Analysis\n")
  cat("================================================================\n\n")

  # --------------------------------------------------------------------------
  # [1/10] Inspect input
  # --------------------------------------------------------------------------
  cat(">>> [1/10] Inspecting isoform interaction network...\n")
  cat(sprintf("    Rows (isoform pairs)    : %d\n",   nrow(isoform_network_df)))
  cat(sprintf("    Unique isoforms (col 1) : %d\n",   n_distinct(isoform_network_df$gene_1)))
  cat(sprintf("    Unique isoforms (col 2) : %d\n",   n_distinct(isoform_network_df$gene_2)))
  cat(sprintf("    pi threshold            : >= %.3f\n\n", pi_threshold))

  # --------------------------------------------------------------------------
  # [2/10] Aggregate to gene level
  # --------------------------------------------------------------------------
  cat(">>> [2/10] Aggregating to gene-level pairs...\n")
  gene_pairs    <- .build_gene_pairs(isoform_network_df, pi_threshold)
  n_interacting <- sum(gene_pairs$any_interact == 1)
  cat(sprintf("    Total gene pairs          : %d\n", nrow(gene_pairs)))
  cat(sprintf("    Pairs with any_interact=1 : %d\n", n_interacting))
  cat(sprintf("    Interaction rate          : %.1f%%\n\n",
              100 * n_interacting / nrow(gene_pairs)))

  # --------------------------------------------------------------------------
  # [3/10] Gene-level degree
  # --------------------------------------------------------------------------
  cat(">>> [3/10] Computing gene-level degree...\n")
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
  cat(sprintf("    Genes in network       : %d\n",     nrow(gene_degree)))
  cat(sprintf("    Median binary degree   : %.0f\n",   median(gene_degree$degree_binary)))
  cat(sprintf("    Max binary degree      : %.0f\n\n", max(gene_degree$degree_binary)))

  # --------------------------------------------------------------------------
  # [4/10] KEGG gene sets
  # --------------------------------------------------------------------------
  cat(">>> [4/10] Loading KEGG gene sets...\n")
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
  # [5/10] Neighbor lookup
  # --------------------------------------------------------------------------
  cat(">>> [5/10] Building neighbor lookup table...\n")
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
  # [6/10] Per-pathway ROC + PR analysis (baseline: all genes, no deduplication)
  # --------------------------------------------------------------------------
  cat(">>> [6/10] Running per-pathway ROC + PR analysis (baseline)...\n")
  cat(sprintf("    Scanning %d pathways\n\n", n_pathways))

  pathway_list <- kegg_genesets %>% group_by(gs_name) %>% group_split()

  pathway_results <- future_map_dfr(
    pathway_list,
    function(pathway_df) {
      .best_cutoff_for_pathway(
        pathway_name      = unique(pathway_df$gs_name),
        pathway_genes     = pathway_df$ensembl_gene,
        neighbor_list     = neighbor_list,
        all_network_genes = all_network_genes
      )
    },
    .options  = furrr_options(
      seed     = TRUE,
      packages = .furrr_pkgs,
      globals  = list(
        .best_cutoff_for_pathway   = .best_cutoff_for_pathway,
        .compute_scores_vectorized = .compute_scores_vectorized,
        neighbor_list              = neighbor_list,
        all_network_genes          = all_network_genes
      )
    ),
    .progress = TRUE
  )

  n_valid   <- sum(!is.na(pathway_results$cutoff_roc))
  n_skipped <- sum( is.na(pathway_results$cutoff_roc))
  cat(sprintf("\n    Valid pathways : %d | Skipped: %d\n", n_valid, n_skipped))
  cat(sprintf("    Median ROC-AUC : %.3f\n",
              median(pathway_results$auc_roc, na.rm = TRUE)))
  cat(sprintf("    Median PR-AUC  : %.3f\n",
              median(pathway_results$auc_pr,  na.rm = TRUE)))
  cat(sprintf("    Median sensitivity (best Youden's J) : %.3f\n",
              median(pathway_results$sensitivity, na.rm = TRUE)))
  cat(sprintf("    Median specificity (best Youden's J) : %.3f\n\n",
              median(pathway_results$specificity, na.rm = TRUE)))

  # --------------------------------------------------------------------------
  # [7/10] AUC-weighted cutoff summaries
  # --------------------------------------------------------------------------
  cat(">>> [7/10] Computing AUC-weighted conservative cutoffs...\n")

  cutoff_roc_summary <- pathway_results %>%
    filter(!is.na(cutoff_roc), auc_roc > 0.5) %>%
    mutate(w = auc_roc - 0.5) %>%
    summarise(
      n_pathways_used  = n(),
      q50              = .weighted_quantile(cutoff_roc, w, 0.50),
      q75              = .weighted_quantile(cutoff_roc, w, 0.75),
      q90              = .weighted_quantile(cutoff_roc, w, 0.90),
      mean_auc_roc     = weighted.mean(auc_roc, w),
      mean_sensitivity = weighted.mean(sensitivity, w, na.rm = TRUE),
      mean_specificity = weighted.mean(specificity, w, na.rm = TRUE)
    )

  cutoff_f1_summary <- pathway_results %>%
    filter(!is.na(cutoff_f1), !is.na(auc_pr)) %>%
    mutate(prevalence = n_pos / (n_pos + n_neg),
           w          = pmax(auc_pr - prevalence, 0)) %>%
    filter(w > 0) %>%
    summarise(
      n_pathways_used = n(),
      q50             = .weighted_quantile(cutoff_f1, w, 0.50),
      q75             = .weighted_quantile(cutoff_f1, w, 0.75),
      q90             = .weighted_quantile(cutoff_f1, w, 0.90),
      mean_auc_pr     = weighted.mean(auc_pr, w)
    )

  cat(sprintf("    Weighted mean AUC — ROC: %.3f | PR: %.3f\n",
              cutoff_roc_summary$mean_auc_roc, cutoff_f1_summary$mean_auc_pr))
  cat(sprintf("    Weighted mean sensitivity : %.3f | specificity : %.3f\n",
              cutoff_roc_summary$mean_sensitivity,
              cutoff_roc_summary$mean_specificity))
  cat(sprintf("    Cutoff q50 — ROC: %.4f | PR/F1: %.4f\n",
              cutoff_roc_summary$q50, cutoff_f1_summary$q50))
  cat(sprintf("    Cutoff q75 — ROC: %.4f | PR/F1: %.4f  <-- recommended\n",
              cutoff_roc_summary$q75, cutoff_f1_summary$q75))
  cat(sprintf("    Cutoff q90 — ROC: %.4f | PR/F1: %.4f\n\n",
              cutoff_roc_summary$q90, cutoff_f1_summary$q90))

  # --------------------------------------------------------------------------
  # [8/10] Pooled validation + direct sens/spec cutoff derivation
  #
  # all_scores is large (~n_genes × n_pathways rows). It is kept alive through
  # step 9 for the strategy comparison, then explicitly discarded.
  # --------------------------------------------------------------------------
  cat(">>> [8/10] Pooled validation across all gene x pathway observations...\n")

  all_scores <- future_map_dfr(
    pathway_list,
    function(pathway_df) {
      .compute_scores_vectorized(
        pathway_df$ensembl_gene, neighbor_list, all_network_genes
      ) %>%
        mutate(pathway = unique(pathway_df$gs_name))
    },
    .options  = furrr_options(
      seed     = TRUE,
      packages = .furrr_pkgs,
      globals  = list(
        .compute_scores_vectorized = .compute_scores_vectorized,
        neighbor_list              = neighbor_list,
        all_network_genes          = all_network_genes
      )
    ),
    .progress = TRUE
  )

  # 8a. Validate the six AUC-weighted quantile cutoffs
  candidate_cutoffs <- c(
    ROC_q50   = cutoff_roc_summary$q50,
    ROC_q75   = cutoff_roc_summary$q75,
    ROC_q90   = cutoff_roc_summary$q90,
    PR_F1_q50 = cutoff_f1_summary$q50,
    PR_F1_q75 = cutoff_f1_summary$q75,
    PR_F1_q90 = cutoff_f1_summary$q90
  )

  validation <- imap_dfr(candidate_cutoffs, function(co, nm) {
    parts  <- strsplit(nm, "_")[[1]]
    method <- paste(parts[-length(parts)], collapse = "/")
    quant  <- parts[length(parts)]
    .validate_cutoff(co, all_scores) %>%
      mutate(method = method, quantile = quant,
             label  = paste0(method, " | ", quant))
  }) %>%
    select(label, method, quantile, everything())

  # 8b. Derive cutoffs directly from pooled sens/spec
  cat("         Deriving direct sens/spec cutoffs from pooled scores...\n")
  sens_spec_cutoffs <- .derive_sens_spec_cutoffs(all_scores)

  # Append to the main validation table so all cutoffs live in one place
  validation <- bind_rows(
    validation,
    imap_dfr(seq_len(nrow(sens_spec_cutoffs)), function(i, ...) {
      .validate_cutoff(sens_spec_cutoffs$cutoff[i], all_scores) %>%
        mutate(
          method   = "SensSpec",
          quantile = recode(
            sens_spec_cutoffs$method[i],
            "Pooled Youden's J"         = "youden",
            "Sensitivity = Specificity" = "equal",
            "Min distance to (0,1)"     = "mindist"
          ),
          label = paste0("SensSpec | ", sens_spec_cutoffs$method[i])
        )
    })
  ) %>%
    select(label, method, quantile, everything())

  cat(sprintf(
    "         Pooled Youden's J cutoff          : %.4f  (sens=%.3f, spec=%.3f)\n",
    sens_spec_cutoffs$cutoff[sens_spec_cutoffs$method == "Pooled Youden's J"],
    sens_spec_cutoffs$sensitivity[sens_spec_cutoffs$method == "Pooled Youden's J"],
    sens_spec_cutoffs$specificity[sens_spec_cutoffs$method == "Pooled Youden's J"]
  ))
  cat(sprintf(
    "         Sensitivity = Specificity cutoff  : %.4f  (sens=%.3f, spec=%.3f)\n",
    sens_spec_cutoffs$cutoff[sens_spec_cutoffs$method == "Sensitivity = Specificity"],
    sens_spec_cutoffs$sensitivity[sens_spec_cutoffs$method == "Sensitivity = Specificity"],
    sens_spec_cutoffs$specificity[sens_spec_cutoffs$method == "Sensitivity = Specificity"]
  ))
  cat(sprintf(
    "         Min distance to (0,1) cutoff      : %.4f  (sens=%.3f, spec=%.3f)\n\n",
    sens_spec_cutoffs$cutoff[sens_spec_cutoffs$method == "Min distance to (0,1)"],
    sens_spec_cutoffs$sensitivity[sens_spec_cutoffs$method == "Min distance to (0,1)"],
    sens_spec_cutoffs$specificity[sens_spec_cutoffs$method == "Min distance to (0,1)"]
  ))

  # --------------------------------------------------------------------------
  # [9/10] Multi-pathway gene handling strategies
  # --------------------------------------------------------------------------
  cat(">>> [9/10] Comparing multi-pathway gene handling strategies...\n")

  strategies <- c("remove", "random", "best")
  multipathway_comparison <- map(strategies, function(strat) {
    .run_strategy(kegg_genesets, neighbor_list, all_network_genes,
                  strategy       = strat,
                  all_scores_ref = all_scores)
  })
  names(multipathway_comparison) <- strategies

  multipathway_comparison[["all"]] <- list(
    strategy           = "all",
    cutoff_roc_summary = cutoff_roc_summary,
    cutoff_f1_summary  = cutoff_f1_summary,
    validation         = validation %>%
      filter(quantile == "q75") %>%
      mutate(strategy = "all"),
    n_genes_used       = n_distinct(kegg_genesets$ensembl_gene),
    n_pathways_used    = n_valid,
    pathway_results    = pathway_results
  )

  # Discard all_scores — too large to carry in the result list
  rm(all_scores)

  # --------------------------------------------------------------------------
  # [10/10] Choose best cutoff and report
  # --------------------------------------------------------------------------
  best_cutoff <- cutoff_roc_summary$q75

  cat("\n================================================================\n")
  cat("  ANALYSIS COMPLETE — Summary\n")
  cat("================================================================\n")
  cat(sprintf("  Gene pairs in network              : %d\n",  nrow(gene_pairs)))
  cat(sprintf("  Interacting pairs (any_interact=1) : %d\n",  n_interacting))
  cat(sprintf("  KEGG pathways scanned              : %d\n",  n_pathways))
  cat(sprintf("  Pathways with valid cutoff         : %d\n",  n_valid))
  cat(sprintf("  Median ROC-AUC                     : %.3f\n",
              median(pathway_results$auc_roc, na.rm = TRUE)))
  cat(sprintf("  Median PR-AUC                      : %.3f\n",
              median(pathway_results$auc_pr,  na.rm = TRUE)))
  cat("  -----------------------------------------------\n")
  cat("  Indirect cutoffs (AUC-weighted quantile of per-pathway bests):\n")
  roc_q75_v <- validation %>% filter(method == "ROC",   quantile == "q75")
  f1_q75_v  <- validation %>% filter(method == "PR/F1", quantile == "q75")
  cat(sprintf(
    "    ROC-based q75   : %.4f  (sens=%.3f, spec=%.3f, prec=%.3f, F1=%.3f)\n",
    cutoff_roc_summary$q75,
    roc_q75_v$sensitivity, roc_q75_v$specificity,
    roc_q75_v$precision,   roc_q75_v$f1))
  cat(sprintf(
    "    PR/F1-based q75 : %.4f  (sens=%.3f, spec=%.3f, prec=%.3f, F1=%.3f)\n",
    cutoff_f1_summary$q75,
    f1_q75_v$sensitivity, f1_q75_v$specificity,
    f1_q75_v$precision,   f1_q75_v$f1))
  cat("  -----------------------------------------------\n")
  cat("  Direct cutoffs (derived from pooled sensitivity/specificity):\n")
  for (i in seq_len(nrow(sens_spec_cutoffs))) {
    r <- sens_spec_cutoffs[i, ]
    cat(sprintf(
      "    %-34s : %.4f  (sens=%.3f, spec=%.3f, prec=%.3f, F1=%.3f)\n",
      r$method, r$cutoff, r$sensitivity, r$specificity,
      r$precision, r$f1))
  }
  cat(sprintf("\n  >>> Returning best cutoff (ROC-based q75): %.4f\n",
              best_cutoff))
  cat("================================================================\n\n")

  result <- list(
    gene_pairs              = gene_pairs,
    gene_degree             = gene_degree,
    pathway_results         = pathway_results,
    validation              = validation,
    cutoff_roc_summary      = cutoff_roc_summary,
    cutoff_f1_summary       = cutoff_f1_summary,
    sens_spec_cutoffs       = sens_spec_cutoffs,
    best_cutoff             = best_cutoff,
    pi_threshold            = pi_threshold,
    kegg_genesets           = kegg_genesets,
    multipathway_comparison = multipathway_comparison
  )

  if (plot) plot_cutoff_diagnostics(result, plot_dir = plot_dir)

  result
}
