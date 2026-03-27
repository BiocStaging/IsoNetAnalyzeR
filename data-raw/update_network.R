# ==============================================================================
# data-raw/update_network.R
#
# Recompute the interaction network cutoff and update internal package data.
#
# USAGE (run once, e.g. from the package root):
#
#   source("data-raw/update_network.R")
#
#   # Basic — uses defaults, no plots
#   update_package_network("data-raw/results_PHYSICAL_Prob_Model_16_02_26.csv")
#
#   # Custom pi threshold + save diagnostic plots
#   update_package_network(
#     network_path = "data-raw/results_PHYSICAL_Prob_Model_16_02_26.csv",
#     pi_threshold = 0.7,
#     plot         = TRUE,
#     plot_dir     = "data-raw/plots"
#   )
#
#   # Regenerate plots from a previous result without rerunning the analysis
#   result <- compute_network_cutoff(
#     isoform_network_df = read.csv("data-raw/results_PHYSICAL_Prob_Model_16_02_26.csv"),
#     pi_threshold       = 0.5
#   )
#   plot_cutoff_diagnostics(result, plot_dir = "data-raw/plots")
#
# OUTPUT (written to R/sysdata.rda):
#   .network_data   — isoform-level igraph object
#   .network_cutoff — best ROC q75 cutoff (scalar), defined at gene level
#   .pi_threshold   — pi threshold used to define any_interact
# ==============================================================================

library("dplyr")
library("tidyr")
library("purrr")
library("msigdbr")
library("ggplot2")
library("cowplot")
library("furrr")
library("future")
library("tibble")
library("progressr")
library("usethis")
library("igraph")

handlers(global = TRUE)
handlers("progress")

# ==============================================================================
# INTERNAL HELPERS
# ==============================================================================

.build_gene_pairs <- function(isoform_network_df, pi_threshold) {
  isoform_network_df %>%
    mutate(
      gene_A = if_else(gene_1 < gene_2, gene_1, gene_2),
      gene_B = if_else(gene_1 < gene_2, gene_2, gene_1)
    ) %>%
    group_by(gene_A, gene_B) %>%
    summarise(
      mean_pi         = mean(pi,                           na.rm = TRUE),
      median_pi       = median(pi,                         na.rm = TRUE),
      max_pi          = max(pi,                            na.rm = TRUE),
      n_isoform_pairs = n(),
      n_interact      = sum(pi >= pi_threshold,            na.rm = TRUE),
      any_interact    = as.integer(any(pi >= pi_threshold, na.rm = TRUE)),
      .groups         = "drop"
    )
}

.build_neighbor_list <- function(gene_pairs) {
  bind_rows(
    gene_pairs %>% filter(any_interact == 1) %>% select(gene = gene_A, partner = gene_B),
    gene_pairs %>% filter(any_interact == 1) %>% select(gene = gene_B, partner = gene_A)
  ) %>%
    group_by(gene) %>%
    summarise(neighbors = list(partner), .groups = "drop") %>%
    deframe()
}

.compute_scores_vectorized <- function(pathway_genes, neighbor_list, all_network_genes) {
  pathway_set <- unique(pathway_genes)
  tibble(gene = all_network_genes) %>%
    mutate(
      neighbors         = neighbor_list[gene],
      n_neighbors       = lengths(neighbors),
      n_pathway_overlap = map_int(neighbors, ~ sum(.x %in% pathway_set)),
      ratio             = if_else(n_neighbors > 0,
                                  n_pathway_overlap / n_neighbors,
                                  NA_real_),
      in_pathway        = gene %in% pathway_set
    ) %>%
    filter(!is.na(ratio))
}

.best_cutoff_for_pathway <- function(pathway_name, pathway_genes,
                                     neighbor_list, all_network_genes) {
  scores <- .compute_scores_vectorized(pathway_genes, neighbor_list, all_network_genes)
  n_pos  <- sum( scores$in_pathway)
  n_neg  <- sum(!scores$in_pathway)

  empty <- tibble(
    pathway    = pathway_name, n_pos = n_pos, n_neg = n_neg,
    cutoff_roc = NA_real_,    cutoff_f1 = NA_real_,
    auc_roc    = NA_real_,    auc_pr    = NA_real_,
    f1_best    = NA_real_
  )

  if (n_pos == 0 || n_neg == 0) return(empty)

  thresholds <- sort(unique(c(0, scores$ratio, Inf)), decreasing = TRUE)

  metrics <- map_dfr(thresholds, function(thr) {
    predicted <- scores$ratio >= thr
    tp  <- sum( predicted &  scores$in_pathway)
    fp  <- sum( predicted & !scores$in_pathway)
    fn  <- sum(!predicted &  scores$in_pathway)
    tn  <- sum(!predicted & !scores$in_pathway)
    tpr       <- tp / (tp + fn)
    fpr       <- fp / (fp + tn)
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 1
    recall    <- tpr
    f1        <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
    tibble(threshold = thr, tpr = tpr, fpr = fpr,
           precision = precision, recall = recall,
           f1 = f1, youden_j = tpr - fpr)
  })

  roc_ordered <- metrics %>% arrange(fpr, tpr)
  auc_roc <- sum(diff(roc_ordered$fpr) *
                   (head(roc_ordered$tpr, -1) + tail(roc_ordered$tpr, -1)) / 2)

  pr_ordered <- metrics %>% arrange(recall)
  auc_pr <- if (nrow(pr_ordered) >= 2)
    sum(diff(pr_ordered$recall) *
          (head(pr_ordered$precision, -1) + tail(pr_ordered$precision, -1)) / 2)
  else NA_real_

  best_roc <- metrics %>% slice_max(youden_j, n = 1, with_ties = FALSE)
  best_f1  <- metrics %>% slice_max(f1,       n = 1, with_ties = FALSE)

  tibble(
    pathway    = pathway_name, n_pos = n_pos, n_neg = n_neg,
    cutoff_roc = best_roc$threshold, cutoff_f1 = best_f1$threshold,
    auc_roc    = auc_roc,            auc_pr    = auc_pr,
    f1_best    = best_f1$f1
  )
}

.weighted_quantile <- function(x, w, probs = 0.75) {
  ord <- order(x)
  x_s <- x[ord]
  w_s <- w[ord] / sum(w[ord])
  cdf <- cumsum(w_s)
  approx(cdf, x_s, xout = probs, rule = 2)$y
}

.validate_cutoff <- function(cutoff, all_scores) {
  predicted <- all_scores$ratio >= cutoff
  tp  <- sum( predicted &  all_scores$in_pathway)
  fp  <- sum( predicted & !all_scores$in_pathway)
  fn  <- sum(!predicted &  all_scores$in_pathway)
  tn  <- sum(!predicted & !all_scores$in_pathway)
  tibble(
    cutoff    = cutoff,
    tp = tp, fp = fp, fn = fn, tn = tn,
    precision = tp / (tp + fp),
    recall    = tp / (tp + fn),
    fpr       = fp / (fp + tn),
    f1        = 2 * tp / (2 * tp + fp + fn)
  )
}

# ==============================================================================
# IGRAPH CONSTRUCTION
# ==============================================================================

#' Construct an isoform-level igraph from the raw interaction data frame
#'
#' Edges carry a `pi` weight and a binary `interact` attribute derived from
#' `pi_threshold`. Vertices carry the parent gene symbol.
#'
#' Note: this graph operates at isoform (transcript) level. The cutoff analysis
#' in [compute_network_cutoff()] aggregates to gene level independently.
#'
#' @param isoform_network_df Data frame loaded from the raw interaction CSV.
#'   Expected columns: `enst_1`, `gene_1`, `enst_2`, `gene_2`, `pi`.
#' @param pi_threshold Numeric scalar. Isoform pairs with `pi >= pi_threshold`
#'   are flagged as interacting (`interact = 1`). Must match the value used in
#'   [compute_network_cutoff()].
#'
#' @return An undirected `igraph` object with vertex attribute `gene` and edge
#'   attributes `pi` and `interact`.
.construct_igraph <- function(isoform_network_df, pi_threshold) {
  vertices <- bind_rows(
    isoform_network_df %>% select(enst = enst_1, gene = gene_1),
    isoform_network_df %>% select(enst = enst_2, gene = gene_2)
  ) %>%
    distinct(enst, .keep_all = TRUE)

  edges <- isoform_network_df %>%
    select(from = enst_1, to = enst_2, pi) %>%
    mutate(interact = if_else(pi >= pi_threshold, 1L, 0L))

  graph_from_data_frame(d = edges, vertices = vertices, directed = FALSE)
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION
# ==============================================================================

#' Compute the network-wide interaction cutoff from KEGG pathway benchmarking
#'
#' Accepts an already-loaded data frame — the caller is responsible for reading
#' the CSV and setting the parallel plan before calling this function. This
#' avoids redundant file I/O and double RAM usage when called from
#' [update_package_network()].
#'
#' @param isoform_network_df Data frame loaded from the raw interaction CSV.
#'   Expected columns: `gene_1`, `gene_2`, `pi`.
#' @param pi_threshold Numeric scalar. Minimum `pi` value to call an isoform
#'   pair interacting. Default `0.5`.
#' @param kegg_genesets Optional pre-fetched MSigDB tibble (columns `gs_name`,
#'   `ensembl_gene`, `gs_url`). If `NULL` (default), fetched fresh from MSigDB.
#'   Pass a cached copy during iterative development to avoid repeated downloads.
#' @param plot Logical. Save diagnostic plots to `plot_dir`. Default `FALSE`.
#' @param plot_dir Character. Directory for plots. Created if absent.
#'
#' @return A named list: `gene_pairs`, `gene_degree`, `pathway_results`,
#'   `validation`, `cutoff_roc_summary`, `cutoff_f1_summary`, `best_cutoff`,
#'   `pi_threshold`, `kegg_genesets`.
#'   Note: `all_scores` is NOT stored — it is ~4M rows and only needed
#'   transiently for validation. Use [plot_cutoff_diagnostics()] which
#'   recomputes the subset it needs from `pathway_results` and `gene_pairs`.
compute_network_cutoff <- function(isoform_network_df,
                                   pi_threshold  = 0.5,
                                   kegg_genesets = NULL,
                                   plot          = FALSE,
                                   plot_dir      = "data-raw/plots") {

  cat("================================================================\n")
  cat("  KEGG Pathway Interaction Cutoff Analysis\n")
  cat("================================================================\n\n")

  # --------------------------------------------------------------------------
  # 1. Summarise the input (no file read — data already in memory)
  # --------------------------------------------------------------------------
  cat(">>> [1/8] Inspecting isoform interaction network...\n")
  cat(sprintf("    Rows (isoform pairs)    : %d\n",   nrow(isoform_network_df)))
  cat(sprintf("    Unique isoforms (col 1) : %d\n",   n_distinct(isoform_network_df$gene_1)))
  cat(sprintf("    Unique isoforms (col 2) : %d\n",   n_distinct(isoform_network_df$gene_2)))
  cat(sprintf("    pi threshold            : >= %.3f\n\n", pi_threshold))

  # --------------------------------------------------------------------------
  # 2. Aggregate to gene level
  # --------------------------------------------------------------------------
  cat(">>> [2/8] Aggregating to gene-level pairs...\n")
  gene_pairs    <- .build_gene_pairs(isoform_network_df, pi_threshold)
  n_interacting <- sum(gene_pairs$any_interact == 1)
  cat(sprintf("    Total gene pairs          : %d\n", nrow(gene_pairs)))
  cat(sprintf("    Pairs with any_interact=1 : %d\n", n_interacting))
  cat(sprintf("    Interaction rate          : %.1f%%\n\n",
              100 * n_interacting / nrow(gene_pairs)))

  # --------------------------------------------------------------------------
  # 3. Gene-level degree
  # --------------------------------------------------------------------------
  cat(">>> [3/8] Computing gene-level degree...\n")
  degree_A <- gene_pairs %>%
    group_by(gene = gene_A) %>%
    summarise(degree_binary = sum(any_interact), degree_pi = sum(mean_pi), .groups = "drop")
  degree_B <- gene_pairs %>%
    group_by(gene = gene_B) %>%
    summarise(degree_binary = sum(any_interact), degree_pi = sum(mean_pi), .groups = "drop")
  gene_degree <- bind_rows(degree_A, degree_B) %>%
    group_by(gene) %>%
    summarise(degree_binary = sum(degree_binary), degree_pi = sum(degree_pi), .groups = "drop")
  cat(sprintf("    Genes in network       : %d\n",     nrow(gene_degree)))
  cat(sprintf("    Median binary degree   : %.0f\n",   median(gene_degree$degree_binary)))
  cat(sprintf("    Max binary degree      : %.0f\n\n", max(gene_degree$degree_binary)))

  # --------------------------------------------------------------------------
  # 4. KEGG gene sets (fetch once, or accept pre-cached copy)
  # --------------------------------------------------------------------------
  cat(">>> [4/8] Loading KEGG gene sets...\n")
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
  cat(sprintf("    Total unique pathways : %d\n",   n_pathways))
  cat(sprintf("    Total unique genes    : %d\n\n", n_kegg_genes))

  # --------------------------------------------------------------------------
  # 5. Neighbor lookup
  # --------------------------------------------------------------------------
  cat(">>> [5/8] Building neighbor lookup table...\n")
  neighbor_list     <- .build_neighbor_list(gene_pairs)
  all_network_genes <- names(neighbor_list)
  n_overlap         <- length(intersect(kegg_genesets$ensembl_gene, all_network_genes))
  cat(sprintf("    Genes with >= 1 interaction : %d\n",        length(all_network_genes)))
  cat(sprintf("    KEGG genes in network       : %d / %d (%.1f%%)\n\n",
              n_overlap, n_kegg_genes, 100 * n_overlap / n_kegg_genes))

  # --------------------------------------------------------------------------
  # 6. Per-pathway ROC + PR analysis
  #    Uses whatever plan() is active — set by the caller (update_package_network)
  # --------------------------------------------------------------------------
  cat(">>> [6/8] Running per-pathway ROC + PR analysis...\n")
  cat(sprintf("    Scanning %d pathways\n\n", n_pathways))

  pathway_list <- kegg_genesets %>% group_by(gs_name) %>% group_split()

  with_progress({
    p <- progressor(along = pathway_list)
    pathway_results <- future_map_dfr(pathway_list, function(pathway_df) {
      result <- .best_cutoff_for_pathway(
        pathway_name      = unique(pathway_df$gs_name),
        pathway_genes     = pathway_df$ensembl_gene,
        neighbor_list     = neighbor_list,
        all_network_genes = all_network_genes
      )
      p()
      result
    }, .options = furrr_options(seed = TRUE))
  })

  n_valid   <- sum(!is.na(pathway_results$cutoff_roc))
  n_skipped <- sum( is.na(pathway_results$cutoff_roc))
  cat(sprintf("\n    Valid pathways : %d | Skipped (no pos or neg): %d\n", n_valid, n_skipped))
  cat(sprintf("    Median ROC-AUC : %.3f\n", median(pathway_results$auc_roc, na.rm = TRUE)))
  cat(sprintf("    Median PR-AUC  : %.3f\n", median(pathway_results$auc_pr,  na.rm = TRUE)))
  cat(sprintf("    ROC-AUC > 0.7  : %d pathways\n",
              sum(pathway_results$auc_roc > 0.7, na.rm = TRUE)))
  cat(sprintf("    PR-AUC  > 0.7  : %d pathways\n\n",
              sum(pathway_results$auc_pr  > 0.7, na.rm = TRUE)))

  # --------------------------------------------------------------------------
  # 7. AUC-weighted cutoff summaries
  # --------------------------------------------------------------------------
  cat(">>> [7/8] Computing AUC-weighted conservative cutoffs...\n")

  cutoff_roc_summary <- pathway_results %>%
    filter(!is.na(cutoff_roc), auc_roc > 0.5) %>%
    mutate(w = auc_roc - 0.5) %>%
    summarise(
      n_pathways_used = n(),
      q50             = .weighted_quantile(cutoff_roc, w, 0.50),
      q75             = .weighted_quantile(cutoff_roc, w, 0.75),
      q90             = .weighted_quantile(cutoff_roc, w, 0.90),
      mean_auc_roc    = weighted.mean(auc_roc, w)
    )

  cutoff_f1_summary <- pathway_results %>%
    filter(!is.na(cutoff_f1), !is.na(auc_pr)) %>%
    mutate(
      prevalence = n_pos / (n_pos + n_neg),
      w          = pmax(auc_pr - prevalence, 0)
    ) %>%
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
  cat(sprintf("    Cutoff q50 — ROC: %.4f | PR/F1: %.4f\n",
              cutoff_roc_summary$q50, cutoff_f1_summary$q50))
  cat(sprintf("    Cutoff q75 — ROC: %.4f | PR/F1: %.4f  <-- recommended\n",
              cutoff_roc_summary$q75, cutoff_f1_summary$q75))
  cat(sprintf("    Cutoff q90 — ROC: %.4f | PR/F1: %.4f\n\n",
              cutoff_roc_summary$q90, cutoff_f1_summary$q90))

  # --------------------------------------------------------------------------
  # 8. Pooled validation
  #    all_scores is ~4M rows — compute transiently and discard after validation
  # --------------------------------------------------------------------------
  cat(">>> [8/8] Pooled validation across all gene x pathway observations...\n")

  all_scores <- future_map_dfr(pathway_list, function(pathway_df) {
    .compute_scores_vectorized(
      pathway_df$ensembl_gene, neighbor_list, all_network_genes
    ) %>%
      mutate(pathway = unique(pathway_df$gs_name))
  }, .options = furrr_options(seed = TRUE))

  validation <- bind_rows(
    .validate_cutoff(cutoff_roc_summary$q50, all_scores) %>% mutate(method = "ROC",   quantile = "q50"),
    .validate_cutoff(cutoff_roc_summary$q75, all_scores) %>% mutate(method = "ROC",   quantile = "q75"),
    .validate_cutoff(cutoff_roc_summary$q90, all_scores) %>% mutate(method = "ROC",   quantile = "q90"),
    .validate_cutoff(cutoff_f1_summary$q50,  all_scores) %>% mutate(method = "PR/F1", quantile = "q50"),
    .validate_cutoff(cutoff_f1_summary$q75,  all_scores) %>% mutate(method = "PR/F1", quantile = "q75"),
    .validate_cutoff(cutoff_f1_summary$q90,  all_scores) %>% mutate(method = "PR/F1", quantile = "q90")
  ) %>%
    mutate(
      label = paste0(method, " | ",
                     case_when(quantile == "q50" ~ "q50 — median",
                               quantile == "q75" ~ "q75 — conservative",
                               quantile == "q90" ~ "q90 — very strict")),
      .before = 1
    )

  # Discard all_scores — too large to carry in the result list
  rm(all_scores)

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
  cat("  Recommended network cutoffs (q75):\n")
  cat(sprintf("    ROC-based   : %.4f  (precision=%.3f, recall=%.3f, FPR=%.3f)\n",
              cutoff_roc_summary$q75,
              validation$precision[validation$method == "ROC"   & validation$quantile == "q75"],
              validation$recall[   validation$method == "ROC"   & validation$quantile == "q75"],
              validation$fpr[      validation$method == "ROC"   & validation$quantile == "q75"]))
  cat(sprintf("    PR/F1-based : %.4f  (precision=%.3f, recall=%.3f, FPR=%.3f)\n",
              cutoff_f1_summary$q75,
              validation$precision[validation$method == "PR/F1" & validation$quantile == "q75"],
              validation$recall[   validation$method == "PR/F1" & validation$quantile == "q75"],
              validation$fpr[      validation$method == "PR/F1" & validation$quantile == "q75"]))
  cat(sprintf("\n  >>> Returning best cutoff (ROC q75): %.4f\n", best_cutoff))
  cat("================================================================\n\n")

  result <- list(
    gene_pairs         = gene_pairs,
    gene_degree        = gene_degree,
    pathway_results    = pathway_results,
    validation         = validation,
    cutoff_roc_summary = cutoff_roc_summary,
    cutoff_f1_summary  = cutoff_f1_summary,
    best_cutoff        = best_cutoff,
    pi_threshold       = pi_threshold,
    kegg_genesets      = kegg_genesets
  )

  if (plot) plot_cutoff_diagnostics(result, plot_dir = plot_dir)

  return(result)
}

# ==============================================================================
# PLOTTING FUNCTION
# ==============================================================================

#' Generate diagnostic plots from a [compute_network_cutoff()] result
#'
#' Can be called independently on any saved result list without rerunning the
#' full analysis. Recomputes the gene-level score summaries it needs directly
#' from `pathway_results` — does not require `all_scores` to be stored.
#'
#' @param result   List returned by [compute_network_cutoff()].
#' @param plot_dir Character. Directory to save PNG files (created if needed).
plot_cutoff_diagnostics <- function(result, plot_dir = "data-raw/plots") {

  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  pathway_results    <- result$pathway_results
  cutoff_roc_summary <- result$cutoff_roc_summary
  cutoff_f1_summary  <- result$cutoff_f1_summary
  kegg_genesets      <- result$kegg_genesets
  neighbor_list      <- .build_neighbor_list(result$gene_pairs)
  all_network_genes  <- names(neighbor_list)

  cat("  Generating diagnostic plots...\n")

  # --- Plot 1: AUC-weighted distribution of per-pathway best cutoffs ---------
  roc_vlines <- tibble(
    label  = c("q50", "q75", "q90"),
    cutoff = c(cutoff_roc_summary$q50, cutoff_roc_summary$q75, cutoff_roc_summary$q90)
  )
  f1_vlines <- tibble(
    label  = c("q50", "q75", "q90"),
    cutoff = c(cutoff_f1_summary$q50, cutoff_f1_summary$q75, cutoff_f1_summary$q90)
  )

  p1a <- pathway_results %>%
    filter(!is.na(cutoff_roc), auc_roc > 0.5) %>%
    mutate(w = auc_roc - 0.5) %>%
    ggplot(aes(x = cutoff_roc, weight = w)) +
    geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
    geom_vline(data = roc_vlines, aes(xintercept = cutoff, colour = label),
               linewidth = 0.9, linetype = "dashed") +
    scale_colour_brewer(palette = "Set1", name = "Quantile") +
    labs(title    = "ROC-based cutoffs",
         subtitle = "Weighted by (AUC_ROC - 0.5)",
         x = "Best cutoff (Youden's J)", y = "Weighted count") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(colour = "grey40"))

  p1b <- pathway_results %>%
    filter(!is.na(cutoff_f1), !is.na(auc_pr)) %>%
    mutate(prevalence = n_pos / (n_pos + n_neg),
           w          = pmax(auc_pr - prevalence, 0)) %>%
    filter(w > 0) %>%
    ggplot(aes(x = cutoff_f1, weight = w)) +
    geom_histogram(bins = 40, fill = "darkorange", alpha = 0.7) +
    geom_vline(data = f1_vlines, aes(xintercept = cutoff, colour = label),
               linewidth = 0.9, linetype = "dashed") +
    scale_colour_brewer(palette = "Set1", name = "Quantile") +
    labs(title    = "PR/F1-based cutoffs",
         subtitle = "Weighted by (AUC_PR - prevalence)",
         x = "Best cutoff (max F1)", y = "Weighted count") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(colour = "grey40"))

  p1 <- cowplot::plot_grid(p1a, p1b, ncol = 2)

  # --- Plot 2: Cutoff vs number of pathway memberships ----------------------
  gene_pathway_count <- kegg_genesets %>%
    group_by(ensembl_gene) %>%
    summarise(n_pathways = n_distinct(gs_name), .groups = "drop")

  # Recompute only the per-gene summary needed for this plot —
  # much smaller than storing all_scores (~n_kegg_genes rows vs ~4M)
  gene_cutoff_summary <- kegg_genesets %>%
    distinct(ensembl_gene, gs_name) %>%
    inner_join(pathway_results %>% select(pathway = gs_name, cutoff_roc, auc_roc),
               by = c("gs_name" = "pathway"),
               relationship = "many-to-many") %>%
    filter(ensembl_gene %in% all_network_genes) %>%
    group_by(gene = ensembl_gene) %>%
    summarise(
      auc_weighted_cutoff    = weighted.mean(cutoff_roc,
                                             w = pmax(auc_roc - 0.5, 0),
                                             na.rm = TRUE),
      n_pathways             = n_distinct(gs_name),
      n_pathways_with_cutoff = sum(!is.na(cutoff_roc)),
      .groups = "drop"
    ) %>%
    filter(!is.na(auc_weighted_cutoff)) %>%
    mutate(pathway_group = cut(n_pathways,
                               breaks = c(0, 1, 2, 5, 10, 20, Inf),
                               labels = c("1", "2", "3-5", "6-10", "11-20", "20+"),
                               right  = TRUE))

  group_medians <- gene_cutoff_summary %>%
    group_by(pathway_group) %>%
    summarise(median_cutoff = median(auc_weighted_cutoff, na.rm = TRUE),
              n_genes = n(), .groups = "drop")

  p2 <- ggplot(gene_cutoff_summary, aes(x = pathway_group, y = auc_weighted_cutoff)) +
    geom_boxplot(aes(fill = pathway_group), outlier.shape = NA, alpha = 0.6, width = 0.55) +
    geom_jitter(aes(colour = pathway_group), width = 0.2, size = 0.8, alpha = 0.4) +
    geom_text(data = group_medians,
              aes(y = median_cutoff, label = round(median_cutoff, 2)),
              vjust = -0.6, size = 3.2, fontface = "bold") +
    geom_text(data = group_medians,
              aes(y = -Inf, label = paste0("n=", n_genes)),
              vjust = -0.5, size = 3, colour = "grey40") +
    geom_hline(yintercept = cutoff_roc_summary$q75,
               linetype = "dashed", colour = "firebrick", linewidth = 0.8) +
    annotate("text", x = Inf, y = cutoff_roc_summary$q75,
             label = paste0("ROC q75: ", round(cutoff_roc_summary$q75, 3)),
             hjust = 1.05, vjust = -0.5, colour = "firebrick", size = 3.2) +
    scale_fill_brewer(palette = "Blues", guide = "none") +
    scale_colour_brewer(palette = "Blues", guide = "none") +
    scale_x_discrete(drop = FALSE) +
    labs(title    = "Best cutoff by number of pathway memberships",
         subtitle = "AUC-weighted ROC cutoff per gene; red line = network threshold (q75)",
         x        = "Number of KEGG pathways gene belongs to",
         y        = "AUC-weighted best cutoff (ROC)",
         caption  = "n = number of genes in each group") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          plot.subtitle = element_text(colour = "grey40"),
          plot.caption  = element_text(colour = "grey40"))

  out1 <- file.path(plot_dir, "plot_cutoff_distribution.png")
  out2 <- file.path(plot_dir, "plot_cutoff_by_pathway_count.png")
  ggsave(out1, p1, width = 14, height = 5, dpi = 150)
  ggsave(out2, p2, width = 10, height = 6, dpi = 150)
  cat(sprintf("    Saved: %s\n", out1))
  cat(sprintf("    Saved: %s\n", out2))
}

# ==============================================================================
# PACKAGE UPDATE FUNCTION
# ==============================================================================

#' Run cutoff analysis and write results to internal package data (sysdata.rda)
#'
#' Reads the CSV once, builds the igraph, runs [compute_network_cutoff()], then
#' writes three internal objects to `R/sysdata.rda`.
#'
#' The parallel plan is set here and is active for all `future_map_dfr` calls
#' inside [compute_network_cutoff()]. Workers are spawned once and reused across
#' both the pathway ROC loop (step 6) and the pooled validation (step 8).
#'
#' @param network_path Character. Path to the raw isoform interaction CSV.
#' @param pi_threshold Numeric scalar. Default `0.5`.
#' @param kegg_genesets Optional pre-fetched MSigDB tibble. Passed to
#'   [compute_network_cutoff()] to avoid repeated downloads.
#' @param plot Logical. Save diagnostic plots. Default `FALSE`.
#' @param plot_dir Character. Directory for plots.
#' @param n_workers Integer. Parallel workers. Defaults to all cores minus one.
#'
#' @return Invisibly returns the result list from [compute_network_cutoff()].
update_package_network <- function(network_path,
                                   pi_threshold  = 0.5,
                                   kegg_genesets = NULL,
                                   plot          = FALSE,
                                   plot_dir      = "data-raw/plots",
                                   n_workers     = parallel::detectCores() - 1) {

  # Set plan once — covers both future_map_dfr calls inside compute_network_cutoff()
  # Workers are spawned here and reused; on.exit ensures cleanup even on error
  plan(multisession, workers = n_workers)
  on.exit(plan(sequential), add = TRUE)
  cat(sprintf(">>> Parallel plan: multisession with %d workers\n\n", n_workers))

  # Read CSV once — passed to both .construct_igraph() and compute_network_cutoff()
  # This is the only read_csv call in the entire pipeline
  cat(">>> Reading raw network CSV...\n")
  isoform_network_df <- read.csv(network_path, encoding = "utf-8")
  cat(sprintf("    Read %d rows from %s\n\n", nrow(isoform_network_df), network_path))

  # Build isoform-level igraph (uses the same data frame, no second read)
  cat(">>> Building isoform-level igraph...\n")
  .network_data <- .construct_igraph(isoform_network_df, pi_threshold)
  cat(sprintf("    igraph: %d vertices, %d edges\n\n",
              vcount(.network_data), ecount(.network_data)))

  # Run full gene-level cutoff analysis (data frame passed directly, no re-read)
  result <- compute_network_cutoff(
    isoform_network_df = isoform_network_df,
    pi_threshold       = pi_threshold,
    kegg_genesets      = kegg_genesets,
    plot               = plot,
    plot_dir           = plot_dir
  )

  .network_cutoff <- result$best_cutoff
  .pi_threshold   <- result$pi_threshold

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
  cat(sprintf("    .network_cutoff : %.4f (gene-level ROC q75)\n", .network_cutoff))
  cat(sprintf("    .pi_threshold   : %.3f\n",                       .pi_threshold))
  cat("    Done. Rebuild the package to make the new data active.\n\n")

  invisible(result)
}
