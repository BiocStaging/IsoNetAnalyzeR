# ==============================================================================
# data-raw/network_helpers.R
#
# Pure data-manipulation helpers used by compute_network_cutoff() and the
# plotting functions.  No file I/O, no parallelism, no ggplot2.
#
# All functions are prefixed with "." and are NOT exported.
#
# Contents
# --------
#   Data aggregation
#     .build_gene_pairs             — isoform pairs → gene-level summary
#     .build_neighbor_list          — gene pairs → named interaction list
#     .compute_scores_vectorized    — pathway gene set → per-gene ratio scores
#
#   Cutoff analysis
#     .weighted_quantile            — AUC-weighted quantile
#     .validate_cutoff              — confusion matrix metrics at one cutoff
#     .derive_sens_spec_cutoffs     — three direct pooled cutoffs from sens/spec
#
#   Multi-pathway gene strategies
#     .apply_multipathway_strategy  — deduplicate gene-set under a strategy
#     .run_strategy                 — full analysis under one strategy
#
#   igraph construction
#     .construct_igraph             — raw data frame → isoform-level igraph
#
#   Per-pathway full curve computation (used by plot_roc_curves / plot_pr_curves)
#     .full_curves_for_pathway      — complete threshold sweep for one pathway
#     .compute_aggregate_curves     — pooled threshold sweep over all pathways
# ==============================================================================

# Shared furrr options: ensures each parallel worker loads the required
# packages. Workers start with a clean R environment and do NOT inherit the
# library() calls made in the main session.
.furrr_pkgs <- c("dplyr", "tidyr", "purrr", "tibble", "progressr")


# ==============================================================================
# DATA AGGREGATION
# ==============================================================================

#' Aggregate isoform-level pairs to gene-level summary
#'
#' For each unordered gene pair (gene_A ≤ gene_B alphabetically) computes
#' per-pair pi summaries and a binary interaction flag.
#'
#' @param isoform_network_df Data frame with columns \code{gene_1}, \code{gene_2},
#'   \code{pi}.
#' @param pi_threshold Numeric scalar. Pairs with \code{pi >= pi_threshold} count
#'   as interacting.
#'
#' @return A tibble with one row per gene pair and columns:
#'   \code{gene_A}, \code{gene_B}, \code{mean_pi}, \code{median_pi},
#'   \code{max_pi}, \code{n_isoform_pairs}, \code{n_interact},
#'   \code{any_interact}.
#' @noRd
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


#' Build a named neighbor list from gene-level pairs
#'
#' Returns a named list where each element is the character vector of a gene's
#' interaction partners (only edges with \code{any_interact == 1}).
#'
#' @param gene_pairs Tibble produced by \code{.build_gene_pairs()}.
#' @return Named list: gene → character vector of partner gene IDs.
#' @noRd
.build_neighbor_list <- function(gene_pairs) {
  bind_rows(
    gene_pairs %>% filter(any_interact == 1) %>% select(gene = gene_A, partner = gene_B),
    gene_pairs %>% filter(any_interact == 1) %>% select(gene = gene_B, partner = gene_A)
  ) %>%
    group_by(gene) %>%
    summarise(neighbors = list(partner), .groups = "drop") %>%
    deframe()
}


#' Compute per-gene ratio scores for a single pathway
#'
#' For each gene in the network, the ratio is the fraction of its interaction
#' partners that belong to the given pathway gene set.
#'
#' @param pathway_genes Character vector of gene IDs belonging to the pathway.
#' @param neighbor_list Named list returned by \code{.build_neighbor_list()}.
#' @param all_network_genes Character vector of genes with ≥ 1 interaction.
#'
#' @return A tibble with columns \code{gene}, \code{neighbors},
#'   \code{n_neighbors}, \code{n_pathway_overlap}, \code{ratio},
#'   \code{in_pathway}. Rows with \code{NA} ratio are dropped.
#' @noRd
.compute_scores_vectorized <- function(pathway_genes, neighbor_list,
                                       all_network_genes) {
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


# ==============================================================================
# CUTOFF ANALYSIS
# ==============================================================================

#' AUC-weighted quantile
#'
#' @param x Numeric vector of values.
#' @param w Numeric vector of weights (same length as \code{x}).
#' @param probs Numeric scalar in \[0, 1\]. Default \code{0.75}.
#' @return Numeric scalar — the weighted quantile.
#' @noRd
.weighted_quantile <- function(x, w, probs = 0.75) {
  ord <- order(x)
  x_s <- x[ord]
  w_s <- w[ord] / sum(w[ord])
  cdf <- cumsum(w_s)
  approx(cdf, x_s, xout = probs, rule = 2)$y
}


#' Validate a single cutoff against pooled scores
#'
#' Computes the full confusion-matrix metrics (including sensitivity and
#' specificity) for one candidate cutoff applied to the pooled
#' gene × pathway score table.
#'
#' @param cutoff Numeric scalar. The ratio threshold to evaluate.
#' @param all_scores Tibble with columns \code{ratio} and \code{in_pathway}
#'   (as produced by \code{.compute_scores_vectorized()}).
#'
#' @return A one-row tibble with columns: \code{cutoff}, \code{tp}, \code{fp},
#'   \code{fn}, \code{tn}, \code{sensitivity}, \code{specificity},
#'   \code{precision}, \code{recall}, \code{fpr}, \code{f1}.
#' @noRd
.validate_cutoff <- function(cutoff, all_scores) {
  predicted <- all_scores$ratio >= cutoff
  tp  <- sum( predicted &  all_scores$in_pathway)
  fp  <- sum( predicted & !all_scores$in_pathway)
  fn  <- sum(!predicted &  all_scores$in_pathway)
  tn  <- sum(!predicted & !all_scores$in_pathway)

  sensitivity <- tp / (tp + fn)   # = recall = TPR
  specificity <- tn / (tn + fp)   # = TNR = 1 - FPR

  tibble(
    cutoff      = cutoff,
    tp = tp, fp = fp, fn = fn, tn = tn,
    sensitivity = sensitivity,
    specificity = specificity,
    precision   = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
    recall      = sensitivity,       # duplicate for API compatibility
    fpr         = 1 - specificity,
    f1          = 2 * tp / (2 * tp + fp + fn)
  )
}


#' Derive three cutoffs directly from pooled sensitivity and specificity
#'
#' Sweeps every unique ratio threshold in \code{all_scores} and computes
#' pooled sensitivity (TPR) and specificity (TNR) at each point.
#' Returns three operating points optimising different criteria:
#'
#' \describe{
#'   \item{Pooled Youden's J}{Maximises \code{sensitivity + specificity − 1}.}
#'   \item{Sensitivity = Specificity}{Minimises \code{|sensitivity − specificity|}.}
#'   \item{Min distance to (0, 1)}{Minimises \eqn{\sqrt{(1-sens)^2 + (1-spec)^2}}.}
#' }
#'
#' @param all_scores Pooled scores tibble (\code{ratio}, \code{in_pathway}).
#' @return A three-row tibble with columns: \code{method}, \code{cutoff},
#'   \code{sensitivity}, \code{specificity}, \code{precision}, \code{f1},
#'   \code{youden_j}, \code{dist_topleft}.
#' @noRd
.derive_sens_spec_cutoffs <- function(all_scores) {

  ratio      <- all_scores$ratio
  in_pathway <- all_scores$in_pathway

  thresholds  <- sort(unique(ratio))
  n_obs       <- length(ratio)

  ord       <- order(ratio)
  ratio_s   <- ratio[ord]
  pathway_s <- in_pathway[ord]

  n_pos_total <- sum(in_pathway)
  n_neg_total <- n_obs - n_pos_total

  cum_tp_asc <- cumsum( pathway_s)
  cum_fp_asc <- cumsum(!pathway_s)

  n_below <- findInterval(thresholds, ratio_s, left.open = TRUE)

  tp <- ifelse(n_below == 0L,
               n_pos_total,
               n_pos_total - cum_tp_asc[n_below])
  fp <- ifelse(n_below == 0L,
               n_neg_total,
               n_neg_total - cum_fp_asc[n_below])
  fn <- n_pos_total - tp
  tn <- n_neg_total - fp

  sens      <- tp / (tp + fn)
  spec      <- tn / (tn + fp)
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
  f1        <- 2 * tp / (2 * tp + fp + fn)

  metrics <- tibble(
    threshold      = thresholds,
    sensitivity    = sens,
    specificity    = spec,
    precision      = precision,
    f1             = f1,
    youden_j       = sens + spec - 1,
    dist_topleft   = sqrt((1 - sens)^2 + (1 - spec)^2),
    sens_spec_diff = abs(sens - spec)
  )

  bind_rows(
    metrics %>%
      slice_max(youden_j,       n = 1, with_ties = FALSE) %>%
      mutate(method = "Pooled Youden's J"),
    metrics %>%
      slice_min(sens_spec_diff, n = 1, with_ties = FALSE) %>%
      mutate(method = "Sensitivity = Specificity"),
    metrics %>%
      slice_min(dist_topleft,   n = 1, with_ties = FALSE) %>%
      mutate(method = "Min distance to (0,1)")
  ) %>%
    select(method, cutoff = threshold, sensitivity, specificity,
           precision, f1, youden_j, dist_topleft)
}


# ==============================================================================
# PER-PATHWAY ROC / PR ANALYSIS
# ==============================================================================

#' Per-pathway ROC + PR analysis — returns best cutoffs and key metrics
#'
#' Runs the full threshold sweep for one pathway and returns:
#' AUC-ROC, AUC-PR, best cutoffs (Youden's J and max F1), and sens/spec
#' at the best Youden's J operating point.
#'
#' @param pathway_name  Character scalar.
#' @param pathway_genes Character vector of pathway gene IDs.
#' @param neighbor_list Named interaction list.
#' @param all_network_genes Character vector.
#'
#' @return A one-row tibble with columns: \code{pathway}, \code{n_pos},
#'   \code{n_neg}, \code{cutoff_roc}, \code{cutoff_f1}, \code{auc_roc},
#'   \code{auc_pr}, \code{f1_best}, \code{sensitivity}, \code{specificity}.
#' @noRd
.best_cutoff_for_pathway <- function(pathway_name, pathway_genes,
                                     neighbor_list, all_network_genes) {
  scores <- .compute_scores_vectorized(pathway_genes, neighbor_list,
                                       all_network_genes)
  n_pos  <- sum( scores$in_pathway)
  n_neg  <- sum(!scores$in_pathway)

  empty <- tibble(
    pathway     = pathway_name, n_pos = n_pos, n_neg = n_neg,
    cutoff_roc  = NA_real_,    cutoff_f1   = NA_real_,
    auc_roc     = NA_real_,    auc_pr      = NA_real_,
    f1_best     = NA_real_,
    sensitivity = NA_real_,    specificity = NA_real_
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
    tnr       <- tn / (tn + fp)
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 1
    recall    <- tpr
    f1        <- if ((precision + recall) > 0)
      2 * precision * recall / (precision + recall) else 0

    tibble(threshold = thr, tpr = tpr, fpr = fpr, tnr = tnr,
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
    pathway     = pathway_name, n_pos = n_pos, n_neg = n_neg,
    cutoff_roc  = best_roc$threshold, cutoff_f1 = best_f1$threshold,
    auc_roc     = auc_roc,            auc_pr    = auc_pr,
    f1_best     = best_f1$f1,
    sensitivity = best_roc$tpr,
    specificity = best_roc$tnr
  )
}


# ==============================================================================
# MULTI-PATHWAY GENE HANDLING
# ==============================================================================

#' Apply a multi-pathway gene deduplication strategy
#'
#' @param kegg_genesets  Raw MSigDB tibble with columns \code{gs_name},
#'   \code{ensembl_gene}.
#' @param neighbor_list  Named interaction list.
#' @param strategy       One of \code{"all"}, \code{"remove"},
#'   \code{"random"}, \code{"best"}.
#' @param seed           Integer seed used when \code{strategy = "random"}.
#'
#' @return Tibble with columns \code{gs_name}, \code{ensembl_gene}.
#' @noRd
.apply_multipathway_strategy <- function(kegg_genesets,
                                         neighbor_list,
                                         strategy = c("all", "remove",
                                                      "random", "best"),
                                         seed = 42L) {
  strategy <- match.arg(strategy)

  gene_pathway_counts <- kegg_genesets %>%
    group_by(ensembl_gene) %>%
    summarise(n_pathways = n_distinct(gs_name), .groups = "drop")

  multi_genes <- gene_pathway_counts %>%
    filter(n_pathways > 1) %>%
    pull(ensembl_gene)

  if (strategy == "all")
    return(kegg_genesets %>% select(gs_name, ensembl_gene))

  if (strategy == "remove")
    return(kegg_genesets %>%
             filter(!ensembl_gene %in% multi_genes) %>%
             select(gs_name, ensembl_gene))

  single_gene_rows <- kegg_genesets %>%
    filter(!ensembl_gene %in% multi_genes) %>%
    select(gs_name, ensembl_gene)

  multi_gene_rows <- kegg_genesets %>%
    filter(ensembl_gene %in% multi_genes) %>%
    select(gs_name, ensembl_gene)

  if (strategy == "random") {
    set.seed(seed)
    assigned <- multi_gene_rows %>%
      group_by(ensembl_gene) %>%
      slice_sample(n = 1) %>%
      ungroup()
    return(bind_rows(single_gene_rows, assigned))
  }

  # strategy == "best": assign to pathway with most within-pathway interactions
  all_pathway_genes_by_name <- kegg_genesets %>%
    group_by(gs_name) %>%
    summarise(members = list(ensembl_gene), .groups = "drop") %>%
    deframe()

  assigned <- multi_gene_rows %>%
    distinct(ensembl_gene, gs_name) %>%
    mutate(
      n_within = map2_int(ensembl_gene, gs_name, function(g, pw) {
        nbrs    <- neighbor_list[[g]]
        pw_mems <- all_pathway_genes_by_name[[pw]]
        if (is.null(nbrs) || is.null(pw_mems)) return(0L)
        sum(nbrs %in% pw_mems)
      })
    ) %>%
    group_by(ensembl_gene) %>%
    arrange(desc(n_within), gs_name) %>%
    slice(1) %>%
    ungroup() %>%
    select(gs_name, ensembl_gene)

  bind_rows(single_gene_rows, assigned)
}


#' Run a complete pathway analysis under one multi-pathway strategy
#'
#' Applies the strategy, re-runs per-pathway analysis, and returns cutoff
#' summaries + validation metrics compatible with the main result list.
#'
#' @param kegg_genesets  Raw MSigDB tibble.
#' @param neighbor_list  Named interaction list.
#' @param all_network_genes Character vector of genes with ≥ 1 interaction.
#' @param strategy       Strategy string.
#' @param all_scores_ref Pooled scores tibble from the baseline "all" run,
#'   used so validation metrics are computed on a common denominator.
#'
#' @return Named list: \code{strategy}, \code{cutoff_roc_summary},
#'   \code{cutoff_f1_summary}, \code{validation}, \code{n_genes_used},
#'   \code{n_pathways_used}, \code{pathway_results}.
#' @noRd
.run_strategy <- function(kegg_genesets, neighbor_list, all_network_genes,
                          strategy, all_scores_ref) {

  cat(sprintf("    [%s] Applying strategy and running per-pathway analysis...\n",
              strategy))

  deduped_sets <- .apply_multipathway_strategy(kegg_genesets, neighbor_list,
                                               strategy = strategy)
  pathway_list <- deduped_sets %>% group_by(gs_name) %>% group_split()

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

  validation <- bind_rows(
    .validate_cutoff(cutoff_roc_summary$q75, all_scores_ref) %>%
      mutate(method = "ROC",   quantile = "q75"),
    .validate_cutoff(cutoff_f1_summary$q75,  all_scores_ref) %>%
      mutate(method = "PR/F1", quantile = "q75")
  ) %>%
    mutate(strategy = strategy, .before = 1)

  list(
    strategy           = strategy,
    cutoff_roc_summary = cutoff_roc_summary,
    cutoff_f1_summary  = cutoff_f1_summary,
    validation         = validation,
    n_genes_used       = n_distinct(deduped_sets$ensembl_gene),
    n_pathways_used    = nrow(pathway_results),
    pathway_results    = pathway_results
  )
}


# ==============================================================================
# IGRAPH CONSTRUCTION
# ==============================================================================

#' Construct an isoform-level igraph from the raw interaction data frame
#'
#' Edges carry a \code{pi} weight and a binary \code{interact} attribute
#' derived from \code{pi_threshold}. Vertices carry the parent gene symbol.
#'
#' @param isoform_network_df Data frame with columns \code{enst_1},
#'   \code{gene_1}, \code{enst_2}, \code{gene_2}, \code{pi}.
#' @param pi_threshold Numeric scalar. Pairs with \code{pi >= pi_threshold} are
#'   flagged as interacting (\code{interact = 1}).
#'
#' @return An undirected \code{igraph} with vertex attribute \code{gene} and
#'   edge attributes \code{pi} and \code{interact}.
#' @noRd
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
# PER-PATHWAY FULL CURVE COMPUTATION (used by network_plots.R)
# ==============================================================================

#' Return the complete ROC / PR curve points for one pathway
#'
#' Mirrors the threshold sweep in \code{.best_cutoff_for_pathway()} but keeps
#' all rows (rather than only the best-threshold summary). Used by
#' \code{plot_roc_curves()} and \code{plot_pr_curves()} to draw spaghetti
#' background curves.
#'
#' @param pathway_name  Character scalar.
#' @param pathway_genes Character vector.
#' @param neighbor_list Named interaction list.
#' @param all_network_genes Character vector.
#'
#' @return A tibble with columns \code{pathway}, \code{threshold}, \code{tpr},
#'   \code{fpr}, \code{precision}, \code{recall}, \code{n_pos}, \code{n_neg},
#'   \code{auc_roc}, \code{auc_pr}. Returns \code{NULL} if no valid
#'   positives/negatives.
#' @noRd
.full_curves_for_pathway <- function(pathway_name, pathway_genes,
                                     neighbor_list, all_network_genes) {

  scores <- .compute_scores_vectorized(pathway_genes, neighbor_list,
                                       all_network_genes)

  n_pos <- sum( scores$in_pathway)
  n_neg <- sum(!scores$in_pathway)

  if (n_pos == 0 || n_neg == 0) return(NULL)

  thresholds <- sort(unique(c(0, scores$ratio, Inf)), decreasing = TRUE)

  curve_df <- purrr::map_dfr(thresholds, function(thr) {
    predicted <- scores$ratio >= thr
    tp  <- sum( predicted &  scores$in_pathway)
    fp  <- sum( predicted & !scores$in_pathway)
    fn  <- sum(!predicted &  scores$in_pathway)
    tn  <- sum(!predicted & !scores$in_pathway)

    tpr       <- tp / (tp + fn)
    fpr       <- fp / (fp + tn)
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 1
    recall    <- tpr

    tibble::tibble(
      pathway   = pathway_name,
      threshold = thr,
      tpr       = tpr,
      fpr       = fpr,
      precision = precision,
      recall    = recall,
      n_pos     = n_pos,
      n_neg     = n_neg
    )
  })

  roc_ord <- curve_df %>% dplyr::arrange(fpr, tpr)
  auc_roc  <- sum(diff(roc_ord$fpr) *
                    (utils::head(roc_ord$tpr, -1) + utils::tail(roc_ord$tpr, -1)) / 2)

  pr_ord <- curve_df %>% dplyr::arrange(recall)
  auc_pr  <- if (nrow(pr_ord) >= 2)
    sum(diff(pr_ord$recall) *
          (utils::head(pr_ord$precision, -1) + utils::tail(pr_ord$precision, -1)) / 2)
  else NA_real_

  curve_df %>% dplyr::mutate(auc_roc = auc_roc, auc_pr = auc_pr)
}


#' Pool scores across all pathways and sweep thresholds to produce a single
#' aggregate ROC/PR curve
#'
#' Exactly what \code{.validate_cutoff()} does for a single threshold, but
#' swept over all unique ratio values to form a complete curve.
#'
#' @param pathway_list List of per-pathway tibbles (from \code{group_split()}).
#' @param neighbor_list Named interaction list.
#' @param all_network_genes Character vector.
#'
#' @return A tibble with columns \code{threshold}, \code{tpr}, \code{fpr},
#'   \code{precision}, \code{recall}.
#' @noRd
.compute_aggregate_curves <- function(pathway_list, neighbor_list,
                                      all_network_genes) {

  all_scores <- purrr::map_dfr(pathway_list, function(pathway_df) {
    .compute_scores_vectorized(
      pathway_df$ensembl_gene, neighbor_list, all_network_genes
    ) %>%
      dplyr::mutate(pathway = unique(pathway_df$gs_name))
  })

  thresholds <- sort(unique(c(0, all_scores$ratio, Inf)), decreasing = TRUE)

  purrr::map_dfr(thresholds, function(thr) {
    predicted <- all_scores$ratio >= thr
    tp  <- sum( predicted &  all_scores$in_pathway)
    fp  <- sum( predicted & !all_scores$in_pathway)
    fn  <- sum(!predicted &  all_scores$in_pathway)
    tn  <- sum(!predicted & !all_scores$in_pathway)

    tibble::tibble(
      threshold = thr,
      tpr       = tp / (tp + fn),
      fpr       = fp / (fp + tn),
      precision = if ((tp + fp) > 0) tp / (tp + fp) else 1,
      recall    = tp / (tp + fn)
    )
  })
}
