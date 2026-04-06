# ==============================================================================
# data-raw/network_plots.R
#
# Diagnostic and curve plots derived from a compute_network_cutoff() result.
#
# Depends on: network_helpers.R  (must be sourced first)
#
# EXPORTED FUNCTIONS:
#   plot_cutoff_diagnostics   — seven-panel diagnostic bundle (PNG files)
#   plot_roc_curves           — per-pathway ROC spaghetti + pooled aggregate
#   plot_pr_curves            — per-pathway PR  spaghetti + pooled aggregate
#
# All functions write PNG files to disk and return the ggplot object invisibly.
# ==============================================================================


# ==============================================================================
# plot_cutoff_diagnostics()
# ==============================================================================

#' Generate the full diagnostic plot bundle from a cutoff analysis result
#'
#' Produces seven PNG files saved to \code{plot_dir}:
#' \enumerate{
#'   \item \code{plot_cutoff_distribution.png} — per-pathway best-cutoff histograms
#'   \item \code{plot_roc_curve.png}           — pooled ROC operating-point scatter
#'   \item \code{plot_pr_curve.png}            — pooled PR operating-point scatter
#'   \item \code{plot_auc_density.png}         — per-pathway AUC density curves
#'   \item \code{plot_cutoff_by_pathway_count.png} — cutoff vs. pathway membership
#'   \item \code{plot_multipathway_comparison.png} — strategy comparison
#'   \item \code{plot_sensitivity_specificity.png} — direct sens/spec derivation
#' }
#'
#' Can be called independently on any saved result list without rerunning the
#' full analysis.
#'
#' @param result   List returned by \code{compute_network_cutoff()}.
#' @param plot_dir Character. Directory for PNG output. Created if absent.
#'   Default \code{"data-raw/plots"}.
#'
#' @return Invisibly returns a named list of the seven ggplot objects.
plot_cutoff_diagnostics <- function(result,
                                    plot_dir = "data-raw/plots") {

  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  pathway_results    <- result$pathway_results
  cutoff_roc_summary <- result$cutoff_roc_summary
  cutoff_f1_summary  <- result$cutoff_f1_summary
  kegg_genesets      <- result$kegg_genesets
  validation         <- result$validation
  sens_spec_cutoffs  <- result$sens_spec_cutoffs
  neighbor_list      <- .build_neighbor_list(result$gene_pairs)
  all_network_genes  <- names(neighbor_list)

  cat("  Generating diagnostic plots...\n")

  # ── Plot 1: Per-pathway best-cutoff histograms ──────────────────────────────
  roc_vlines <- tibble(
    label  = c("q50", "q75", "q90"),
    cutoff = c(cutoff_roc_summary$q50, cutoff_roc_summary$q75,
               cutoff_roc_summary$q90)
  )
  f1_vlines <- tibble(
    label  = c("q50", "q75", "q90"),
    cutoff = c(cutoff_f1_summary$q50, cutoff_f1_summary$q75,
               cutoff_f1_summary$q90)
  )

  p1a <- pathway_results %>%
    filter(!is.na(cutoff_roc), auc_roc > 0.5) %>%
    mutate(w = auc_roc - 0.5) %>%
    ggplot(aes(x = cutoff_roc, weight = w)) +
    geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
    geom_vline(data = roc_vlines,
               aes(xintercept = cutoff, colour = label),
               linewidth = 0.9, linetype = "dashed") +
    scale_colour_brewer(palette = "Set1", name = "Quantile") +
    labs(title    = "ROC-based cutoffs",
         subtitle = "AUC-weighted histogram (AUC_ROC > 0.5 pathways only)",
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
    geom_vline(data = f1_vlines,
               aes(xintercept = cutoff, colour = label),
               linewidth = 0.9, linetype = "dashed") +
    scale_colour_brewer(palette = "Set1", name = "Quantile") +
    labs(title    = "PR/F1-based cutoffs",
         subtitle = "AUC-weighted histogram (AUC_PR > prevalence pathways only)",
         x = "Best cutoff (max F1)", y = "Weighted count") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(colour = "grey40"))

  p1 <- cowplot::plot_grid(p1a, p1b, ncol = 2)

  # ── Plot 2: Pooled ROC operating-point scatter ──────────────────────────────
  roc_ops <- pathway_results %>%
    filter(!is.na(cutoff_roc), auc_roc > 0.5) %>%
    mutate(fpr = 1 - specificity, tpr = sensitivity)

  pooled_roc_pts <- roc_ops %>%
    arrange(fpr) %>%
    bind_rows(tibble(fpr = c(0, 1), tpr = c(0, 1),
                     auc_roc = NA_real_, pathway = NA_character_))

  op_points_roc <- validation %>%
    mutate(fpr = 1 - specificity, tpr = sensitivity,
           key = paste0(method, "/", quantile))

  p2 <- ggplot(pooled_roc_pts, aes(x = fpr, y = tpr)) +
    geom_point(data = roc_ops, aes(colour = auc_roc),
               size = 1.8, alpha = 0.55) +
    scale_colour_gradient(low = "#AED6F1", high = "#1A5276",
                          name = "Per-pathway\nAUC-ROC",
                          limits = c(0.5, 1)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dotted", colour = "grey60") +
    geom_point(data = op_points_roc,
               aes(x = fpr, y = tpr, shape = method, fill = quantile),
               size = 4, colour = "black", stroke = 0.6) +
    scale_shape_manual(values = c("ROC" = 21, "PR/F1" = 24),
                       name = "Method") +
    scale_fill_manual(values = c(q50 = "#FEE08B", q75 = "#FC8D59",
                                 q90 = "#D73027"),
                      name = "Quantile") +
    ggrepel::geom_label_repel(
      data = op_points_roc,
      aes(x = fpr, y = tpr, label = label),
      size = 2.8, box.padding = 0.35, max.overlaps = 20, fill = "white"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title    = "Pooled operating points — ROC space",
         subtitle = "Each dot = one KEGG pathway at its best Youden's J cutoff",
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(colour = "grey40"))

  # ── Plot 3: Pooled PR operating-point scatter ───────────────────────────────
  pr_ops2 <- pathway_results %>%
    filter(!is.na(f1_best), sensitivity > 0) %>%
    mutate(
      recall_pt    = sensitivity,
      precision_pt = ifelse((2 * sensitivity - f1_best) > 0,
                            f1_best * sensitivity / (2 * sensitivity - f1_best),
                            NA_real_),
      prevalence   = n_pos / (n_pos + n_neg)
    ) %>%
    filter(!is.na(precision_pt))

  op_points_pr <- validation %>%
    mutate(recall_pt    = sensitivity,
           precision_pt = precision,
           key          = paste0(method, "/", quantile))

  p3 <- ggplot(pr_ops2, aes(x = recall_pt, y = precision_pt)) +
    geom_point(aes(colour = auc_pr), size = 1.8, alpha = 0.55) +
    scale_colour_gradient(low = "#F9E4B7", high = "#B7410E",
                          name = "Per-pathway\nAUC-PR",
                          limits = c(0, 1)) +
    geom_point(data = op_points_pr,
               aes(x = recall_pt, y = precision_pt,
                   shape = method, fill = quantile),
               size = 4, colour = "black", stroke = 0.6) +
    scale_shape_manual(values = c("ROC" = 21, "PR/F1" = 24), name = "Method") +
    scale_fill_manual(values = c(q50 = "#FEE08B", q75 = "#FC8D59",
                                 q90 = "#D73027"),
                      name = "Quantile") +
    ggrepel::geom_label_repel(
      data = op_points_pr,
      aes(x = recall_pt, y = precision_pt, label = label),
      size = 2.8, box.padding = 0.35, max.overlaps = 20, fill = "white"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title    = "Pooled operating points — Precision-Recall space",
         subtitle = "Each dot = one KEGG pathway at its best F1 cutoff",
         x = "Recall (Sensitivity)", y = "Precision") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(colour = "grey40"))

  # ── Plot 4: Per-pathway AUC density ────────────────────────────────────────
  auc_long <- pathway_results %>%
    filter(!is.na(auc_roc) | !is.na(auc_pr)) %>%
    select(pathway, auc_roc, auc_pr) %>%
    pivot_longer(c(auc_roc, auc_pr), names_to = "metric", values_to = "auc") %>%
    filter(!is.na(auc)) %>%
    mutate(metric = recode(metric, auc_roc = "AUC-ROC", auc_pr = "AUC-PR"))

  p4 <- ggplot(auc_long, aes(x = auc, fill = metric, colour = metric)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
    scale_fill_manual(values = c("AUC-ROC" = "steelblue",
                                 "AUC-PR"  = "darkorange"),
                      name = "Metric") +
    scale_colour_manual(values = c("AUC-ROC" = "steelblue",
                                   "AUC-PR"  = "darkorange"),
                        name = "Metric") +
    annotate("text", x = 0.5, y = Inf, label = "Random\n(0.5)",
             hjust = 1.1, vjust = 1.5, size = 3.2, colour = "grey50") +
    labs(title    = "Per-pathway AUC distributions",
         subtitle = "All KEGG pathways with at least one positive and one negative gene",
         x = "AUC", y = "Density") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(colour = "grey40"))

  auc_means <- auc_long %>%
    group_by(metric) %>%
    summarise(mean_auc = mean(auc, na.rm = TRUE), .groups = "drop")

  p4 <- p4 +
    geom_vline(data = auc_means,
               aes(xintercept = mean_auc, colour = metric),
               linetype = "solid", linewidth = 1) +
    geom_text(data = auc_means,
              aes(x = mean_auc, y = Inf,
                  label = sprintf("mean=%.2f", mean_auc),
                  colour = metric),
              hjust = -0.1, vjust = 1.5, size = 3, show.legend = FALSE)

  # ── Plot 5: Cutoff vs number of pathway memberships ─────────────────────────
  gene_cutoff_summary <- kegg_genesets %>%
    distinct(ensembl_gene, gs_name) %>%
    inner_join(pathway_results %>%
                 select(gs_name = pathway, cutoff_roc, auc_roc),
               by = "gs_name", relationship = "many-to-many") %>%
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
                               labels = c("1", "2", "3-5", "6-10",
                                          "11-20", "20+"),
                               right  = TRUE))

  group_medians <- gene_cutoff_summary %>%
    group_by(pathway_group) %>%
    summarise(median_cutoff = median(auc_weighted_cutoff, na.rm = TRUE),
              n_genes = n(), .groups = "drop")

  p5 <- ggplot(gene_cutoff_summary,
               aes(x = pathway_group, y = auc_weighted_cutoff)) +
    geom_boxplot(aes(fill = pathway_group), outlier.shape = NA,
                 alpha = 0.6, width = 0.55) +
    geom_jitter(aes(colour = pathway_group), width = 0.2,
                size = 0.8, alpha = 0.4) +
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
         x       = "Number of KEGG pathways gene belongs to",
         y       = "AUC-weighted best cutoff (ROC)",
         caption = "n = number of genes in each group") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          plot.subtitle = element_text(colour = "grey40"),
          plot.caption  = element_text(colour = "grey40"))

  # ── Plot 6: Multi-pathway gene strategy comparison ───────────────────────────
  mp_comp <- result$multipathway_comparison

  strategy_labels <- c(
    all    = "All genes\n(baseline)",
    remove = "Remove\nmulti-pathway genes",
    random = "Random\npathway assignment",
    best   = "Best-interaction\npathway assignment"
  )

  mp_summary <- imap_dfr(mp_comp, function(strat_res, strat_name) {
    roc_val <- strat_res$validation %>%
      filter(method %in% c("ROC", "ROC/q75") |
               (quantile == "q75" & grepl("ROC", method))) %>%
      filter(row_number() == 1)

    tibble(
      strategy         = strat_name,
      strategy_label   = strategy_labels[strat_name],
      n_genes_used     = strat_res$n_genes_used,
      n_pathways_used  = strat_res$n_pathways_used,
      cutoff_roc_q75   = strat_res$cutoff_roc_summary$q75,
      cutoff_f1_q75    = strat_res$cutoff_f1_summary$q75,
      mean_auc_roc     = strat_res$cutoff_roc_summary$mean_auc_roc,
      mean_sensitivity = strat_res$cutoff_roc_summary$mean_sensitivity,
      mean_specificity = strat_res$cutoff_roc_summary$mean_specificity,
      sensitivity_val  = if (nrow(roc_val) > 0) roc_val$sensitivity[[1]] else NA_real_,
      specificity_val  = if (nrow(roc_val) > 0) roc_val$specificity[[1]] else NA_real_,
      precision_val    = if (nrow(roc_val) > 0) roc_val$precision[[1]]   else NA_real_,
      f1_val           = if (nrow(roc_val) > 0) roc_val$f1[[1]]          else NA_real_
    )
  }) %>%
    mutate(strategy_label = factor(
      strategy_label,
      levels = strategy_labels[c("all", "remove", "random", "best")]
    ))

  mp_long <- mp_summary %>%
    select(strategy_label, sensitivity_val, specificity_val,
           precision_val, f1_val) %>%
    pivot_longer(-strategy_label, names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric,
                           sensitivity_val = "Sensitivity",
                           specificity_val = "Specificity",
                           precision_val   = "Precision",
                           f1_val          = "F1"))

  p6a <- ggplot(mp_long, aes(x = strategy_label, y = value, fill = metric)) +
    geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.3f", value)),
              position = position_dodge(0.75), vjust = -0.4, size = 2.8) +
    scale_fill_brewer(palette = "Set2", name = "Metric") +
    scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
    labs(title    = "Validation metrics at ROC q75 cutoff",
         subtitle = "By multi-pathway gene handling strategy",
         x = NULL, y = "Score") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(size = 9),
          plot.subtitle = element_text(colour = "grey40"))

  mp_cutoff <- mp_summary %>%
    select(strategy_label, cutoff_roc_q75, cutoff_f1_q75) %>%
    pivot_longer(-strategy_label, names_to = "method", values_to = "cutoff") %>%
    mutate(method = recode(method,
                           cutoff_roc_q75 = "ROC q75",
                           cutoff_f1_q75  = "PR/F1 q75"))

  p6b <- ggplot(mp_cutoff,
                aes(x = strategy_label, y = cutoff,
                    colour = method, group = method)) +
    geom_line(linewidth = 0.9, linetype = "dashed") +
    geom_point(size = 4) +
    geom_text(aes(label = sprintf("%.3f", cutoff)),
              vjust = -0.9, size = 3, show.legend = FALSE) +
    scale_colour_manual(values = c("ROC q75" = "#1A5EA8",
                                   "PR/F1 q75" = "#E76F51"),
                        name = "Cutoff\nmethod") +
    labs(title    = "Recommended cutoff across strategies",
         subtitle = "q75 thresholds (ROC and PR/F1)",
         x = NULL, y = "Cutoff value") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(size = 9),
          plot.subtitle = element_text(colour = "grey40"))

  p6 <- cowplot::plot_grid(p6a, p6b, ncol = 1, rel_heights = c(1.3, 1))

  # ── Plot 7: Sensitivity and specificity vs threshold ─────────────────────────
  sens_spec_ops <- pathway_results %>%
    filter(!is.na(sensitivity), !is.na(specificity), !is.na(cutoff_roc)) %>%
    select(pathway, cutoff_roc, sensitivity, specificity, auc_roc)

  direct_vlines <- sens_spec_cutoffs %>%
    mutate(linetype_lab = method, colour_lab = method)

  p7a <- ggplot(sens_spec_ops, aes(x = cutoff_roc)) +
    geom_point(aes(y = sensitivity, colour = auc_roc),
               size = 1.5, alpha = 0.5, shape = 16) +
    geom_point(aes(y = specificity, colour = auc_roc),
               size = 1.5, alpha = 0.5, shape = 17) +
    scale_colour_gradient(low = "#AED6F1", high = "#1A5276",
                          name = "Per-pathway\nAUC-ROC") +
    geom_vline(data = direct_vlines,
               aes(xintercept = cutoff, linetype = method),
               colour = "firebrick", linewidth = 0.9) +
    scale_linetype_manual(
      values = c("Pooled Youden's J"        = "solid",
                 "Sensitivity = Specificity" = "dashed",
                 "Min distance to (0,1)"     = "dotted"),
      name = "Direct cutoff"
    ) +
    annotate("text", x = -Inf, y = 0.98, label = "\u25b2 Specificity",
             hjust = -0.1, size = 3.2, colour = "grey30") +
    annotate("text", x = -Inf, y = 0.93, label = "\u25cf Sensitivity",
             hjust = -0.1, size = 3.2, colour = "grey30") +
    labs(title    = "Per-pathway sensitivity and specificity at best Youden's J",
         subtitle = paste("Circles = sensitivity, triangles = specificity;",
                          "vertical lines = direct pooled cutoffs"),
         x = "Cutoff (ratio threshold at best Youden's J per pathway)",
         y = "Metric value") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(colour = "grey40"))

  direct_ops <- sens_spec_cutoffs %>%
    mutate(fpr = 1 - specificity)

  p7b <- ggplot(sens_spec_ops,
                aes(x = 1 - specificity, y = sensitivity,
                    colour = auc_roc)) +
    geom_point(size = 1.8, alpha = 0.55) +
    scale_colour_gradient(low = "#AED6F1", high = "#1A5276",
                          name = "Per-pathway\nAUC-ROC",
                          limits = c(0.5, 1)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dotted", colour = "grey60") +
    geom_point(data = direct_ops,
               aes(x = fpr, y = sensitivity),
               colour = "firebrick", size = 5, shape = 18,
               inherit.aes = FALSE) +
    ggrepel::geom_label_repel(
      data = direct_ops,
      aes(x = fpr, y = sensitivity,
          label = sprintf("%s\nc=%.3f", method, cutoff)),
      colour = "firebrick", size = 2.8, box.padding = 0.4,
      max.overlaps = 10, fill = "white", inherit.aes = FALSE
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title    = "Sensitivity vs specificity trade-off",
         subtitle = paste("Each dot = one pathway at its best Youden's J;",
                          "diamonds = direct pooled cutoffs"),
         x = "1 \u2212 Specificity (FPR)",
         y = "Sensitivity (TPR)") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle = element_text(colour = "grey40"))

  p7 <- cowplot::plot_grid(p7a, p7b, ncol = 2, rel_widths = c(1.1, 1))

  # ── Save ─────────────────────────────────────────────────────────────────────
  plots <- list(
    plot_cutoff_distribution     = list(p = p1, w = 14, h = 5),
    plot_roc_curve               = list(p = p2, w =  9, h = 7),
    plot_pr_curve                = list(p = p3, w =  9, h = 7),
    plot_auc_density             = list(p = p4, w = 10, h = 5),
    plot_cutoff_by_pathway_count = list(p = p5, w = 10, h = 6),
    plot_multipathway_comparison = list(p = p6, w = 12, h = 9),
    plot_sensitivity_specificity = list(p = p7, w = 16, h = 6)
  )

  for (nm in names(plots)) {
    out <- file.path(plot_dir, paste0(nm, ".png"))
    ggsave(out, plots[[nm]]$p,
           width = plots[[nm]]$w, height = plots[[nm]]$h, dpi = 150)
    cat(sprintf("    Saved: %s\n", out))
  }

  invisible(plots)
}


# ==============================================================================
# plot_roc_curves()
# ==============================================================================

#' Plot all per-pathway ROC curves with aggregate pooled curve on top
#'
#' Each grey line is one KEGG pathway's ROC curve. Colour saturation encodes
#' per-pathway AUC-ROC (darker = higher AUC). Three coloured operating points
#' mark the q50/q75/q90 quantile cutoffs on the pooled aggregate curve.
#'
#' @param result          List returned by \code{compute_network_cutoff()}.
#' @param n_sample        Integer. Number of pathways to draw as background
#'   lines (random sample for readability). Default \code{300}. Pass
#'   \code{Inf} to draw all valid pathways.
#' @param alpha_pathway   Numeric in \[0, 1\]. Opacity of background lines.
#'   Default \code{0.08}.
#' @param plot_dir        Character. Output directory.
#' @param filename        Character. Output filename.
#'
#' @return A ggplot object, invisibly.
plot_roc_curves <- function(result,
                            n_sample      = 300,
                            alpha_pathway = 0.08,
                            plot_dir      = "data-raw/plots",
                            filename      = "plot_roc_all_curves.png") {

  cat("  Building full per-pathway ROC curves...\n")

  pathway_results    <- result$pathway_results
  cutoff_roc_summary <- result$cutoff_roc_summary
  kegg_genesets      <- result$kegg_genesets
  neighbor_list      <- .build_neighbor_list(result$gene_pairs)
  all_network_genes  <- names(neighbor_list)
  pathway_list       <- kegg_genesets %>%
    dplyr::group_by(gs_name) %>%
    dplyr::group_split()

  # ── Per-pathway full curves (sampled for readability) ─────────────────────
  valid_pathways <- pathway_results %>%
    dplyr::filter(!is.na(auc_roc)) %>%
    dplyr::pull(pathway)

  sample_n <- min(n_sample, length(valid_pathways))
  set.seed(42L)
  sampled   <- sample(valid_pathways, sample_n)

  sampled_list <- Filter(
    function(df) unique(df$gs_name) %in% sampled,
    pathway_list
  )

  # .full_curves_for_pathway() already attaches auc_roc and auc_pr to every
  # row, so no join is needed here. Joining would create auc_roc.x / auc_roc.y
  # and break the aes(colour = auc_roc) mapping downstream.
  per_pathway_curves <- purrr::map_dfr(sampled_list, function(pathway_df) {
    .full_curves_for_pathway(
      pathway_name      = unique(pathway_df$gs_name),
      pathway_genes     = pathway_df$ensembl_gene,
      neighbor_list     = neighbor_list,
      all_network_genes = all_network_genes
    )
  })

  cat(sprintf("    %d pathways drawn as background\n", sample_n))

  # ── Aggregate ROC curve (all pathways pooled) ──────────────────────────────
  cat("  Computing pooled aggregate ROC curve...\n")
  agg_curve <- .compute_aggregate_curves(
    pathway_list, neighbor_list, all_network_genes
  ) %>%
    dplyr::arrange(fpr, tpr)

  # ── Operating points at q50 / q75 / q90 ───────────────────────────────────
  cutoffs <- c(
    q50 = cutoff_roc_summary$q50,
    q75 = cutoff_roc_summary$q75,
    q90 = cutoff_roc_summary$q90
  )

  op_points <- purrr::imap_dfr(cutoffs, function(thr, label) {
    idx <- which.min(abs(agg_curve$threshold - thr))
    agg_curve[idx, ] %>%
      dplyr::mutate(label = label, cutoff = thr)
  })

  # ── Build plot ─────────────────────────────────────────────────────────────
  # KEY FIX: The aggregate curve is drawn with a fixed colour argument outside
  # aes() so it does NOT inherit from the per-pathway colour scale. The
  # op_points use a separate fill scale. This avoids the "object 'auc_roc' not
  # found" error that occurs when geom_path/geom_point layers try to evaluate
  # a colour aesthetic that is present in the global scale but absent in their
  # data frame.
  quantile_colours <- c(q50 = "#E41A1C", q75 = "#377EB8", q90 = "#4DAF4A")

  p <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                colour = "grey60", linewidth = 0.4) +
    # Per-pathway background: colour = auc_roc from per_pathway_curves
    geom_path(
      data    = per_pathway_curves %>%
        dplyr::arrange(pathway, fpr, tpr),
      mapping = aes(x = fpr, y = tpr, group = pathway, colour = auc_roc),
      alpha   = alpha_pathway,
      linewidth = 0.35
    ) +
    scale_colour_gradient(
      low      = "grey80",
      high     = "#08306B",
      name     = "Per-pathway\nAUC-ROC",
      na.value = "transparent",
      guide    = guide_colourbar(barheight = 6, order = 1)
    ) +
    # Aggregate curve: fixed colour — no aes mapping, no conflict with scale
    geom_path(
      data      = agg_curve,
      mapping   = aes(x = fpr, y = tpr),
      colour    = "#333333",
      linewidth = 1.1
    ) +
    # Quantile operating points: separate fill scale via shape 21
    geom_point(
      data        = op_points,
      mapping     = aes(x = fpr, y = tpr, fill = label),
      shape       = 21,
      size        = 4,
      colour      = "white",
      stroke      = 1.2,
      inherit.aes = FALSE
    ) +
    ggrepel::geom_label_repel(
      data        = op_points,
      mapping     = aes(x = fpr, y = tpr,
                        label = sprintf("%s  c=%.3f", label, cutoff),
                        fill  = label),
      colour      = "white",
      size        = 3,
      fontface    = "bold",
      box.padding = 0.4,
      min.segment.length = 0,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = quantile_colours,
      name   = "ROC quantile\ncutoff",
      guide  = guide_legend(
        order = 2,
        override.aes = list(size = 4, shape = 21,
                            colour = "white", stroke = 1.2)
      )
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title    = "Per-pathway ROC curves",
      subtitle = sprintf(
        "%d sampled pathways (grey) · pooled aggregate curve (black) · quantile operating points",
        sample_n
      ),
      x = "False Positive Rate (1 \u2212 Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.subtitle    = element_text(colour = "grey40"),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )

  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(plot_dir, filename)
  ggsave(out, p, width = 8, height = 7, dpi = 150)
  cat(sprintf("    Saved: %s\n", out))

  invisible(p)
}


# ==============================================================================
# plot_pr_curves()
# ==============================================================================

#' Plot all per-pathway PR curves with aggregate pooled curve on top
#'
#' Each grey line is one KEGG pathway's precision-recall curve. Colour
#' saturation encodes per-pathway AUC-PR. Three coloured operating points mark
#' the q50/q75/q90 PR/F1 quantile cutoffs on the pooled aggregate curve. A
#' dashed horizontal line marks the mean pathway prevalence (random baseline).
#'
#' @param result          List returned by \code{compute_network_cutoff()}.
#' @param n_sample        Integer. Number of pathways to draw. Default \code{300}.
#' @param alpha_pathway   Numeric in \[0, 1\]. Default \code{0.08}.
#' @param plot_dir        Character. Output directory.
#' @param filename        Character. Output filename.
#'
#' @return A ggplot object, invisibly.
plot_pr_curves <- function(result,
                           n_sample      = 300,
                           alpha_pathway = 0.08,
                           plot_dir      = "data-raw/plots",
                           filename      = "plot_pr_all_curves.png") {

  cat("  Building full per-pathway PR curves...\n")

  pathway_results   <- result$pathway_results
  cutoff_f1_summary <- result$cutoff_f1_summary
  kegg_genesets     <- result$kegg_genesets
  neighbor_list     <- .build_neighbor_list(result$gene_pairs)
  all_network_genes <- names(neighbor_list)
  pathway_list      <- kegg_genesets %>%
    dplyr::group_by(gs_name) %>%
    dplyr::group_split()

  # ── Per-pathway full curves (sampled) ─────────────────────────────────────
  valid_pathways <- pathway_results %>%
    dplyr::filter(!is.na(auc_pr)) %>%
    dplyr::pull(pathway)

  sample_n <- min(n_sample, length(valid_pathways))
  set.seed(42L)
  sampled   <- sample(valid_pathways, sample_n)

  sampled_list <- Filter(
    function(df) unique(df$gs_name) %in% sampled,
    pathway_list
  )

  # .full_curves_for_pathway() already attaches auc_roc and auc_pr to every
  # row, so no join is needed here.
  per_pathway_curves <- purrr::map_dfr(sampled_list, function(pathway_df) {
    .full_curves_for_pathway(
      pathway_name      = unique(pathway_df$gs_name),
      pathway_genes     = pathway_df$ensembl_gene,
      neighbor_list     = neighbor_list,
      all_network_genes = all_network_genes
    )
  })

  cat(sprintf("    %d pathways drawn as background\n", sample_n))

  # ── Aggregate PR curve ─────────────────────────────────────────────────────
  cat("  Computing pooled aggregate PR curve...\n")
  agg_curve <- .compute_aggregate_curves(
    pathway_list, neighbor_list, all_network_genes
  ) %>%
    dplyr::arrange(recall, dplyr::desc(precision))

  # ── Operating points ───────────────────────────────────────────────────────
  cutoffs <- c(
    q50 = cutoff_f1_summary$q50,
    q75 = cutoff_f1_summary$q75,
    q90 = cutoff_f1_summary$q90
  )

  op_points <- purrr::imap_dfr(cutoffs, function(thr, label) {
    idx <- which.min(abs(agg_curve$threshold - thr))
    agg_curve[idx, ] %>%
      dplyr::mutate(label = label, cutoff = thr)
  })

  # ── Mean prevalence baseline ───────────────────────────────────────────────
  mean_prevalence <- pathway_results %>%
    dplyr::filter(!is.na(n_pos), !is.na(n_neg)) %>%
    dplyr::mutate(prev = n_pos / (n_pos + n_neg)) %>%
    dplyr::summarise(mean_prev = mean(prev)) %>%
    dplyr::pull(mean_prev)

  # ── Build plot ─────────────────────────────────────────────────────────────
  # Same approach as plot_roc_curves(): aggregate curve drawn with fixed colour
  # outside aes() to avoid the 'auc_pr not found' error.
  quantile_colours <- c(q50 = "#E41A1C", q75 = "#377EB8", q90 = "#4DAF4A")

  p <- ggplot() +
    geom_hline(yintercept = mean_prevalence, linetype = "dashed",
               colour = "grey50", linewidth = 0.5) +
    annotate("text",
             x     = 0.98,
             y     = mean_prevalence + 0.012,
             label = sprintf("Random baseline\n(mean prevalence = %.3f)",
                             mean_prevalence),
             hjust = 1, size = 2.8, colour = "grey50") +
    # Per-pathway spaghetti: colour = auc_pr from per_pathway_curves
    geom_path(
      data    = per_pathway_curves %>%
        dplyr::arrange(pathway, recall, dplyr::desc(precision)),
      mapping = aes(x = recall, y = precision,
                    group = pathway, colour = auc_pr),
      alpha   = alpha_pathway,
      linewidth = 0.35
    ) +
    scale_colour_gradient(
      low      = "grey80",
      high     = "#7F0000",
      name     = "Per-pathway\nAUC-PR",
      na.value = "transparent",
      guide    = guide_colourbar(barheight = 6, order = 1)
    ) +
    # Aggregate PR curve: fixed colour, no aes conflict
    geom_path(
      data      = agg_curve,
      mapping   = aes(x = recall, y = precision),
      colour    = "#333333",
      linewidth = 1.1
    ) +
    geom_point(
      data        = op_points,
      mapping     = aes(x = recall, y = precision, fill = label),
      shape       = 21,
      size        = 4,
      colour      = "white",
      stroke      = 1.2,
      inherit.aes = FALSE
    ) +
    ggrepel::geom_label_repel(
      data        = op_points,
      mapping     = aes(x = recall, y = precision,
                        label = sprintf("%s  c=%.3f", label, cutoff),
                        fill  = label),
      colour      = "white",
      size        = 3,
      fontface    = "bold",
      box.padding = 0.4,
      min.segment.length = 0,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = quantile_colours,
      name   = "PR/F1 quantile\ncutoff",
      guide  = guide_legend(
        order = 2,
        override.aes = list(size = 4, shape = 21,
                            colour = "white", stroke = 1.2)
      )
    ) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title    = "Per-pathway Precision-Recall curves",
      subtitle = sprintf(
        "%d sampled pathways (grey) · pooled aggregate curve (black) · quantile operating points",
        sample_n
      ),
      x = "Recall (Sensitivity)",
      y = "Precision"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.subtitle    = element_text(colour = "grey40"),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )

  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(plot_dir, filename)
  ggsave(out, p, width = 8, height = 7, dpi = 150)
  cat(sprintf("    Saved: %s\n", out))

  invisible(p)
}
