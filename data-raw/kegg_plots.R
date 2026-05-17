# ==============================================================================
# data-raw/kegg_plots.R
#
# KEGG pathway characterisation and gene membership distribution plots.
# These are exploratory / thesis figures — they do NOT modify package internals.
#
# Depends on: network_helpers.R  (must be sourced first)
#
# EXPORTED FUNCTIONS:
#   plot_kegg_characterisation    — five-panel KEGG characterisation bundle
#   plot_gene_pathway_membership  — three-panel gene membership distribution
#
# USAGE:
#   source("data-raw/network_helpers.R")
#   source("data-raw/kegg_plots.R")
#
#   result <- <your saved result list>
#   plot_kegg_characterisation(result)
#   plot_gene_pathway_membership(result)
#
# EXTRA DEPENDENCIES (beyond network_helpers.R):
#   install.packages(c("ggrepel", "forcats", "stringr",
#                      "scales", "RColorBrewer"))
# ==============================================================================


# ==============================================================================
# INTERNAL HELPERS
# ==============================================================================

#' Infer a human-readable biological category from a KEGG pathway name
#'
#' KEGG Legacy names look like:  KEGG_RIBOSOME
#' KEGG Medicus names look like: KEGG_MEDICUS_REFERENCE_DIABETES_MELLITUS
#'
#' @param gs_name Character vector of pathway names.
#' @return Character vector of category labels.
#' @noRd
.infer_category <- function(gs_name) {
  nm <- tolower(gs_name)

  dplyr::case_when(
    # Metabolism
    stringr::str_detect(nm,
                        "metabol|biosynthes|degradation|oxidative_phosphorylation|
       fatty_acid|glycol|tca|citrate|purine|pyrimidine|
       amino_acid|carbon|drug_metabol") ~ "Metabolism",

    # Signalling
    stringr::str_detect(nm,
                        "signaling|signal_transduction|jak|mapk|pi3k|mtor|
       wnt|notch|hedgehog|tgf|vegf|egfr|fgf|igf|
       calcium|camp|cgmp|hippo|nf.kb|toll|jak.stat") ~ "Signalling",

    # Immune / Infection
    stringr::str_detect(nm,
                        "immun|inflam|cytokine|chemokine|interferon|
       t_cell|b_cell|nk_cell|complement|autoimmun|
       infection|virus|bacteria|parasite|hiv|hepatitis|
       tuberculosis|malaria|influenza|covid|herpes|
       leishmania|leishmaniasis|salmonella|staphylo|
       legionella|measles|papillomavirus|prion") ~ "Immune / Infection",

    # Cancer
    stringr::str_detect(nm,
                        "cancer|carcinoma|glioma|leukemia|melanoma|
       colorectal|pancreatic|prostate|breast|bladder|
       renal|thyroid|endometrial|small_cell_lung|
       non_small_cell|oncogene|tumor|tumour") ~ "Cancer",

    # Cell Cycle / DNA
    stringr::str_detect(nm,
                        "cell_cycle|dna_repair|dna_replication|
       mismatch_repair|nucleotide_excision|
       base_excision|homologous_recomb|p53|
       apoptosis|autophagy|senescence") ~ "Cell Cycle / DNA",

    # Nervous System
    stringr::str_detect(nm,
                        "neurodegener|alzheimer|parkinson|huntington|
       amyotrophic|prion|synapse|long.term|
       serotonin|dopamine|glutamate|gaba|
       neurotrophin|axon|circadian|addiction") ~ "Nervous System",

    # Endocrine
    stringr::str_detect(nm,
                        "diabetes|insulin|adipocyte|thyroid_hormone|
       cortisol|estrogen|progesterone|androgen|
       oxytocin|prolactin|gnrh|adipogenesis") ~ "Endocrine",

    # Development / ECM
    stringr::str_detect(nm,
                        "stem_cell|pluripoten|differentia|embryo|
       axon_guidance|focal_adhesion|ecm|
       extracellular_matrix|gap_junction|
       tight_junction|adherens") ~ "Development / ECM",

    # RNA / Protein Processing
    stringr::str_detect(nm,
                        "ribosom|spliceosom|proteasom|ubiquitin|
       rna_transport|mrna|protein_export|
       phagosom|lysosom|endocytosis|
       peroxisom|mitophagy") ~ "RNA / Protein Processing",

    TRUE ~ "Other"
  )
}


#' Build the master per-pathway summary table used by all KEGG plots
#'
#' @param result List returned by \code{compute_network_cutoff()}.
#' @return A tibble with one row per pathway and columns: \code{pathway},
#'   \code{n_genes_kegg}, \code{n_genes_network}, \code{coverage},
#'   all columns from \code{result$pathway_results}, \code{category},
#'   \code{short_name}.
#' @noRd
.build_pathway_summary <- function(result) {
  kegg      <- result$kegg_genesets
  pr        <- result$pathway_results
  net_genes <- names(.build_neighbor_list(result$gene_pairs))

  pathway_genes <- kegg %>%
    group_by(pathway = gs_name) %>%
    summarise(
      n_genes_kegg    = n_distinct(ensembl_gene),
      n_genes_network = n_distinct(ensembl_gene[ensembl_gene %in% net_genes]),
      coverage        = n_genes_network / n_genes_kegg,
      .groups         = "drop"
    )

  pathway_genes %>%
    left_join(pr, by = "pathway") %>%
    mutate(
      category   = .infer_category(pathway),
      short_name = pathway %>%
        stringr::str_remove("^KEGG_MEDICUS_REFERENCE_") %>%
        stringr::str_remove("^KEGG_LEGACY_") %>%
        stringr::str_remove("^KEGG_") %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_to_title()
    )
}


# ==============================================================================
# PLOT 1 — Pathway size distribution by category
# ==============================================================================

.plot_pathway_sizes <- function(psum, plot_dir, cat_palette) {
  cat("  [1/5] Pathway size by category...\n")

  cat_order <- psum %>%
    count(category, sort = TRUE) %>%
    pull(category)

  psum_ord <- psum %>%
    mutate(category = factor(category, levels = cat_order))

  cat_ann <- psum_ord %>%
    group_by(category) %>%
    summarise(
      n         = n(),
      median_sz = median(n_genes_kegg),
      .groups   = "drop"
    ) %>%
    mutate(label = sprintf("n=%d\nmedian=%d", n, median_sz))

  p <- ggplot(psum_ord,
              aes(x = category, y = n_genes_kegg, colour = category)) +
    geom_jitter(width = 0.25, size = 1.2, alpha = 0.55) +
    geom_boxplot(width = 0.35, outlier.shape = NA,
                 fill = NA, colour = "grey30", linewidth = 0.5) +
    geom_text(data = cat_ann,
              mapping = aes(x = category, y = -Inf, label = label),
              inherit.aes = FALSE,
              vjust = -0.3, size = 2.7, colour = "grey40") +
    scale_y_log10(breaks = c(5, 10, 25, 50, 100, 200, 500),
                  labels = scales::comma) +
    scale_colour_manual(values = cat_palette, guide = "none") +
    coord_flip() +
    labs(title    = "Pathway size distribution by category",
         subtitle = "Each dot = one KEGG pathway  ·  y-axis log-scaled",
         x = NULL, y = "Number of genes in pathway",
         caption = "n = pathways in category, median = median gene count") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle      = element_text(colour = "grey40"),
          plot.caption       = element_text(colour = "grey40"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank())

  out <- file.path(plot_dir, "01_pathway_size_by_category.png")
  ggsave(out, p, width = 10, height = 7, dpi = 150)
  cat(sprintf("    Saved: %s\n", out))
  invisible(p)
}


# ==============================================================================
# PLOT 2 — Network coverage per pathway, ordered within category
# ==============================================================================

.plot_network_coverage <- function(psum, plot_dir, cat_palette) {
  cat("  [2/5] Network coverage by category...\n")

  psum_v <- psum %>% filter(!is.na(coverage))

  cat_order <- psum_v %>%
    group_by(category) %>%
    summarise(med = median(coverage), .groups = "drop") %>%
    arrange(desc(med)) %>%
    pull(category)

  psum_v <- psum_v %>%
    mutate(category = factor(category, levels = cat_order))

  cat_summ <- psum_v %>%
    group_by(category) %>%
    summarise(
      med_cov  = median(coverage),
      pct_half = round(100 * mean(coverage >= 0.5)),
      n        = n(),
      .groups  = "drop"
    ) %>%
    mutate(category = factor(category, levels = cat_order))

  p <- ggplot(psum_v, aes(x = category, y = coverage, colour = category)) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "firebrick", linewidth = 0.5) +
    geom_jitter(width = 0.22, size = 0.9, alpha = 0.5) +
    geom_boxplot(width = 0.35, outlier.shape = NA,
                 fill = NA, colour = "grey30", linewidth = 0.5) +
    geom_text(data = cat_summ,
              mapping = aes(x = category, y = -0.02,
                            label = sprintf("\u226550%%: %d%%", pct_half)),
              inherit.aes = FALSE,
              hjust = 1, size = 2.7, colour = "grey40") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(-0.05, 1.02)) +
    scale_colour_manual(values = cat_palette, guide = "none") +
    coord_flip() +
    labs(title    = "Network coverage per pathway by category",
         subtitle = "Fraction of pathway genes present in the interaction network",
         x = NULL,
         y = "Coverage (genes in network / genes in pathway)",
         caption = paste("Red dashed line = 50% coverage  ·",
                         "'\u226550%' = % of pathways in category above threshold")) +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle      = element_text(colour = "grey40"),
          plot.caption       = element_text(colour = "grey40"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank())

  out <- file.path(plot_dir, "02_network_coverage_by_category.png")
  ggsave(out, p, width = 10, height = 7, dpi = 150)
  cat(sprintf("    Saved: %s\n", out))
  invisible(p)
}


# ==============================================================================
# PLOT 3 — Pathway size vs network coverage scatter (plain + enriched)
# ==============================================================================

.plot_size_vs_coverage <- function(psum, plot_dir, cat_palette) {
  cat("  [3/6] Pathway size vs network coverage scatter...\n")

  psum_v <- psum %>% filter(!is.na(coverage), !is.na(n_genes_kegg))

  # Shared axes used by both versions
  x_scale <- scale_x_log10(
    breaks = c(5, 10, 25, 50, 100, 200, 500),
    labels = scales::comma
  )
  y_scale <- scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.02)
  )

  # ── 3a: Plain version ────────────────────────────────────────────────────────
  p_plain <- ggplot(psum_v, aes(x = n_genes_kegg, y = coverage,
                                colour = category)) +
    geom_point(size = 1.8, alpha = 0.55) +
    x_scale +
    y_scale +
    scale_colour_manual(values = cat_palette, name = "Category") +
    labs(
      title = "Pathway size vs network coverage",
      x     = "Pathway size (number of genes, log scale)",
      y     = "Network coverage (fraction of pathway genes in network)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "right")

  out_plain <- file.path(plot_dir, "03a_size_vs_coverage_plain.png")
  ggsave(out_plain, p_plain, width = 8, height = 6, dpi = 150)
  cat(sprintf("    Saved: %s\n", out_plain))

  # ── 3b: Enriched version ─────────────────────────────────────────────────────
  # Quadrant reference lines: minimum size = 15, minimum coverage = 50%
  size_threshold     <- 15
  coverage_threshold <- 0.5

  n_passing <- psum_v %>%
    filter(n_genes_kegg >= size_threshold, coverage >= coverage_threshold) %>%
    nrow()
  n_total <- nrow(psum_v)

  p_enriched <- ggplot(psum_v,
                       aes(x = n_genes_kegg, y = coverage,
                           colour = category, size = auc_pr)) +
    # Quadrant shading — top-right = ideal pathways
    annotate("rect",
             xmin = size_threshold, xmax = Inf,
             ymin = coverage_threshold, ymax = Inf,
             fill = "#E8F5E9", alpha = 0.4) +
    # Quadrant label
    annotate("text",
             x = max(psum_v$n_genes_kegg) * 0.95,
             y = coverage_threshold + 0.02,
             label = sprintf("Candidates\nn = %d / %d (%.0f%%)",
                             n_passing, n_total,
                             100 * n_passing / n_total),
             hjust = 1, vjust = 0, size = 3.2,
             colour = "#388E3C", fontface = "italic") +
    # Points
    geom_point(alpha = 0.6) +
    x_scale +
    y_scale +
    scale_colour_manual(values = cat_palette, name = "Category") +
    scale_size_continuous(
      name   = "AUC-PR",
      range  = c(0.8, 5),
      breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    ) +
    labs(
      title    = "Pathway size vs network coverage",
      subtitle = "Each dot = one KEGG pathway  \u00b7  x-axis log-scaled  \u00b7  point size = AUC-PR",
      x       = "Pathway size (number of genes, log scale)",
      y       = "Network coverage (fraction of pathway genes in network)",
      caption = sprintf(
        "Green shading = candidate pathways (\u2265%d genes & \u226550%% coverage)",
        size_threshold
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.subtitle    = element_text(colour = "grey40"),
      plot.caption     = element_text(colour = "grey40"),
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )

  out_enriched <- file.path(plot_dir, "03b_size_vs_coverage_enriched.png")
  ggsave(out_enriched, p_enriched, width = 12, height = 8, dpi = 150)
  cat(sprintf("    Saved: %s\n", out_enriched))

  invisible(list(plain = p_plain, enriched = p_enriched))
}


# ==============================================================================
# PLOT 4 — Coverage vs AUC-ROC scatter
# ==============================================================================

.plot_coverage_vs_auc <- function(psum, plot_dir, cat_palette) {
  cat("  [4/6] Coverage vs AUC scatter...\n")

  psum_v <- psum %>% filter(!is.na(coverage), !is.na(auc_roc))

  label_df <- bind_rows(
    psum_v %>% slice_max(auc_roc, n = 6, with_ties = FALSE),
    psum_v %>% slice_min(auc_roc, n = 6, with_ties = FALSE),
    psum_v %>%
      filter(coverage < 0.2, auc_roc > 0.8) %>%
      slice_max(auc_roc, n = 4, with_ties = FALSE)
  ) %>%
    distinct(pathway, .keep_all = TRUE)

  p <- ggplot(psum_v,
              aes(x = coverage, y = auc_roc,
                  colour = category, size = n_genes_kegg)) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_point(alpha = 0.55) +
    geom_smooth(data = psum_v,
                mapping = aes(x = coverage, y = auc_roc),
                method = "loess", formula = y ~ x,
                colour = "grey20", fill = "grey85",
                linewidth = 0.9, se = TRUE, alpha = 0.3,
                inherit.aes = FALSE) +
    ggrepel::geom_label_repel(
      data = label_df,
      mapping = aes(x = coverage, y = auc_roc, label = short_name),
      size = 2.5, box.padding = 0.35, max.overlaps = 20,
      segment.colour = "grey50", inherit.aes = FALSE
    ) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_size_continuous(name = "Pathway size\n(# genes)",
                          range = c(1, 6),
                          breaks = c(10, 50, 100, 200)) +
    scale_colour_manual(values = cat_palette, name = "Category") +
    labs(title    = "Network coverage vs discriminative power (AUC-ROC)",
         subtitle = "Each dot = one KEGG pathway  ·  LOESS trend across all pathways",
         x       = "Network coverage (fraction of pathway genes in network)",
         y       = "AUC-ROC (best Youden's J per pathway)",
         caption = "Dashed line = random classifier (AUC = 0.5)") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle    = element_text(colour = "grey40"),
          plot.caption     = element_text(colour = "grey40"),
          panel.grid.minor = element_blank(),
          legend.position  = "right")

  out <- file.path(plot_dir, "04_coverage_vs_auc.png")
  ggsave(out, p, width = 11, height = 8, dpi = 150)
  cat(sprintf("    Saved: %s\n", out))
  invisible(p)
}


# ==============================================================================
# PLOT 4 — AUC-ROC and AUC-PR by category
# ==============================================================================

.plot_auc_by_category <- function(psum, plot_dir, cat_palette) {
  cat("  [5/6] AUC by category...\n")

  psum_v <- psum %>% filter(!is.na(auc_roc))

  cat_order <- psum_v %>%
    group_by(category) %>%
    summarise(med = median(auc_roc, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med)) %>%
    pull(category)

  plong <- psum_v %>%
    select(pathway, category, n_genes_kegg, auc_roc, auc_pr) %>%
    pivot_longer(cols = c(auc_roc, auc_pr),
                 names_to  = "metric",
                 values_to = "auc") %>%
    filter(!is.na(auc)) %>%
    mutate(
      category = factor(category, levels = cat_order),
      metric   = factor(metric,
                        levels = c("auc_roc", "auc_pr"),
                        labels = c("AUC-ROC", "AUC-PR"))
    )

  cat_ann <- plong %>%
    group_by(category, metric) %>%
    summarise(n = n(), med = median(auc, na.rm = TRUE), .groups = "drop") %>%
    mutate(category = factor(category, levels = cat_order))

  p <- ggplot(plong,
              aes(x = category, y = auc, fill = category, colour = category)) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "grey55", linewidth = 0.4) +
    geom_violin(alpha = 0.25, linewidth = 0.3, trim = TRUE) +
    geom_boxplot(width = 0.25, outlier.shape = NA,
                 alpha = 0.6, linewidth = 0.45, colour = "grey20") +
    geom_text(data = cat_ann,
              mapping = aes(x = category, y = -0.03,
                            label = sprintf("n=%d\nmed=%.2f", n, med)),
              inherit.aes = FALSE,
              size = 2.4, colour = "grey40") +
    facet_wrap(~ metric, ncol = 2) +
    scale_y_continuous(limits = c(-0.06, 1.02)) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_colour_manual(values = cat_palette, guide = "none") +
    coord_flip() +
    labs(title    = "Discriminative power by pathway category",
         subtitle = paste("Violin + boxplot of per-pathway AUC",
                          "·  dashed line = random classifier"),
         x = NULL, y = "AUC",
         caption = "n = number of pathways  ·  med = median AUC") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle      = element_text(colour = "grey40"),
          plot.caption       = element_text(colour = "grey40"),
          strip.text         = element_text(face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank())

  out <- file.path(plot_dir, "05_auc_by_category.png")
  ggsave(out, p, width = 14, height = 8, dpi = 150)
  cat(sprintf("    Saved: %s\n", out))
  invisible(p)
}


# ==============================================================================
# PLOT 5 — Gene sharing between pathways (Jaccard similarity heatmap)
# ==============================================================================

.plot_jaccard_heatmap <- function(psum, result, plot_dir, n_pathways = 80) {
  cat(sprintf("  [6/6] Jaccard heatmap (top %d pathways by size)...\n",
              n_pathways))

  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    stop("Package 'RColorBrewer' is required for the Jaccard heatmap: ",
         "install.packages('RColorBrewer')", call. = FALSE)

  kegg   <- result$kegg_genesets
  n_cats <- n_distinct(psum$category)
  per_cat <- max(3L, floor(n_pathways / n_cats))

  selected <- psum %>%
    group_by(category) %>%
    slice_max(n_genes_kegg, n = per_cat, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(category, desc(n_genes_kegg)) %>%
    pull(pathway)

  if (length(selected) < n_pathways) {
    extras <- psum %>%
      filter(!pathway %in% selected) %>%
      slice_max(n_genes_kegg, n = n_pathways - length(selected),
                with_ties = FALSE) %>%
      pull(pathway)
    selected <- c(selected, extras)
  }

  gene_sets <- kegg %>%
    filter(gs_name %in% selected) %>%
    group_by(gs_name) %>%
    summarise(genes = list(unique(ensembl_gene)), .groups = "drop") %>%
    deframe()

  pathways_sel <- names(gene_sets)
  n_sel        <- length(pathways_sel)

  cat(sprintf("    Computing %d \u00d7 %d Jaccard matrix...\n", n_sel, n_sel))

  jac_mat <- matrix(0, nrow = n_sel, ncol = n_sel,
                    dimnames = list(pathways_sel, pathways_sel))

  for (i in seq_len(n_sel)) {
    for (j in seq_len(n_sel)) {
      if (i == j) { jac_mat[i, j] <- 1; next }
      if (j < i)  { jac_mat[i, j] <- jac_mat[j, i]; next }
      a <- gene_sets[[i]]; b <- gene_sets[[j]]
      jac_mat[i, j] <- length(intersect(a, b)) / length(union(a, b))
    }
  }

  hc  <- hclust(as.dist(1 - jac_mat), method = "ward.D2")
  ord <- hc$labels[hc$order]

  label_map <- psum %>%
    filter(pathway %in% pathways_sel) %>%
    select(pathway, short_name, category)

  jac_long <- as.data.frame(as.table(jac_mat)) %>%
    rename(pathway_x = Var1, pathway_y = Var2, jaccard = Freq) %>%
    mutate(
      pathway_x = factor(pathway_x, levels = ord),
      pathway_y = factor(pathway_y, levels = ord)
    ) %>%
    left_join(label_map %>% rename(pathway_x = pathway,
                                   short_x   = short_name,
                                   cat_x     = category),
              by = "pathway_x") %>%
    left_join(label_map %>% rename(pathway_y = pathway,
                                   short_y   = short_name,
                                   cat_y     = category),
              by = "pathway_y")

  n_cats_local <- n_distinct(label_map$category)
  heatmap_pal  <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(8, "Set2")
  )(n_cats_local)
  cat_colours <- stats::setNames(
    heatmap_pal,
    unique(label_map$category[match(ord, label_map$pathway)])
  )

  p_heat <- ggplot(jac_long,
                   aes(x = pathway_x, y = pathway_y, fill = jaccard)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", mid = "#FDAE61", high = "#D73027",
                         midpoint = 0.25, limits = c(0, 1),
                         name = "Jaccard\nsimilarity") +
    scale_x_discrete(labels = stats::setNames(
      label_map$short_name[match(ord, label_map$pathway)], ord
    )) +
    scale_y_discrete(labels = stats::setNames(
      label_map$short_name[match(ord, label_map$pathway)], ord
    )) +
    labs(title    = sprintf(
      "Pairwise gene-set overlap between %d KEGG pathways", n_sel),
      subtitle = "Jaccard similarity  ·  rows and columns clustered by Ward's D2",
      x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.subtitle   = element_text(colour = "grey40"),
          axis.text.x     = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y     = element_text(size = 6),
          panel.grid      = element_blank(),
          legend.position = "right")

  out <- file.path(plot_dir, "06_gene_sharing_jaccard.png")
  ggsave(out, p_heat, width = 14, height = 13, dpi = 150)
  cat(sprintf("    Saved: %s\n", out))
  invisible(p_heat)
}


# ==============================================================================
# plot_kegg_characterisation() — main entry point
# ==============================================================================

#' Generate KEGG pathway characterisation plots
#'
#' Produces seven PNG files that characterise the KEGG pathway data used in
#' the network validation, without modifying any package internals.
#'
#' Output files (written to \code{plot_dir}):
#' \enumerate{
#'   \item \code{01_pathway_size_by_category.png}
#'   \item \code{02_network_coverage_by_category.png}
#'   \item \code{03a_size_vs_coverage_plain.png}
#'   \item \code{03b_size_vs_coverage_enriched.png}
#'   \item \code{04_coverage_vs_auc.png}
#'   \item \code{05_auc_by_category.png}
#'   \item \code{06_gene_sharing_jaccard.png}
#' }
#'
#' @param result    List returned by \code{compute_network_cutoff()}.
#' @param plot_dir  Character. Output directory. Default
#'   \code{"data-raw/plots/kegg_characterisation"}.
#' @param n_jaccard Integer. Pathways included in the Jaccard heatmap. Default
#'   \code{80}.
#'
#' @return Invisibly returns a named list of the six ggplot objects.
plot_kegg_characterisation <- function(
    result,
    plot_dir  = "data-raw/plots/kegg_characterisation",
    n_jaccard = 80L) {

  required <- c("kegg_genesets", "pathway_results", "gene_pairs")
  missing  <- setdiff(required, names(result))
  if (length(missing))
    stop("result is missing required elements: ",
         paste(missing, collapse = ", "), call. = FALSE)

  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  # 10-colour palette (extends Set2 which tops out at 8)
  cat_palette <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(8, "Set2")
  )(10L)

  cat("================================================================\n")
  cat("  KEGG Pathway Characterisation Plots\n")
  cat("================================================================\n\n")

  cat("  Building pathway summary table...\n")
  psum <- .build_pathway_summary(result)
  cat(sprintf("  %d pathways  ·  %d categories\n\n",
              nrow(psum), n_distinct(psum$category)))

  p1 <- .plot_pathway_sizes(psum, plot_dir, cat_palette)
  p2 <- .plot_network_coverage(psum, plot_dir, cat_palette)
  p3 <- .plot_size_vs_coverage(psum, plot_dir, cat_palette)
  p4 <- .plot_coverage_vs_auc(psum, plot_dir, cat_palette)
  p5 <- .plot_auc_by_category(psum, plot_dir, cat_palette)
  p6 <- .plot_jaccard_heatmap(psum, result, plot_dir,
                              n_pathways = n_jaccard)

  cat("\n================================================================\n")
  cat(sprintf("  Done. All plots saved to: %s\n", plot_dir))
  cat("================================================================\n\n")

  invisible(list(
    pathway_sizes            = p1,
    network_coverage         = p2,
    size_vs_coverage_plain   = p3$plain,
    size_vs_coverage_enriched = p3$enriched,
    coverage_vs_auc          = p4,
    auc_by_category          = p5,
    jaccard_heatmap          = p6
  ))
}


# ==============================================================================
# plot_gene_pathway_membership() — standalone function
# ==============================================================================

#' Plot the distribution of gene pathway-membership counts
#'
#' Three-panel figure:
#' \itemize{
#'   \item Panel A: histogram of how many genes appear in exactly 1, 2, 3 …
#'     pathways, coloured by membership class.
#'   \item Panel B: cumulative percentage of genes vs pathway-membership count.
#'   \item Panel C: top 25 highest-membership genes (lollipop), coloured by
#'     network presence, sized by number of biological categories spanned.
#' }
#'
#' @param result    List returned by \code{compute_network_cutoff()}.
#' @param plot_dir  Character. Output directory.
#'   Default \code{"data-raw/plots/kegg_characterisation"}.
#' @param filename  Character. Output filename.
#' @param x_cap     Integer. Maximum membership count shown on the histogram
#'   x-axis (long tail annotated). Default \code{30}.
#'
#' @return A ggplot object (the assembled cowplot grid), invisibly.
plot_gene_pathway_membership <- function(
    result,
    plot_dir = "data-raw/plots/kegg_characterisation",
    filename = "06_gene_pathway_membership.png",
    x_cap    = 30L) {

  cat("  [6] Gene pathway-membership frequency distribution...\n")

  kegg      <- result$kegg_genesets
  net_genes <- names(.build_neighbor_list(result$gene_pairs))

  gene_membership <- kegg %>%
    group_by(gene = ensembl_gene) %>%
    summarise(
      n_pathways   = n_distinct(gs_name),
      n_categories = n_distinct(.infer_category(gs_name)),
      .groups      = "drop"
    ) %>%
    mutate(
      in_network = gene %in% net_genes,
      hub_class  = case_when(
        n_pathways == 1  ~ "1  (singleton)",
        n_pathways <= 3  ~ "2-3",
        n_pathways <= 10 ~ "4-10",
        n_pathways <= 20 ~ "11-20",
        TRUE             ~ "20+"
      ),
      hub_class = factor(hub_class,
                         levels = c("1  (singleton)", "2-3",
                                    "4-10", "11-20", "20+"))
    )

  freq_tbl <- gene_membership %>%
    count(n_pathways, hub_class, name = "n_genes") %>%
    arrange(n_pathways)

  n_total    <- nrow(gene_membership)
  n_tail     <- sum(gene_membership$n_pathways > x_cap)
  med_memb   <- median(gene_membership$n_pathways)
  mean_memb  <- mean(gene_membership$n_pathways)
  pct_multi  <- round(100 * mean(gene_membership$n_pathways > 1), 1)

  hub_colours <- c(
    "1  (singleton)" = "#B0C4DE",
    "2-3"            = "#6BAED6",
    "4-10"           = "#2171B5",
    "11-20"          = "#FD8D3C",
    "20+"            = "#D94801"
  )

  # Panel A: histogram
  hist_data <- freq_tbl %>% filter(n_pathways <= x_cap)

  p_hist <- ggplot(hist_data,
                   aes(x = n_pathways, y = n_genes, fill = hub_class)) +
    geom_col(width = 0.8, colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = med_memb,  linetype = "solid",
               colour = "grey30", linewidth = 0.7) +
    geom_vline(xintercept = mean_memb, linetype = "dashed",
               colour = "grey30", linewidth = 0.7) +
    annotate("text", x = med_memb + 0.4, y = Inf,
             label = sprintf("median = %g", med_memb),
             hjust = 0, vjust = 1.5, size = 3, colour = "grey25") +
    annotate("text", x = mean_memb + 0.4, y = Inf,
             label = sprintf("mean = %.1f", mean_memb),
             hjust = 0, vjust = 3.2, size = 3, colour = "grey25") +
    { if (n_tail > 0)
      annotate("text", x = x_cap, y = max(hist_data$n_genes) * 0.6,
               label = sprintf("+ %d genes\nwith >%d pathways\n(not shown)",
                               n_tail, x_cap),
               hjust = 1, size = 2.8, colour = "grey40", fontface = "italic")
    } +
    scale_fill_manual(values = hub_colours, name = "Membership\nclass") +
    scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30),
                       limits = c(0.5, x_cap + 0.5)) +
    scale_y_continuous(labels = scales::comma) +
    labs(title    = "How many KEGG pathways does each gene belong to?",
         subtitle = sprintf(
           "%d genes total  ·  %g%% in more than one pathway  ·  x-axis capped at %d",
           n_total, pct_multi, x_cap),
         x = "Number of KEGG pathways", y = "Number of genes") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle    = element_text(colour = "grey40"),
          panel.grid.minor = element_blank(),
          legend.position  = "right")

  # Panel B: cumulative distribution
  cumul_data <- freq_tbl %>%
    arrange(n_pathways) %>%
    mutate(cumul_genes = cumsum(n_genes),
           cumul_pct   = 100 * cumul_genes / n_total)

  threshold_pts <- purrr::map_dfr(c(50, 75, 90), function(pct) {
    idx <- which(cumul_data$cumul_pct >= pct)[1]
    if (is.na(idx)) return(NULL)
    tibble::tibble(
      n_pathways = cumul_data$n_pathways[idx],
      cumul_pct  = cumul_data$cumul_pct[idx],
      label      = sprintf("%d%%\nat k=%d", pct, cumul_data$n_pathways[idx])
    )
  })

  p_cumul <- ggplot(cumul_data, aes(x = n_pathways, y = cumul_pct)) +
    geom_line(colour = "#2171B5", linewidth = 1) +
    geom_point(data = threshold_pts,
               aes(x = n_pathways, y = cumul_pct),
               colour = "#D94801", size = 3, inherit.aes = FALSE) +
    ggrepel::geom_label_repel(
      data = threshold_pts,
      mapping = aes(x = n_pathways, y = cumul_pct, label = label),
      size = 2.8, colour = "#D94801", box.padding = 0.3,
      min.segment.length = 0, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = c(50, 75, 90),
               linetype = "dotted", colour = "grey60", linewidth = 0.4) +
    scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30),
                       limits = c(1, x_cap)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       labels = function(x) paste0(x, "%")) +
    labs(title    = "Cumulative % of genes by membership count",
         subtitle = "What fraction of genes have membership \u2264 k?",
         x = "Pathway membership count (k)",
         y = "Cumulative % of genes") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle    = element_text(colour = "grey40"),
          panel.grid.minor = element_blank())

  # Panel C: top 25 hub genes lollipop
  top_genes <- gene_membership %>%
    slice_max(n_pathways, n = 25, with_ties = FALSE) %>%
    arrange(n_pathways) %>%
    mutate(gene = factor(gene, levels = gene))

  p_lollipop <- ggplot(top_genes,
                       aes(x = gene, y = n_pathways, colour = in_network)) +
    geom_segment(aes(xend = gene, y = 0, yend = n_pathways),
                 linewidth = 0.7) +
    geom_point(aes(size = n_categories)) +
    scale_colour_manual(
      values = c("TRUE" = "#2171B5", "FALSE" = "#BDBDBD"),
      labels = c("TRUE" = "In network", "FALSE" = "Not in network"),
      name   = NULL
    ) +
    scale_size_continuous(
      name   = "# biological\ncategories",
      range  = c(2, 7),
      breaks = c(1, 3, 5, 7)
    ) +
    coord_flip() +
    labs(title    = "Top 25 genes by KEGG pathway membership",
         subtitle = "Point size = number of distinct biological categories spanned",
         x = NULL, y = "Number of pathways") +
    theme_minimal(base_size = 12) +
    theme(plot.subtitle    = element_text(colour = "grey40"),
          panel.grid.minor = element_blank(),
          legend.position  = "right")

  # Assemble and save
  top_row    <- cowplot::plot_grid(p_hist, p_cumul, ncol = 2,
                                   rel_widths = c(1.4, 1),
                                   labels = c("A", "B"), label_size = 11)
  bottom_row <- cowplot::plot_grid(p_lollipop, ncol = 1,
                                   labels = "C", label_size = 11)
  final <- cowplot::plot_grid(top_row, bottom_row,
                              nrow = 2, rel_heights = c(1, 1.2))

  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(plot_dir, filename)
  ggsave(out, final, width = 14, height = 12, dpi = 150)
  cat(sprintf("    Saved: %s\n", out))

  invisible(final)
}
