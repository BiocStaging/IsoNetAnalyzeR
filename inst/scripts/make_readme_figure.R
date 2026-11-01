## make_readme_figure.R
## Draw the CLTCL1 figure shown in README.md (man/figures/README-cltcl1.png).
##
## It repeats the vignette's analysis of the bundled GSE141432 example seeds
## and plots CLTCL1 the way the vignette does.
##
## Usage, from the package root, with the package installed:
##   Rscript inst/scripts/make_readme_figure.R man/figures/README-cltcl1.png

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L)
    stop("Usage: Rscript make_readme_figure.R <output.png>")
output <- args[[1L]]

library(IsoNetAnalyzeR)

seeds <- read.csv(system.file("extdata", "gse141432_t2d_vs_lean_seeds.csv",
                              package = "IsoNetAnalyzeR"))
obj <- IsoNetAnalyzeR("GSE141432")
obj <- run_differential_pulldown(
    obj, "t2d", seeds,
    order = 1, confirm = FALSE, verbose = FALSE,
    seed_padj_col = "within_gene_padj", seed_fdr = 0.05
)
cltcl1 <- get_gene_result(obj, "t2d", "ENSG00000070371")
p <- plot_network(cltcl1, direction = "overlay", level = "gene",
                  node_size = 3.5, label_size = 3)

set.seed(1)  # ggrepel places labels at random when the plot is drawn
ggplot2::ggsave(output, plot = p, width = 7.5, height = 6, dpi = 150,
                bg = "white")
message("Wrote ", output)
