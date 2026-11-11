## make_example_data.R
## Build the example seed table shipped in inst/extdata for examples and the
## vignette.
##
## Usage, from the package root:
##   Rscript inst/scripts/make_example_data.R <annotated_seeds.csv>
##
## gse141432_t2d_vs_lean_seeds.csv
##    Isoform seeds for the Obese-T2D versus Lean contrast of GEO series
##    GSE141432 (human white adipose tissue; 8 lean, 8 obese normoglycaemic
##    and 8 obese type 2 diabetic donors). The input is the annotated seed table
##    of the IsoNetAnalyzeR 0.6.0 retained-significant-seed analysis
##    (analysis_id gse141432_retained_significant_seed_v060_2026_08_28,
##    contrast Obese_T2D_minus_Lean, dIF floor 0.10). Genes passed a
##    gene-level FDR of 0.05 with a maximum absolute dIF of at least 0.1;
##    transcripts were then selected by a within-gene Benjamini-Hochberg rule,
##    requiring a significant retained seed and mapped seeds in both
##    directions. Only the columns needed by the package are kept.
##
## GEO holds only a gene-level count matrix for GSE141432
## (GSE141432_countmatrix.txt.gz), so transcript usage was quantified from the
## sequencing reads upstream of this package; the seed table cannot be
## obtained with GEOquery.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L)
    stop("Usage: Rscript make_example_data.R <annotated_seeds.csv>")

sha256 <- function(path) digest::digest(file = path, algo = "sha256")
expected_seeds <- paste0(
    "1e054bf6490db9d8042f0060d789da64214add7a055e73fe85770fcb082bd561"
)
if (!identical(sha256(args[[1L]]), expected_seeds))
    stop("Input is not the frozen GSE141432 annotated seed table.")

seeds <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE)
cols <- c("gene_id", "gene_name", "isoform_id", "direction", "signed_dif",
          "within_gene_padj")
stopifnot(all(cols %in% names(seeds)), nrow(seeds) == 137L,
          length(unique(seeds$gene_id)) == 65L)
seeds <- seeds[, cols]
seeds$signed_dif <- signif(seeds$signed_dif, 6L)
seeds$within_gene_padj <- signif(seeds$within_gene_padj, 6L)
out_seeds <- file.path("inst", "extdata", "gse141432_t2d_vs_lean_seeds.csv")
utils::write.csv(seeds, out_seeds, row.names = FALSE)

cat("Input seeds SHA-256 :", sha256(args[[1L]]), "\n")
cat(out_seeds, ":", nrow(seeds), "rows,", length(unique(seeds$gene_id)),
    "genes, SHA-256", sha256(out_seeds), "\n")
