## make_example_data.R
## Build the example data shipped in inst/extdata for examples and the
## vignette.
##
## Usage, from the package root:
##   Rscript inst/scripts/make_example_data.R <annotated_seeds.csv>
##
## 1. gse141432_t2d_vs_lean_seeds.csv
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
## 2. msigdb_hallmark_ensembl.rds
##    MSigDB Hallmark gene sets (Liberzon et al. 2015, Cell Systems 1:417)
##    as a named list of Ensembl gene identifiers, retrieved with msigdbr.
##    MSigDB Hallmark gene sets are distributed under CC BY 4.0.

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

hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- hallmark[!is.na(hallmark$ensembl_gene) &
                     nzchar(hallmark$ensembl_gene), ]
pathways <- lapply(split(hallmark$ensembl_gene, hallmark$gs_name), unique)
attr(pathways, "source") <- list(
    collection = "MSigDB Hallmark (H)",
    db_version = unique(hallmark$db_version),
    msigdbr = as.character(utils::packageVersion("msigdbr")),
    identifier = "Ensembl gene",
    licence = "CC BY 4.0"
)
out_hallmark <- file.path("inst", "extdata", "msigdb_hallmark_ensembl.rds")
saveRDS(pathways, out_hallmark, compress = "xz")

cat("Input seeds SHA-256 :", sha256(args[[1L]]), "\n")
cat(out_seeds, ":", nrow(seeds), "rows,", length(unique(seeds$gene_id)),
    "genes, SHA-256", sha256(out_seeds), "\n")
cat(out_hallmark, ":", length(pathways), "sets, MSigDB",
    attr(pathways, "source")$db_version, ", SHA-256", sha256(out_hallmark),
    "\n")
