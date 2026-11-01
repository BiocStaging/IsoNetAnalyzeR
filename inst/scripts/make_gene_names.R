## make_gene_names.R
## Write the GENCODE v41 gene symbols used to label plots
## (inst/extdata/gene_names_gencode_v41.csv.gz).
##
## The source is the project's GENCODE v41 isoform annotation
## (switchList_v41.Rdata, an IsoformSwitchAnalyzeR switchAnalyzeRlist), the
## same object that produced transcript_biotypes_gencode_v41.csv.gz; see
## inst/PROVENANCE.md. Gene IDs lose their version suffix and the _PAR_Y tag
## of pseudoautosomal copies, so every stable gene ID has one symbol. Only
## genes of the bundled network are kept; user-supplied networks can carry
## their own gene_name vertex attribute.
##
## Usage, from the package root:
##   Rscript inst/scripts/make_gene_names.R <switchList_v41.Rdata> \
##       R/sysdata.rda inst/extdata/gene_names_gencode_v41.csv.gz

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
    stop("Usage: Rscript make_gene_names.R <switchList_v41.Rdata> ",
         "<sysdata.rda> <output>")
input <- args[[1L]]
sysdata <- args[[2L]]
output <- args[[3L]]

sha256 <- function(path) digest::digest(file = path, algo = "sha256")

expected_input <- paste0(
    "b4c6225bd65a9100670e9b92ddf86ea2eeb57abe3bf020b01f939452c9686931"
)
if (!identical(sha256(input), expected_input))
    stop("Input is not the project's GENCODE v41 switchList_v41.Rdata.")

env <- new.env()
objects <- load(input, envir = env)
stopifnot(length(objects) == 1L)
features <- env[[objects]]$isoformFeatures

genes <- unique(data.frame(
    gene_id = sub("\\..*$", "", features$gene_id),
    gene_name = trimws(features$gene_name)
))
genes <- genes[!is.na(genes$gene_name) & nzchar(genes$gene_name), ]
stopifnot(!anyDuplicated(genes$gene_id))

bundled <- new.env()
load(sysdata, envir = bundled)
network_genes <- unique(sub("\\..*$", "",
                            igraph::V(bundled$.network_data)$gene))
stopifnot(all(network_genes %in% genes$gene_id))
genes <- genes[genes$gene_id %in% network_genes, ]
genes <- genes[order(genes$gene_id), ]

con <- gzfile(output, "w")
utils::write.csv(genes, con, row.names = FALSE)
close(con)

message(sprintf("Wrote %d gene symbols to %s (SHA-256 %s).",
                nrow(genes), output, sha256(output)))
