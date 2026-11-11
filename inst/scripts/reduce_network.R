## reduce_network.R
## Derive the bundled interaction-only network from the full scored network
## shipped in IsoNetAnalyzeR 0.6.5.
##
## The 0.6.5 R/sysdata.rda stored every scored isoform pair (2,853,266 edges),
## but all analysis code first restricts the graph to predicted interaction
## edges (interact == 1). This script keeps those edges and all 58,164
## vertices, so isoforms without a predicted interaction remain
## mapped-isolated rather than becoming unmapped.
##
## The source table lists six transcript pairs twice, with slightly
## different pi, and five transcripts paired with themselves. Each repeated
## pair becomes one edge carrying the mean pi of its copies, and the
## self-pairs are removed, giving a simple graph of 158,447 edges. Every
## other edge keeps its position and attributes. Vertex order, vertex
## attributes, the neighbourhood cutoff and the pi threshold are kept.
##
## Usage, from the package root:
##   Rscript inst/scripts/reduce_network.R <0.6.5 sysdata.rda> R/sysdata.rda
##
## The 0.6.5 input is available from the v0.6.5 git tag:
##   git show v0.6.5:R/sysdata.rda > sysdata_0.6.5.rda
##
## The full network was built by data-raw/update_network.R from the
## Isoform Analysis group's predicted interaction table
## (results_PHYSICAL_700_FULL_threshold_22_06_26.csv); see inst/PROVENANCE.md.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
    stop("Usage: Rscript reduce_network.R <input sysdata.rda> <output .rda>")
input <- args[[1L]]
output <- args[[2L]]

sha256 <- function(path) digest::digest(file = path, algo = "sha256")

expected_input <- paste0(
    "3130209f090b9b83bbf2ebf328f29a0ad9a8fcd89922a277d0966325bc163ebb"
)
if (!identical(sha256(input), expected_input))
    stop("Input is not the IsoNetAnalyzeR 0.6.5 R/sysdata.rda.")

full <- new.env()
load(input, envir = full)
g <- full$.network_data

keep <- which(igraph::edge_attr(g, "interact") == 1L)
predicted <- igraph::subgraph_from_edges(g, keep, delete.vertices = FALSE)

el <- igraph::as_edgelist(predicted, names = FALSE)
pair <- paste(pmin(el[, 1L], el[, 2L]), pmax(el[, 1L], el[, 2L]))
self_pair <- el[, 1L] == el[, 2L]
repeated <- duplicated(pair)
igraph::E(predicted)$pi <- stats::ave(igraph::E(predicted)$pi, pair,
                                      FUN = mean)
.network_data <- igraph::delete_edges(predicted,
                                      which(repeated | self_pair))
.network_cutoff <- full$.network_cutoff
.pi_threshold <- full$.pi_threshold

stopifnot(
    identical(igraph::V(.network_data)$name, igraph::V(g)$name),
    identical(igraph::V(.network_data)$gene, igraph::V(g)$gene),
    all(igraph::E(.network_data)$interact == 1L),
    igraph::is_simple(.network_data),
    sum(repeated & !self_pair) == 6L, sum(self_pair) == 5L,
    igraph::ecount(.network_data) ==
        length(unique(pair[!self_pair]))
)

save(.network_data, .network_cutoff, .pi_threshold,
     file = output, compress = "xz")

cat(sprintf("Input  : %d vertices, %d edges, SHA-256 %s\n",
            igraph::vcount(g), igraph::ecount(g), expected_input))
cat(sprintf("Predicted interaction edges: %d (%d repeated pairs merged, %d self-pairs removed)\n",
            length(keep), sum(repeated & !self_pair), sum(self_pair)))
cat(sprintf("Output : %d vertices, %d edges, SHA-256 %s\n",
            igraph::vcount(.network_data), igraph::ecount(.network_data),
            sha256(output)))
cat(sprintf("network_cutoff = %.7f; pi_threshold = %s\n", .network_cutoff,
            if (is.null(.pi_threshold)) "NULL" else .pi_threshold))
