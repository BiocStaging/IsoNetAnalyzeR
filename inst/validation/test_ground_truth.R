## test_ground_truth.R
## Recomputes, without any IsoNetAnalyzeR code, the values that the tests on
## the bundled data assert (tests/testthat/test-bundled-network.R and
## tests/testthat/test-duplicate-isoforms.R). It uses only the data files, the
## Ensembl 107 EnsDb read with ensembldb, and a plain igraph implementation of
## the virtual pulldown as documented in ?virtual_pulldown.
##
## Usage, from the package root:
##   git show v0.6.5:R/sysdata.rda > sysdata_0.6.5.rda
##   Rscript inst/validation/test_ground_truth.R sysdata_0.6.5.rda

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L)
    stop("Usage: Rscript test_ground_truth.R <0.6.5 sysdata.rda>")
suppressPackageStartupMessages(library(igraph))

n_failed <- 0L
check <- function(label, value, expected) {
    ok <- isTRUE(all.equal(value, expected, check.attributes = FALSE))
    if (!ok) n_failed <<- n_failed + 1L
    cat(sprintf("%-58s %-18s %s\n", label, paste(format(value), collapse = " "),
                if (ok) "agrees" else
                    paste("DIFFERS from", paste(expected, collapse = " "))))
}
pair_key <- function(g) {
    el <- as_edgelist(g)
    paste(pmin(el[, 1L], el[, 2L]), pmax(el[, 1L], el[, 2L]))
}

## Virtual pulldown as documented in ?virtual_pulldown: keep interact == 1
## edges, take the vertices within `order` hops of the seeds, keep the seeds
## and every other vertex with at least `cutoff` of its edges inside that
## candidate set, and return the subgraph induced by the kept vertices.
pulldown <- function(net, seeds, order, cutoff) {
    net <- subgraph_from_edges(net, which(E(net)$interact == 1L),
                               delete.vertices = FALSE)
    hops <- apply(distances(net, v = seeds), 2L, min)
    candidates <- V(net)$name[hops <= order]
    inside <- degree(induced_subgraph(net, candidates))[candidates]
    total <- degree(net, v = candidates)
    keep <- candidates[candidates %in% seeds | inside / total >= cutoff]
    induced_subgraph(net, keep)
}

## -- The 0.6.5 network, the source of the bundled one -------------------------
full <- new.env()
load(args[[1L]], envir = full)
predicted <- subgraph_from_edges(full$.network_data,
                                 which(E(full$.network_data)$interact == 1L),
                                 delete.vertices = FALSE)
pairs <- pair_key(predicted)
self_pair <- ends(predicted, E(predicted))[, 1L] ==
    ends(predicted, E(predicted))[, 2L]
check("0.6.5: transcripts", vcount(predicted), 58164L)
check("0.6.5: predicted interaction edges", ecount(predicted), 158458L)
check("0.6.5: pairs listed twice", sum(duplicated(pairs) & !self_pair), 6L)
check("0.6.5: transcripts paired with themselves", sum(self_pair), 5L)
check("0.6.5: distinct pairs of two transcripts",
      length(unique(pairs[!self_pair])), 158447L)
prpf <- paste("ENST00000304992", "ENST00000378572")
check("0.6.5: pi of the two PRPF18-PRPF8 copies",
      sort(E(predicted)$pi[pairs == prpf]), c(0.5594, 0.5911))

## -- The bundled network ----------------------------------------------------------
bundled <- new.env()
load(file.path("R", "sysdata.rda"), envir = bundled)
g <- bundled$.network_data
check("bundled: same transcripts as 0.6.5",
      identical(V(g)$name, V(predicted)$name), TRUE)
check("bundled: edges are the distinct 0.6.5 pairs",
      setequal(pair_key(g), unique(pairs[!self_pair])) &&
          !anyDuplicated(pair_key(g)), TRUE)
check("bundled: neighbourhood cutoff", round(bundled$.network_cutoff, 7L),
      0.2608696)

## -- Ensembl 107 annotation (AnnotationHub AH104864) ---------------------------
edb <- AnnotationHub::AnnotationHub()[["AH104864"]]
tx <- ensembldb::transcripts(edb, columns = c("tx_id", "tx_biotype"),
                             return.type = "data.frame")
biotype <- tx$tx_biotype[match(V(g)$name, tx$tx_id)]
check("Ensembl 107: network transcripts without biotype", sum(is.na(biotype)),
      0L)
check("Ensembl 107: protein_coding / NMD / NSD",
      as.vector(table(biotype)[c("protein_coding", "nonsense_mediated_decay",
                                 "non_stop_decay")]),
      c(46433L, 11662L, 69L))
eligible <- induced_subgraph(g, which(!biotype %in% c(
    "nonsense_mediated_decay", "non_stop_decay")))
check("eligible network: transcripts", vcount(eligible), 46433L)
check("eligible network: edges", ecount(eligible), 146256L)
genes <- ensembldb::genes(edb, columns = c("gene_id", "gene_name"),
                          return.type = "data.frame")
check("Ensembl 107: gene ID of CLTCL1",
      genes$gene_id[genes$gene_name == "CLTCL1"], "ENSG00000070371")
check("Ensembl 107: gene ID of APAF1",
      genes$gene_id[genes$gene_name == "APAF1"], "ENSG00000120868")

## -- Example seeds and protein-identical map -----------------------------------
seeds <- utils::read.csv(file.path("inst", "extdata",
                                   "gse141432_t2d_vs_lean_seeds.csv"))
check("seeds: rows", nrow(seeds), 137L)
check("seeds: genes", length(unique(seeds$gene_id)), 65L)
map <- utils::read.csv(file.path("inst", "extdata",
                                 "protein_identical_isoform_map.csv"))
check("map: rows", nrow(map), 16630L)
check("map: distinct kept transcripts", length(unique(map$kept_ENST)), 9760L)
check("map: distinct redundant transcripts",
      length(unique(map$redundant_ENST)), 16630L)
check("map: kept transcript in the network", sum(map$kept_ENST %in% V(g)$name),
      10592L)
check("map: kept transcript not in the network",
      sum(!map$kept_ENST %in% V(g)$name), 6038L)
check("map: redundant transcripts that are also kept",
      sum(map$redundant_ENST %in% map$kept_ENST), 0L)

## -- Example neighbourhoods -------------------------------------------------------
seed_of <- function(gene, direction)
    seeds$isoform_id[seeds$gene_name == gene & seeds$direction == direction]
example <- c(seed_of("CLTCL1", "up"), seed_of("CLTCL1", "down"),
             seed_of("APAF1", "up"), seed_of("APAF1", "down"))
check("example seeds are eligible network transcripts",
      all(example %in% V(eligible)$name), TRUE)
check("example seeds are not redundant transcripts",
      any(example %in% map$redundant_ENST), FALSE)
check("CLTCL1 up: edges at cutoff 0.2609",
      ecount(pulldown(eligible, seed_of("CLTCL1", "up"), 1L, 0.2609)), 8L)
check("CLTCL1 down: edges at cutoff 0.2609",
      ecount(pulldown(eligible, seed_of("CLTCL1", "down"), 1L, 0.2609)), 33L)
apaf1_up <- pulldown(eligible, seed_of("APAF1", "up"), 1L,
                     bundled$.network_cutoff)
apaf1_down <- pulldown(eligible, seed_of("APAF1", "down"), 1L,
                       bundled$.network_cutoff)
check("APAF1 up: transcripts, edges",
      c(vcount(apaf1_up), ecount(apaf1_up)), c(3, 2))
check("APAF1 up: DIABLO ENST00000464942 included",
      "ENST00000464942" %in% V(apaf1_up)$name, FALSE)
check("APAF1 down: transcripts, edges",
      c(vcount(apaf1_down), ecount(apaf1_down)), c(9, 11))

cat(sprintf("\n%d value(s) differ from the tests.\n", n_failed))
