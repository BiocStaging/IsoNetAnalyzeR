# IsoNetAnalyzeR data provenance

## Interaction network

`R/sysdata.rda` holds the bundled isoform interaction network
(`.network_data`), the neighbourhood cutoff (`.network_cutoff` = 0.2608696)
and the pi threshold (`.pi_threshold`, `NULL`: the `interact` flag was taken
from the source table).

The interactions derive from STRING physical protein-protein interactions
with a combined score of at least 700 (Szklarczyk et al. 2023, "The STRING
database in 2023", Nucleic Acids Res 51:D638-D646, doi:10.1093/nar/gkac1000).
STRING reports interactions between genes. A Bayesian latent-variable model
developed by Isabella Østerlund (Isoform Analysis Group, DTU) treats each
gene-level interaction as evidence that at least one isoform pair of the two
genes interacts, and estimates each isoform pair's interaction probability
(`pi`) by MCMC sampling from the posterior. For each gene pair, the n isoform
pairs with the highest `pi` are marked as interacting (`interact = 1`), where
n is the sum of the gene pair's `pi` values rounded to the nearest integer.
STRING data are distributed under the Creative Commons Attribution 4.0
licence (CC BY 4.0, <https://string-db.org/cgi/download>); this section is
the attribution for the network derived from them.

The network was built by `data-raw/update_network.R` from the Isoform
Analysis group's predicted isoform interaction table
`results_PHYSICAL_700_FULL_threshold_22_06_26.csv` (columns `enst_1`,
`gene_1`, `enst_2`, `gene_2`, `pi`). Up to version 0.6.5 the package shipped
every scored pair: 58,164 transcripts and 2,853,266 edges, of which 158,458
are predicted interactions (`interact == 1`). That file has SHA-256
`3130209f090b9b83bbf2ebf328f29a0ad9a8fcd89922a277d0966325bc163ebb` and remains
available from the v0.6.5 tag.

From version 0.99.0, `inst/scripts/reduce_network.R` keeps the 158,458
predicted interaction edges and all 58,164 transcripts, with their `name`,
`gene`, `pi` and `interact` attributes and the original vertex and edge
order. Every analysis restricts the network to `interact == 1` edges before
use, so results are unchanged. This was verified against 0.6.5 for the
interaction graph itself, 65 GSE141432 genes at six order/cutoff settings,
500 random seed sets, lookups, components, enrichment and diagnostics. Only
reported counts of all scored edges change: the session edge count, the
eligibility audit's `edges_before`/`edges_after`, and the diagnostics field
`full_edges`. The reduced file has SHA-256
`41df4ce0aa5f83878382238de5444697a3d7d9ebc733d05fd8a1b7434c2abf9b`.

## Gene symbols

`inst/extdata/gene_names_gencode_v41.csv.gz` lists the GENCODE v41 gene name
(usually the HGNC symbol) of each of the 10,554 genes in the bundled network.
It comes from the same GENCODE v41 isoform annotation as the transcript
biotypes below (source object SHA-256
`b4c6225bd65a9100670e9b92ddf86ea2eeb57abe3bf020b01f939452c9686931`) and is
written by `inst/scripts/make_gene_names.R`. Gene IDs are stored without
version suffix. The exported table has SHA-256
`33db6f5a982389191eb26d378d638948a4fb858e55e255a308fd73d48d1a59a7`.
The symbols only label plots and tables; no analysis uses them.

## Example data

`inst/extdata/gse141432_t2d_vs_lean_seeds.csv` and
`inst/extdata/msigdb_hallmark_ensembl.rds` are produced by
`inst/scripts/make_example_data.R`; see `inst/scripts/README.md`.

## Protein-identical isoform map

`inst/extdata/protein_identical_isoform_map.csv` is built from
`data-raw/duplicates_long_kept_redundant_enst_ensp_ensg.csv`, supplied by
Isabella Østerlund. The source SHA-256 is
`cfe6fb57ed06c7231b9c77d7246b148be3cfdb1ac6aad4f9895cc422efaa36f5`.

The source contains 16,630 unique redundant ENSTs mapping unambiguously to
9,760 kept ENST representatives. It contains no missing identifiers, duplicate
rows, transitive kept/redundant chains, or redundant-to-multiple-kept mappings.
The supplied table has no redundant-gene column. `redundant_ENSG` is therefore
derived as `kept_ENSG`; every kept and redundant transcript was present in the
checksum-pinned GSE141432 `SummarizedExperiment`, with zero gene mismatches.
The evidence SE SHA-256 is
`8ab52523b1564b061c9424e5189dcf6014f0f4e0161d1eb1e6f87c855dbe7b35`.

The map is transcriptome-wide. Of its rows, 10,592 have a kept transcript in
the embedded IsoNetAnalyzeR network and 6,038 do not. All redundant transcripts
are absent from that network, and every applicable kept transcript agrees with
the network's gene attribute. Transcriptome-only rows remain in the package for
count aggregation.

The map was validated against the 0.6.5 `R/sysdata.rda` (SHA-256
`3130209f090b9b83bbf2ebf328f29a0ad9a8fcd89922a277d0966325bc163ebb`). The
0.99.0 network has the same vertices, so network applicability is unchanged. Machine-readable build
details and the generated map checksum are stored in
`inst/extdata/protein_identical_isoform_map_provenance.csv`.


## Transcript biotypes for decay exclusion

`inst/extdata/transcript_biotypes_gencode_v41.csv.gz` contains 250,161 unique
stable transcript IDs and transcript-level biotypes derived from the local
GENCODE v41 `switchList_v41$isoformFeatures` annotation supplied for this project.
Source object SHA-256:
`b4c6225bd65a9100670e9b92ddf86ea2eeb57abe3bf020b01f939452c9686931`.
The exported compressed table SHA-256 is
`f56d86c8447410de0de802f6b58d462e483deca7191b47191711989168ad8a84`.

All 58,164 raw bundled network vertices match this annotation: 46,433 are
`protein_coding`, 11,662 `nonsense_mediated_decay`, and 69 `non_stop_decay`.
This table is independent of Isabella's protein-identical deduplication map.
Version suffixes are normalized for matching; the pinned annotation is not a
claim that transcript structures remain identical across annotation releases.
The raw bundled network is unchanged. Filtering is applied to each new session.
