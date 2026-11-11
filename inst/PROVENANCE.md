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

From version 0.99.1 the bundled network is a simple graph. The source table
lists six transcript pairs twice, with `pi` values up to 0.032 apart
(PRPF18-PRPF8, and DIABLO with BIRC2, APAF1 and another DIABLO transcript),
and five DIABLO transcripts paired with themselves. Because the neighbourhood
cutoff compares vertex degrees, which count a repeated edge twice and a
self-loop as two, these entries changed the pulldowns that reached them: in
the GSE141432 example, the APAF1 up-labelled neighbourhood at order 1
included a DIABLO transcript only through the repeated edges.
`inst/scripts/reduce_network.R` now keeps one edge per pair with the mean
`pi` of its copies and removes the self-pairs, leaving 158,447 edges. All
other edges keep their position and attributes. The file has SHA-256
`c2a0eeb08c60e566b839204473002a582a38b6d7795eb6b67410da3ca3597b52`.
Supplied networks are simplified the same way when a session is created.

## Example data

`inst/extdata/gse141432_t2d_vs_lean_seeds.csv` is produced by
`inst/scripts/make_example_data.R`; see `inst/scripts/README.md`. It is the
output of a differential transcript usage analysis of GEO series GSE141432.
GEO holds only a gene-level count matrix for that series, so the table
cannot be retrieved with GEOquery. Gene sets for the examples and the
vignette come from the msigdb package (MSigDB 2023.1, ExperimentHub).

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


## Transcript annotation for decay exclusion and gene symbols

Transcript biotypes and gene symbols are not shipped. They come from Ensembl
release 107, the gene build of GENCODE v41, as the EnsDb in AnnotationHub
record AH104864 ("Ensembl 107 EnsDb for Homo sapiens"), loaded when a
session needs them and cached by AnnotationHub.

The bundled network was built on GENCODE v41 transcripts. Up to version
0.99.0 the package shipped a GENCODE v41 biotype table (SHA-256
`f56d86c8447410de0de802f6b58d462e483deca7191b47191711989168ad8a84`) and a
gene-symbol table (SHA-256
`33db6f5a982389191eb26d378d638948a4fb858e55e255a308fd73d48d1a59a7`), both
exported from the project's GENCODE v41 `switchList_v41$isoformFeatures`
annotation (SHA-256
`b4c6225bd65a9100670e9b92ddf86ea2eeb57abe3bf020b01f939452c9686931`).
All 58,164 network transcripts are in AH104864 with the same biotype and
gene as in that table: 46,433 `protein_coding`, 11,662
`nonsense_mediated_decay` and 69 `non_stop_decay`, so decay exclusion is
unchanged. The gene symbols agree for 10,496 of the 10,554 network genes;
the other 58 have no symbol in Ensembl, and their gene ID is shown instead,
as before. Version suffixes are normalized for matching; the annotation is
not a claim that transcript structures remain identical across releases.
Filtering is applied to each new session; the bundled network is unchanged.
