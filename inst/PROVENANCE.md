# IsoNetAnalyzeR data provenance

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

The map was validated against `R/sysdata.rda` SHA-256
`3130209f090b9b83bbf2ebf328f29a0ad9a8fcd89922a277d0966325bc163ebb`.
The embedded network is not modified by this feature. Machine-readable build
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
