# IsoNetAnalyzeR data provenance

## Protein-identical isoform map

`inst/extdata/protein_identical_isoform_map.csv` is built from
`data-raw/duplicates_long_kept_redundant_enst_ensp_ensg.csv`, supplied by the
network authors. The source SHA-256 is
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
