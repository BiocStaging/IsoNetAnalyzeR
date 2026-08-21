# IsoNetAnalyzeR 0.5.0 validation protocol

## Automated package checks

1. Run the complete testthat suite from a clean session.
2. Build the source package with the vignette enabled.
3. Install the source tarball into an empty temporary library.
4. Run `R CMD check --no-manual` and require `Status: OK`.
5. Run `R CMD check --as-cran`; archive the full log and explain any external
   incoming-feasibility notes separately from package defects.
6. Confirm version 0.5.0, the canonical tarball SHA-256, and the unchanged
   `R/sysdata.rda` SHA-256.

## Required regression cases

- Installed loading of the 16,630-row packaged map, its schema and provenance,
  9,760 representatives, and 10,592 network-applicable mappings.
- Version normalization, direct-vertex precedence, compatible and incompatible
  custom graphs, passthrough, and `resolve_duplicates = FALSE`.
- Same-direction representative collapse, opposite-direction whole-gene
  exclusion, all-genes-excluded failure, submitted-versus-resolved counts,
  resolved dominance scores, and persistence of all public API audits.
- Matrix/data-frame and single-assay `SummarizedExperiment` aggregation,
  absent kept rows, multiple redundants, deterministic metadata selection,
  invalid counts, overflow, feature ordering, and gene/library sum invariants.

- Three-state coverage overall and per supplied direction: unmapped, mapped
  without a predicted edge, and mapped with a predicted edge.
- Diagnostic low-mapping warnings, typed returned coverage, and persistence in
  `obj@metadata$analysis_network_coverage[[label]]`.
- Missing identifiers/directions, exact duplicates, direction and gene
  conflicts, network mapping disagreement, fractional integers, non-finite
  values, and invalid stored cutoffs.
- Partially missing multi-isoform lookups.
- Directed, unnamed, gene-incomplete, non-binary, and all-isolated custom
  networks.
- Clear reporting when every gene raises a genuine processing error.
- Valid fully unmapped, mapped-isolated, and one-direction results.
- ENSG bookkeeping at gene level and ENST bookkeeping at isoform level.
- Connected and singleton components at isoform and gene levels, including
  filtered component IDs.
- Explicit seed flags for disconnected multi-seed graphs and legacy saved
  node-stat tables without the flag.
- Plot assignment and `ggsave()` on explicit output paths without creating
  `Rplots.pdf`.
- Execution of the synthetic vignette workflow.

## Scientific invariants

- The embedded network file and its checksum are unchanged from repaired
  commit `f24b467`.
- The complete duplicate map is used for count aggregation independently of
  network applicability. Count aggregation preserves per-gene and global totals
  but requires a complete downstream statistical refit.
- Existing repaired graph inventories remain unchanged except where malformed
  inputs are now rejected before analysis.
- The fields `n_edges_gained` and `n_edges_lost` retain their public names but
  are interpreted only as edges unique to the up-seed and down-seed static
  neighbourhoods.
- A mapped isolated transcript is a valid one-node result in the pinned
  network. It is not evidence of complete biological interaction loss.
- `mixed` and topology classes are descriptive network summaries, not
  mechanistic or disease classifications.

## Publication boundary

No publication dataset, seed-selection method, threshold, candidate list,
sensitivity result, or golden output is included in the package. Publication
reproduction is owned by the companion analysis and is validated separately.
PairedGSEA redevelopment remains out of scope.

## Release artifacts

- one canonical source tarball and SHA-256;
- unchanged embedded-network checksum;
- test, build, clean-install, and check logs;
- annotated `v0.5.0` tag pointing at the validated release commit.
