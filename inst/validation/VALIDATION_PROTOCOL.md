# IsoNetAnalyzeR 0.99.0 validation protocol

## Automated package checks

1. Run the complete testthat suite from a clean session.
2. Build the source package with the vignette enabled.
3. Install the source tarball into an empty temporary library.
4. Run `R CMD check --no-manual` and require `Status: OK`.
5. Run `R CMD check --as-cran`; archive the full log and explain any external
   incoming-feasibility notes separately from package defects.
6. Confirm version 0.99.0, the canonical tarball SHA-256, and the
   `R/sysdata.rda` SHA-256 recorded in `inst/PROVENANCE.md`.
7. Run `BiocCheck::BiocCheck(new-package = TRUE)` on the source tarball.
8. Whenever the network file changes, rerun the network equivalence
   comparison against the previous release and archive its report.

## Required regression cases

- Default exclusion of NMD and NSD vertices and incident edges, independent
  inclusion flags, transcript-specific custom annotations and unknown-biotype
  rejection. An excluded alias must not bypass the policy by remapping.
- Retained-significant-seed validation after biotype filtering, identifier
  resolution and mapping. Require both mapped directions and one significant
  retained seed; allow nonsignificant complementary and mapped isolated seeds.
  Invalid p-values and excluded significant drivers must be rejected.
- Differential overlay plots distinguish purple focal seeds from directional
  neighbour membership, including the shared blue-neighbour case.
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
  `get_metadata(obj, "analysis_network_coverage")[[label]]`.
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
- Execution of the vignette and every man-page example.

## Scientific invariants

- The embedded network is the interaction-only reduction of the 0.6.5 file
  (`inst/scripts/reduce_network.R`). It keeps every vertex and every
  `interact == 1` edge, so analysis results are identical to 0.6.5; only
  counts of all scored edges change.
- The complete duplicate map is used for count aggregation independently of
  network applicability. Count aggregation preserves per-gene and global totals
  but requires a complete downstream statistical refit.
- Under identical seed inputs and eligibility policy, 0.6.5 analysis graphs
  retain the 0.6.0 semantics; the newer plotting palette changes presentation.
  Rebuilding a historical unfiltered session with the newer default decay
  policy can change vertices, edges and internal-edge fraction denominators.
- The packaged biotype and protein-identical maps retain their pinned checksums.
  Network eligibility does not change upstream RNA counts, isoform fractions,
  statistical tests or multiple-testing adjustments.
- The fields `n_edges_gained` and `n_edges_lost` retain their public names but
  are interpreted only as edges unique to the up-seed and down-seed static
  neighbourhoods.
- A mapped isolated transcript is a valid one-node result in the pinned
  network. It is not evidence of complete biological interaction loss.
- `mixed` and topology classes are descriptive network summaries, not
  mechanistic or disease classifications.

## Publication boundary

The package ships one example seed table (GSE141432, obese type 2 diabetic
versus lean adipose tissue) for the vignette and examples. Its origin and
selection rules are documented in `inst/scripts/`. Study-wide sensitivity
analyses, candidate lists and golden outputs remain in the companion analysis
and are validated separately. PairedGSEA redevelopment remains out of scope.

## Release artifacts

- one canonical source tarball and SHA-256;
- embedded-network checksum matching `inst/PROVENANCE.md`;
- test, build, clean-install, and check logs;
- annotated local `v0.99.0` tag pointing at the validated release commit;
- release notes and manual publication instructions; remote publication is a
  separate action by the release maintainer.
