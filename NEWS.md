# IsoNetAnalyzeR 0.5.5

- New sessions exclude annotated NMD and NSD transcripts from protein networks by default. Independent `include_nmd`/`include_nsd` flags restore either class explicitly.
- Pinned transcript-biotype annotations and exclusion audits distinguish biological eligibility from statistical seed selection and network coverage.
- Custom networks require known transcript biotypes when exclusions are enabled. Existing saved sessions are unchanged; rebuild them to apply the new policy.
- No maximum seed count is imposed by the package; upstream studies should declare effect/significance rules and whether both directions are required.

# IsoNetAnalyzeR 0.5.0

## Protein-identical isoform support

- Ship a provenance-checked, transcriptome-wide map from 16,630 redundant
  transcript identifiers to 9,760 protein representatives. The complete map is
  available through `isonet_duplicate_map()`; 10,592 mappings are applicable
  to the embedded network and 6,038 additional mappings remain available for
  upstream count aggregation.
- Add `resolve_duplicate_isoforms()` for version normalization and
  graph-compatible duplicate resolution. Differential pulldown, direct
  pulldown, and isoform lookup resolve identifiers by default and retain
  row-level audits. Set `resolve_duplicates = FALSE` for the exact identifier
  behaviour used by 0.4.0.
- Detect pseudo-switches created when oppositely directed submitted isoforms
  resolve to the same protein representative. Differential analysis drops the
  affected gene before graph expansion; same-direction collisions collapse to
  one representative.
- Distinguish cleaned submitted counts from resolved-representative counts.
  Dominance scores now explicitly use the resolved analysis table.

## Count aggregation and validation

- Add `aggregate_duplicate_counts()` for matrices, data frames, and optional
  single-assay `SummarizedExperiment` inputs. It preserves nonmapped features,
  feature order, gene and library totals, and emits a complete arithmetic and
  mapping audit.
- Expand coverage reports to show both resolution paths and network outcomes,
  with submitted and unique-representative totals reported separately.
- Add installed-resource, custom-graph, collision, audit-persistence, count
  arithmetic, overflow, and optional-container regression tests. The embedded
  `R/sysdata.rda` network is unchanged.

# IsoNetAnalyzeR 0.4.0

## Usability and validation

- Add a persistent three-state network coverage summary: unmapped, mapped
  without a predicted edge, and mapped with a predicted edge, overall and per
  supplied direction. Mapping below 10% now emits diagnostic warnings.
- Validate non-finite and fractional integer arguments, stored cutoffs, input
  duplicates and conflicts, transcript-to-gene mappings, and the complete
  custom-network contract.
- Support networks with no predicted edges and report complete per-gene
  processing failure clearly while preserving valid unmapped, isolated, and
  one-direction results.
- Correct partially missing multi-seed lookups.
- Return plot objects visibly without explicit printing, enabling assignment
  and headless `ggsave()` without creating `Rplots.pdf`.

## Documentation

- Define the package as analysis-agnostic and keep seed derivation and
  publication-specific analyses outside the package.
- Replace proxy-direction examples with a caller-prepared synthetic workflow,
  qualify compositional dIF guidance, and describe results as static predicted
  seed neighbourhoods.
- Document reproducible enrichment using seeded runs and frozen pathway lists.

# IsoNetAnalyzeR 0.3.2

## Correctness fixes

- Preserve mapped isolated seed isoforms and genes through cross-gene graph
  construction, component analysis, component plotting, single-direction
  plotting, and differential overlay plotting. A mapped seed with no retained
  interactions is now represented as a valid one-node/zero-edge result.
- Correct gene-level component bookkeeping so ENSG vertices are compared with
  gene identifiers rather than ENST seed identifiers when counting direction
  contributions.
- Plot a stored component using its authoritative node-table membership, so
  component IDs remain correct when component filters were used.
- Store explicit per-node seed flags in virtual-pulldown statistics, retaining
  a legacy fallback for serialized objects created before version 0.3.2.
- Give singleton summaries finite structural conventions (density, average
  degree, and clustering equal to zero; centralization remains undefined).

## Documentation and maintenance

- Correct package author/maintainer metadata, Rd markup, and tidy-evaluation
  global-variable declarations.
- Transfer release maintenance to Dimitrios S. Kanakoglou while retaining
  Nikolas Mumm as the original package author.
- Correct the GitHub installation URL, vignette name, and cross-gene summary
  column count in the README and vignette.
- Add regression and end-to-end tests covering gene/isoform component counts,
  isolated seeds, explicit filters, plotting, serialization, invalid public
  inputs, fully unmapped analyses, typed empty outputs, and single-direction
  sparse analyses.
- Run plotting regressions on an explicit temporary graphics device so package
  tests cannot leave `Rplots.pdf` in the source tree.
- Replace soft-deprecated igraph plotting aliases with their current API names;
  singleton and filtered plotting now complete without deprecation warnings and
  retain identical graph semantics.
