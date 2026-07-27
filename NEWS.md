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
