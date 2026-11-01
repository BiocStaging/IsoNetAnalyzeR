# IsoNetAnalyzeR 0.99.0

Version prepared for submission to Bioconductor.

## Bundled network

- The bundled network now contains the 158,458 predicted interaction edges
  (`interact == 1`) and all 58,164 transcripts, instead of all 2,853,266
  scored pairs. Every analysis already restricted the network to these
  edges, so pulldowns, lookups, components, enrichment and diagnostics are
  unchanged; this was verified against 0.6.5. Pulldowns are about ten times
  faster and the source package is under Bioconductor's 5 MB limit.
- Counts of all scored edges change accordingly: the session edge count, the
  eligibility audit fields `edges_before`/`edges_after`, and the diagnostics
  field `full_edges` now count predicted interaction edges.
  `plot_network(obj, all_edges = TRUE)` on a whole session no longer draws
  non-interacting pairs.
- `inst/scripts/reduce_network.R` derives the new network from the 0.6.5
  file; `inst/PROVENANCE.md` documents both checksums.

## New features

- `get_metadata()` returns session metadata such as coverage, seed-support
  and eligibility audits, replacing direct slot access.
- `get_pulldown_graph()` returns a gene's up- or down-seeded pulldown graph.
- `get_network()` returns the session network.
- Network plots label nodes with gene symbols. `plot_network()` and
  `plot_component()` gain `node_labels = c("symbol", "id", "none")` to switch
  to transcript or gene identifiers or to hide labels, and gene plot titles
  show the symbol next to the gene ID. Gene labels in the cross-gene scatter
  and boxplots also use symbols.
- `isonet_gene_names()` returns the GENCODE v41 gene symbols of the bundled
  network (`inst/extdata/gene_names_gencode_v41.csv.gz`, produced by
  `inst/scripts/make_gene_names.R`), for annotating result tables.

## Printed output

- Printed summaries keep their column alignment, which `cli` had collapsed;
  the coverage report labels share one width.
- Session creation reports counts with thousands separators and states that
  excluded NMD and NSD transcripts take their edges with them, which explains
  the two edge counts. The session summary calls the bundled network
  "bundled".
- "Access with" hints use the name of the variable passed as `obj`, and the
  enrichment summaries show "user-supplied" instead of `NA` for gene sets
  passed through `pathways`.

## Plot fixes

- The `direction = "overlay"` plot of `plot_network()` follows the style of
  the paper figures: dashed down-only edges, thicker shared edges, a short
  centred title with the gene symbol, and legends below the plot. The
  vignette draws it at gene level, as the paper does.

- `plot_component(colour_by = "direction")` draws one legend for colour and
  one for seed status instead of two legends both titled "Direction";
  equidistant neighbours are labelled "Equidistant".
- The `plot_cross_gene_overview()` subtitle no longer describes the up/down
  comparison as a change between conditions.
- The gene-connectivity panel of `plot_network_diagnostics()` leaves room for
  the count above its tallest bar.

## Documentation and data

- New vignette, "Introduction to IsoNetAnalyzeR", analysing isoform switches
  between obese type 2 diabetic and lean adipose tissue (GEO GSE141432).
  It credits Isabella Østerlund's network model and protein-identical map
  and has an Acknowledgements section.
- Example data: `gse141432_t2d_vs_lean_seeds.csv` (137 seeds of 65 genes)
  and the MSigDB Hallmark gene sets (`msigdb_hallmark_ensembl.rds`), with
  generating scripts in `inst/scripts/`.
- `inst/PROVENANCE.md` describes how the network was derived from STRING
  physical interactions and attributes STRING's CC BY 4.0 licence.
- The README shows the CLTCL1 example (`inst/scripts/make_readme_figure.R`)
  and explains how to install from GitHub until the package is on
  Bioconductor.
- Every exported function and class has runnable examples.
- The help page for `run_network_diagnostics()` no longer truncates the
  description of the default sticky-node rule.

## Package infrastructure

- `biocViews`, `URL` and `BugReports` added; ORCID iDs, copyright holders
  and the thesis supervisors listed in `Authors@R`; licence declared as
  `MIT + file LICENSE` with the full text in `LICENSE.md`.
- Requires R 4.6.0 or later, the R version of Bioconductor devel.
- `msigdbr` moved to Suggests; `withr` added to Imports to draw the
  diagnostics plot subsample with a local seed.
- Condition messages no longer build text with `paste()`; message text is
  unchanged.

# IsoNetAnalyzeR 0.6.5

- Add a concise README with direct user-guide links, a self-contained installed instructions file, and a separate network-eligibility explanation with references.
- Keep the complete MIT licence in a single `LICENSE` file; update R licence metadata and the README link while preserving the copyright notice and terms.
- Consolidate the default NMD/NSD exclusions from 0.5.5 and the optional retained-significant-seed validation from 0.6.0 into the 0.6.5 release. These scientific rules and the pinned reference resources are unchanged from the checked 0.6.0 analysis archive.
- Separate focal-seed role from neighbourhood membership in differential overlay plots (ISO-PLOT-001). Focal seeds now use a reddish-purple fill (`#CC79A7`) with directional seed shapes; non-seed nodes and edges use green/up, orange/down and blue/both-direction colours. This corrects the plotting behaviour that was not included in the earlier 0.6.0 archive.
- Add an executable retained-support example to the vignette and full workflow, clarify nonsignificant complementary seeds and the exploratory three-column interface, and update installation and validation instructions for 0.6.5.
- Preserve the raw embedded network, protein-identical transcript map and GENCODE v41 biotype table. Existing checksum-pinned study archives and RNA results remain separate from this release.

# IsoNetAnalyzeR 0.6.0

- Add a retained-seed support gate to `run_differential_pulldown()` via `seed_padj_col` and `seed_fdr`. After decay filtering, ID resolution and mapping, every retained gene must contain at least one significant retained seed and at least one up and one down seed. Non-significant complementary seeds remain allowed; there is no seed maximum.
- Store gene/row support audits in `metadata$analysis_seed_support`. An excluded, unmapped or collapsed-away transcript cannot supply significance for the final seed set.
- RNA testing and effect-size selection remain upstream. The GSE141432 v0.6 analysis enables the gate explicitly; three-column exploratory inputs remain compatible when `seed_padj_col = NULL`.

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
