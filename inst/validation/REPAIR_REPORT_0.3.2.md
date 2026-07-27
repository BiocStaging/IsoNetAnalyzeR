# IsoNetAnalyzeR 0.3.2 repair and validation report

Prepared: 2026-07-27

## Executive summary

IsoNetAnalyzeR 0.3.1 contained two high-impact reporting/completeness bugs,
three additional correctness defects in downstream component and node
statistics, and several packaging, documentation, and maintenance problems.
The repaired 0.3.2 branch preserves the core virtual-pulldown interaction
graph construction while correcting how seeds, components, component
identities, and singleton statistics are carried into summaries and plots.

The most important scientific boundary is:

- the ENSG/ENST defect changed gene-level component summary fields, not the
  underlying isoform interaction edges;
- the isolated-seed defect removed valid mapped loss-of-interaction results
  from downstream cross-gene components and plots;
- non-singleton component node/edge inventories remain unchanged in
  repair-only real-data validation;
- historical gene-component summaries and component figures must not be used
  without regeneration.

The original upstream checkout is retained separately. The repair branch is
based on commit `1385b25` and the package version is 0.3.2. The embedded
network object is unchanged by this repair and has SHA-256
`3130209f090b9b83bbf2ebf328f29a0ad9a8fcd89922a277d0966325bc163ebb`.

## Review scope

The review covered:

- S4 object validity and compatibility;
- virtual pulldown and differential pulldown flow;
- isoform- and gene-level cross-gene graph construction;
- connected-component summaries and filtering;
- seed status and per-node statistics;
- singleton and sparse-network behavior;
- component, lookup, overlay, and filtered-network plotting;
- package metadata, generated documentation, examples, and vignette behavior;
- package tests, source build, installation, and R package checks;
- repair-only validation using frozen real-data pulldown sessions.

pairedGSEA and the study-specific analysis pipeline are outside the
IsoNetAnalyzeR package repair. Their separately identified issues are recorded
in the analysis repository and are not presented as package defects here.

## Repair matrix

| ID | Severity | Previous behavior | Correction | Scientific consequence |
|---|---|---|---|---|
| ISO-BUG-001 | High | Gene-level component vertices use ENSG identifiers, but `n_genes_up` and `n_genes_down` were counted against ENST seed vectors. | Component-table construction now receives the analysis level. Gene-level counts use the ENSG seed direction keyed by `gene_id`; isoform-level counts continue to use ENST membership. | Gene-level component counts and `mixed` classifications are corrected. Core isoform edges are unchanged. |
| ISO-BUG-002 | High | Mapped degree-zero seeds survived `virtual_pulldown()` but were deleted during cross-gene graph construction, component analysis, and plotting. | Unconditional isolate deletion was removed. Plot cleanup now retains isolated seed vertices while dropping only isolated non-seeds. `min_edges` remains an explicit opt-in filter. | Valid mapped loss-of-interaction observations and singleton components are restored. Historical component inventories and plots are incomplete. |
| ISO-BUG-003 | Medium | `plot_component()` reconstructed components from the unfiltered graph, so a stored component ID could identify different nodes after filters were applied. | The plot is induced from the authoritative stored component node table. | Filtered component IDs now refer to the exact component selected during analysis. |
| ISO-BUG-004 | Medium | `get_isoform_stats()` inferred seed status from mean distance equal to zero. Disconnected multi-seed graphs can give seed vertices infinite mean distance. | `analyze_subgraph()` records `is_seed` explicitly from order-zero vertices. Accessors use the stored flag and retain a legacy fallback for old serialized objects. | Seed status is correct for disconnected and multi-seed results without breaking older RDS objects. |
| ISO-BUG-005 | Low | Singleton density, clustering, degree, or centralization calculations could return `NaN`. | Density, mean degree, and clustering use zero for singletons; centralization remains undefined and is represented by `NA`. | Singleton tables are machine-readable and no longer mix valid zeros with `NaN`. |
| ISO-DOC-001 | Low | README installation text, vignette naming, and the documented cross-gene column count were stale or incorrect. | Installation text, vignette calls, example workflow, and the 32-column description were aligned with the package. | Users receive a coherent executable workflow description. |
| ISO-DOC-002 | Low | Generated class documentation contained corrupted intersection/union text. | Roxygen source uses portable wording and the generated Rd file was repaired. | Installed help renders correctly. |
| ISO-MAINT-001 | Low | Metadata did not parse robustly for all build/check tools, the package stayed at 0.3.1 after repairs, and testthat was undeclared. | Metadata was completed, version advanced to 0.3.2, testthat edition 3 was declared, and maintenance ownership moved to Dimitrios S. Kanakoglou while Nikolas Mumm remains an author. | The repaired package builds and its tests are part of package checking. |
| ISO-MAINT-002 | Low | Generated tidy-evaluation columns in network-diagnostic plots caused code-analysis notes. | `log10_x` and `n` were registered as global variables. | Package checks no longer report these false-positive notes. |
| ISO-MAINT-003 | Low | Singleton and filtered plotting paths used igraph aliases soft-deprecated in igraph 2.x. | Exact current equivalents `subgraph_from_edges()` and `delete_vertices()` are used. | Deprecation warnings are removed without changing graph contents. |

## Detailed changes

### ISO-BUG-001: gene-level ENSG/ENST bookkeeping

Affected code:

- `run_component_analysis()`
- `.build_component_table()`

At isoform level, component vertex names and the stored `up_isoforms` and
`down_isoforms` vectors are all ENST identifiers. At gene level, however,
component vertex names are ENSG identifiers. The previous implementation used
the ENST membership test for both levels. A valid gene-level mixed component
could therefore contain both up and down seeds while reporting zero
contributing genes in both directions.

The repaired function receives `level`. At gene level it indexes
`seed_direction` by component vertex name and tests the current
`GeneDiffPulldown@gene_id`. At isoform level the original ENST membership logic
is retained.

Regression tests construct equivalent ENSG- and ENST-level components and
verify `n_genes_up`, `n_genes_down`, and `mixed` independently.

### ISO-BUG-002: isolated mapped seeds

Affected code:

- `.build_cross_gene_graph()`
- `run_component_analysis()`
- `plot_component()`
- `.build_plot_graph()`
- `.plot_gene_overlay()`

An isoform can be present in the network but have no edge with
`interact == 1`. That is different from an unmapped isoform. The former is a
valid mapped singleton and represents absence of predicted interactions under
the chosen threshold; the latter cannot be evaluated.

The old downstream code applied unconditional degree-zero deletion after the
pulldown step. The repair preserves mapped seed isolates across isoform and
gene contraction, component decomposition, overlays, lookups, and plots.
Users can still request `min_edges = 1` when they intentionally want to omit
singletons.

Tests distinguish:

- interacting mapped seeds;
- mapped isolated seeds;
- completely unmapped seeds;
- a seed present in both directions;
- singleton removal through the explicit `min_edges` filter.

### ISO-BUG-003: filtered component identity

Affected code:

- `plot_component()`

Component IDs are assigned only after component filtering and size/edge
sorting. Re-decomposing the original unfiltered graph in the plotting function
could shift those IDs. The repaired plot path reads the stored node table for
the requested component and induces that exact graph. A regression filters
out the singleton component, plots component 1, and verifies that only the
stored connected nodes are present.

### ISO-BUG-004: explicit seed flags

Affected code:

- `analyze_subgraph()`
- `get_isoform_stats()`
- `.safe_analyze()`

Distance from multiple seeds is not a reliable seed predicate in a
disconnected graph. The repaired node-stat schema contains a logical
`is_seed` column. The empty schema contains the same field so sparse and
all-failure results can be combined safely. When an older serialized result
lacks the column, the accessor uses the historical fallback rather than
invalidating the object.

Tests cover disconnected seeds, typed empty results, serialization, reload,
and a manually constructed legacy node table.

### ISO-BUG-005: singleton metric conventions

Affected code:

- component structural summaries;
- `GeneDiffPulldown` graph summaries;
- virtual-pulldown node summaries.

For a one-node graph:

- density is `0`;
- average degree is `0`;
- clustering is `0`;
- degree centralization is `NA` because it is undefined.

These conventions distinguish absence of edges from a failed calculation.

### Plotting and graphics-device behavior

Package plotting tests previously allowed plot side effects to open the
default graphics device and leave `Rplots.pdf` in the source tree. Plot tests
now run inside an explicit temporary PDF device, close and remove that device
on exit, and assert that no default `Rplots.pdf` appears. This is a test
harness correction; user-facing plot return values remain ggplot objects.

## Public interface and compatibility

No exported function or S4 class was removed.

Observable changes are:

- per-isoform statistics include logical `is_seed`;
- gene-level component summaries contain corrected direction counts;
- mapped singleton components appear by default;
- filtered component plots use the stored component membership;
- singleton numeric conventions are explicit.

Old serialized objects without `is_seed` remain readable. Users who require
the former edge-containing-only component inventory can request
`min_edges = 1`.

## Real-data validation

Repair-only validation reused frozen differential pulldown graphs so that the
component repairs could be isolated from differential modeling and network
expansion.

Across the historical GSE141432 validations:

| Contrast | Old components | Repaired components | Restored singletons | Non-singleton inventory unchanged |
|---|---:|---:|---:|---|
| Lean vs Obese NGT | 2 | 6 | 4 | Yes |
| Lean vs Obese T2D | 31 | 42 | 11 | Yes |
| Obese NGT vs Obese T2D | 1 | 1 | 0 | Yes |

Across GSE141432 and the historical GSE231509 validation set, 24 omitted
gene-level singleton components were restored. Every previously present
non-singleton component retained the same node/edge-count inventory.
GSE231509 is cited only as package-validation evidence and is not part of the
new biological analysis.

The authoritative July workflow also demonstrates that the repaired primary
GSE141432 comparison (`Obese_T2D_minus_Obese_NGT`) changes from no strict
post-mapping bidirectional seed to one retained gene with two mapped isolated
seed rows after removal of an undocumented transcript-effect filter. That
pipeline change is separate from the IsoNetAnalyzeR package repairs.

## Package-check evidence

The authoritative pre-handoff standard check is:

- package: IsoNetAnalyzeR 0.3.2;
- tarball SHA-256:
  `33bb03050eb94461bcb89bd77345e60a511d228fc18069020c6e8381f1da5c1a`;
- log:
  `final_pipeline/artifacts/checks_no_manual/IsoNetAnalyzeR.Rcheck/00check.log`;
- result: `Status: OK`.

That check predates the maintainer handoff and expanded input/sparse tests.
The final release evidence must therefore be regenerated from the committed
branch. The prior `--as-cran` failure was produced without LaTeX, qpdf, recent
HTML Tidy, or the R V8 package. The release check is run only after installing
those tools and must successfully build the PDF manual before any remaining
finding is classified.

## Affected historical outputs

The following output classes require regeneration or explicit
superseded status:

- gene-level component overviews containing `n_genes_up`,
  `n_genes_down`, or `mixed`;
- component counts or figures produced before isolated-seed retention;
- filtered component figures produced before stored-membership plotting;
- per-node seed summaries inferred only from mean distance;
- any downstream ranking or narrative based on those fields;
- pairedGSEA paired-ORA outputs affected by the separate universe/join bug.

The detailed path-level decision is provided in the companion
`affected_outputs.csv`. Historical data are retained for provenance and are
not overwritten.

## Remaining limitations

- Network coverage depends on exact transcript identifiers in the embedded
  interaction network. Unmapped transcripts must always be reported.
- Gene-level adjacency uses an any/max collapse: a gene pair is connected if
  at least one contributing isoform pair passes the interaction threshold.
- A singleton is valid evidence of mapped interaction loss, but is not
  sufficient by itself for a comparative rewiring claim.
- Topology classifications are heuristics and require biological review.
- Tests concentrate on the repaired and high-risk paths. Public enrichment and
  diagnostic helpers outside this repair boundary receive additional coverage
  only when it does not delay the release gate.

## Release acceptance

The final 0.3.2 release candidate is accepted only when:

- the committed test suite has no failures or warnings;
- no plotting test creates `Rplots.pdf`;
- source build and clean-library installation succeed;
- `R CMD check --no-manual` returns `Status: OK`;
- `R CMD check --as-cran` successfully builds and validates the PDF manual;
- real-data repair invariants remain satisfied;
- the source commit, tarball, logs, reports, and transfer artifacts have
  recorded SHA-256 values.
