# IsoNetAnalyzeR audited change log

Release target: local development version 0.3.2. Original upstream clone
`../source` is not modified.

| Issue | File / symbol | Previous behavior | Change and rationale | Test / output consequence |
|---|---|---|---|---|
| ISO-BUG-001 | `R/component_analysis.R`, `.build_component_table()` | Used ENST membership checks for both isoform- and gene-level graphs. | Added explicit `level`; gene level uses the ENSG vertex direction keyed by `gr@gene_id`. | Existing two regression tests verify ENSG and ENST paths. Real gene component summaries must be regenerated. |
| ISO-BUG-002 | `R/graph_utils.R`, `.build_cross_gene_graph()` | Deleted every isolated vertex after merge/contraction. | Removed unconditional deletion; mapped singleton seeds remain graph vertices. | End-to-end test inventories connected plus singleton components. Repair-only real-data regeneration should restore 24 gene-level singleton components. |
| ISO-BUG-002 | `R/component_analysis.R`, `run_component_analysis()` | Deleted isolated vertices a second time. | Removed deletion; explicit `min_edges` remains available to exclude singletons intentionally. | Tests compare unfiltered and `min_edges = 1` results. |
| ISO-BUG-002/003 | `R/plot_component.R`, `plot_component()` | Deleted isolates and reconstructed component IDs from an unfiltered decomposition. | Retains isolates and induces the requested graph from the authoritative stored node table. | Singleton component returns a valid ggplot; filtered IDs remain exact. |
| ISO-BUG-002 | `R/plot_network.R`, `.build_plot_graph()` and callers | Plot cleanup removed all isolates. | Added seed-aware cleanup: retain isolated seeds and remove only isolated non-seeds. Callers pass actual order-zero seeds before contraction. | Lookup and single-direction singleton plots render. |
| ISO-BUG-002 | `R/plot_network.R`, `.plot_gene_overlay()` | Overlay cleanup removed all isolates and then failed on a seed-only graph. | Retain vertices whose `seed_direction` is not `non_seed`; omit an empty edge scale and revise the empty error. | Isolated up/down overlay renders without an irrelevant scale warning. |
| ISO-BUG-004 | `R/virtual_pulldown_analysis.R`, `analyze_subgraph()` | Returned distance but no explicit seed flag. | Adds `is_seed` directly from order-zero vertex membership. | Disconnected multi-seed test remains correctly flagged. |
| ISO-BUG-004 | `R/diff_pulldown_analysis_class.R`, `get_isoform_stats()` | Recomputed `is_seed` as `distance_from_seed == 0`. | Uses stored flag; legacy objects without the column keep the old fallback. | Serialization and legacy-schema compatibility are tested. |
| ISO-BUG-004 | `R/run_differential_pulldown.R`, `.safe_analyze()` | Empty schema lacked `is_seed`. | Empty results now match the new node-stat schema. | Rbind/accessor behavior remains stable. |
| ISO-BUG-005 | `R/component_analysis.R`, `R/diff_pulldown_gene_class.R`, `R/virtual_pulldown_analysis.R` | Singleton metrics could be `NaN`. | Defines density/degree/clustering as zero; centralization as `NA`. | End-to-end tests assert finite conventions. |
| ISO-MAINT-001 | `DESCRIPTION` | Incomplete explicit author metadata and no testthat config. | Completes metadata, declares testthat edition 3, bumps 0.3.1 to 0.3.2. | Package build/check includes tests. |
| ISO-MAINT-002 | `R/plot_network_diagnostics.R` | `log10_x` and `n` were absent from global-variable registration. | Registers both tidy-evaluation names. | Removes check notes without changing plots. |
| ISO-MAINT-003 | `R/plot_network.R`, `.build_plot_graph()` and `.plot_gene_overlay()` | Used the soft-deprecated igraph aliases `subgraph.edges()` and `delete.vertices()`, producing four warnings in singleton plot tests under igraph 2.3.0. | Uses `subgraph_from_edges()` and `delete_vertices()` with the same arguments. These are the current names for the same operations. | No-warning singleton and filtered-graph regressions verify unchanged vertices, edges, attributes, and isolated-seed retention. Scientific and visual outputs are unchanged. |
| ISO-DOC-002 | `R/diff_pulldown_analysis_class.R`, generated class Rd | Intersection/union symbols rendered incorrectly/corrupted. | Uses portable text markup and repairs Rd. | Documentation builds successfully. |
| ISO-BUG-002/004 | `man/analyze_subgraph.Rd`, `man/get_isoform_stats.Rd`, `man/plot_network.Rd` | Generated help omitted explicit seed flags and described all isolated plot nodes as removed. | Regenerated Rd from the corrected source documentation. | Installed help matches singleton and legacy behavior. |
| ISO-DOC-001 | `README.md`, `vignettes/vignette_introduction.Rmd`, `inst/examples/full_workflow.R` | URL, vignette call, and 35-column count were wrong or stale. | Removes URL space, names `vignette_introduction`, and states 32 columns consistently. | README, vignette, and full-workflow example now match package artifacts. |
| documentation | vignette chunk options | Install/example chunks could be tangled during checks. | Marks non-executed illustrative chunks `purl = FALSE`. | Vignette build is reproducible; full workflow validation lives in tests. |
| release | `NEWS.md`, `inst/validation/*` | No local release/audit record. | Adds versioned release notes, issue register, change log, and protocol. | Every edited symbol is traceable to an issue and test. |

No pairedGSEA file is part of this change log. The pairedGSEA correction lives
in the analysis workspace by explicit user instruction.
