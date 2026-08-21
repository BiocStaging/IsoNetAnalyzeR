# IsoNetAnalyzeR issue register

This register covers problems found during the 2026-07 audit. IDs are stable
and are referenced by the change log and validation report.

| ID | Class | Severity | Status | Problem and impact |
|---|---|---|---|---|
| ISO-BUG-001 | bug | high | fixed | Gene-level component summaries compared ENSG vertex names with ENST seed lists. `n_genes_up`, `n_genes_down`, and `mixed` could be wrong at gene level. Core isoform pulldowns and graph edges were unaffected. |
| ISO-BUG-002 | bug | high | fixed | Mapped degree-zero seeds were preserved by `virtual_pulldown()` but deleted by cross-gene construction, component analysis, component plotting, generic plotting, and overlay cleanup. Valid single-node inventory results disappeared. |
| ISO-BUG-003 | bug | medium | fixed | `plot_component()` rebuilt and sorted unfiltered components rather than using the stored filtered membership. Component IDs could refer to a different graph when filters had been used. |
| ISO-BUG-004 | bug | medium | fixed | `get_isoform_stats()` inferred seed status from an average distance equal to zero. Multiple disconnected seeds can have infinite average distance even though each is an order-zero seed. |
| ISO-BUG-005 | bug | low | fixed | Singleton structural summaries could contain `NaN` density/clustering/centralization values. Density, mean degree, and clustering now use zero; centralization remains undefined. |
| ISO-DOC-001 | documentation defect | low | fixed | README install URL contained a trailing space, named a nonexistent vignette, and the README/example workflow described 35 rather than 32 cross-gene columns. |
| ISO-DOC-002 | documentation defect | low | fixed | Generated class Rd contained corrupted intersection/union markup; source documentation also used symbols that did not render reliably. |
| ISO-MAINT-001 | maintenance defect | low | fixed | Author/maintainer metadata was incomplete for robust build/check tooling and testthat was undeclared. |
| ISO-MAINT-002 | maintenance defect | low | fixed | Tidy-evaluation variables in network diagnostic plots were not registered, causing check notes in some environments. |
| ISO-MAINT-003 | maintenance defect | low | fixed | Singleton and filtered plot paths used igraph aliases soft-deprecated in igraph 2.x, producing four package-test warnings under igraph 2.3.0. Current API names are now used with identical arguments and graph semantics. |
| ISO-METH-001 | methodological limitation | informational | documented | Gene-level adjacency is an any/max-threshold collapse: a gene-gene edge exists when at least one contributing isoform pair passes the interaction threshold. It does not represent an average isoform interaction. |
| ISO-METH-002 | methodological limitation | medium | managed | Network interpretation depends on transcript-ID coverage. Versioned and protein-identical identifiers now resolve only when compatible with the active graph; unresolved inputs and post-mapping one-direction-only genes remain explicitly reported and are not package failures. |
| ISO-METH-003 | methodological limitation | medium | documented | Singletons and tiny components are valid static-network inventory results. They do not establish biological interaction loss, and mixed components do not establish rewiring. |
| ISO-METH-004 | methodological limitation | high | managed | Count aggregation can recover redundant seeds but can also make oppositely directed protein-identical transcripts collide. Version 0.5.0 drops the complete affected gene before graph analysis and records every submitted identifier, direction, and representative. Historical 0.4.0 retained seed tables were collision-free; this guard is prospective. |
| ISO-SPARSE-001 | expected sparse-data behavior | informational | documented | A direction with no mapped seed produces no subgraph. A mapped seed with no predicted interaction produces a one-node/zero-edge subgraph. These cases are now distinguished. |

## Reproduction summaries

ISO-BUG-001 is reproduced by a gene-level component whose vertex is an ENSG
identifier while its originating up/down seed vectors contain ENST identifiers.
The old comparison returned zero contributing genes.

ISO-BUG-002 is reproduced by a network vertex present in the network but with
no `interact == 1` edge. `virtual_pulldown()` returns that seed as a one-node
graph; downstream calls formerly removed it via unconditional
`delete.vertices(degree == 0)`.

ISO-BUG-003 is reproduced by running component analysis with `min_edges = 1`
and then plotting a stored component ID. The prior plot path decomposed the
unfiltered graph, allowing omitted singleton components to change the index.

ISO-BUG-004 is reproduced with two disconnected order-zero seed vertices.
Distances from each seed to the other are infinite, so column means are
infinite even though both vertices are explicit seeds.

ISO-MAINT-003 is reproduced under igraph 2.3.0 by the singleton
single-direction, overlay, lookup, and post-serialization lookup plot paths.
The former `delete.vertices()` and `subgraph.edges()` names are exact
soft-deprecation wrappers around `delete_vertices()` and
`subgraph_from_edges()`; replacing the aliases is output-neutral.

## Semantics of gene-level networks

Isoform graphs are contracted by parent gene and simplified. Consequently,
gene-level adjacency means that at least one isoform-level edge contributed to
that pair (equivalently, any/max thresholding across contributing isoform
pairs). ISO-BUG-001 was downstream summary bookkeeping; it did not alter this
adjacency rule or create/remove core isoform edges.
