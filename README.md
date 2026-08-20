# IsoNetAnalyzeR

IsoNetAnalyzeR builds and compares predicted interaction neighbourhoods around caller-supplied isoform seed sets. The package is analysis-agnostic: it does not test differential transcript usage, select seeds, infer biological states, or claim condition-specific protein interactions.

## Installation

Release 0.5.0 can be installed from its pinned tag:

```r
remotes::install_github("dskanakoglou/IsoNetAnalyzeR@v0.5.0",
                        build_vignettes = TRUE)
```

The bundled static network contains 58,164 transcript vertices and 2,853,266 total edges, of which 158,458 have `interact == 1`. A custom undirected `igraph` can be supplied instead; vertices require unique `name` values and a non-missing `gene` attribute, and edges require a binary `interact` attribute.

## Input contract

`run_differential_pulldown()` accepts a CSV, Excel file, or data frame with three required columns:

| Column | Meaning |
|---|---|
| `gene_id` | Gene assigned to the submitted transcript |
| `isoform_id` | Submitted transcript ID; version suffixes and packaged protein-identical duplicates are resolved when graph-compatible |
| `direction` | Caller-defined seed set: `"up"` or `"down"` |

The names `up` and `down` are labels for two directional seed sets. They do not turn the static network into measured condition-specific interaction networks.

```r
library(IsoNetAnalyzeR)

seeds <- data.frame(
  gene_id = c("ENSG00000000419", "ENSG00000000419"),
  isoform_id = c("ENST00000371588", "ENST00000371582"),
  direction = c("up", "down")
)

obj <- IsoNetAnalyzeR("example")
obj <- run_differential_pulldown(
  obj, "seed_neighbourhoods", seeds,
  order = 1, confirm = FALSE
)

summary <- get_cross_gene_summary(obj, "seed_neighbourhoods")
gene_result <- get_gene_result(
  obj, "seed_neighbourhoods", "ENSG00000000419"
)
p <- plot_network(gene_result)
```

By default, the package normalizes versioned ENST identifiers and remaps a
protein-identical redundant transcript only when the representative is a valid
vertex in the active graph, the redundant transcript is absent, and the gene
annotation agrees. Exact vertices always win. Incompatible mappings remain
untouched and are included in the audit. Pass `resolve_duplicates = FALSE` to
retain the exact-ID behaviour of 0.4.0.

The pre-run coverage report separates submitted IDs by resolution path
(direct, version-normalized, duplicate-remapped, or unresolved) and network
outcome:

- unmapped;
- mapped with no predicted edge in the pinned network;
- mapped with at least one predicted edge.

Submitted-ID and unique-representative totals are reported separately. Coverage
is retained under `obj@metadata$analysis_network_coverage[[label]]`, while the
cleaned submitted table and row-level resolution audit are retained under
`obj@metadata$analysis_duplicate_resolution[[label]]`. A mapped isolated
transcript is a valid single-node result; it is not evidence that the transcript
has no biological interaction.

If same-direction submitted isoforms resolve to one representative, they are
collapsed before analysis. If opposite directions collide, the complete gene
is excluded with a warning to prevent a protein-identical pseudo-switch.

## Aggregating counts before statistical fitting

`aggregate_duplicate_counts()` applies the full transcriptome-wide map before
differential modelling. It accepts a count matrix or data frame plus optional
feature metadata, or an optional single-assay `SummarizedExperiment`:

```r
aggregated <- aggregate_duplicate_counts(counts, feature_data)
counts_deduplicated <- aggregated$object
feature_data_deduplicated <- aggregated$feature_data
aggregated$audit
aggregated$summary
```

Per-gene and global count totals are preserved, but representative numerators
and the fitted feature set change. A complete statistical refit is therefore
required; adjusted values from the unaggregated fit must not be reused.

## Choosing seed isoforms

Seed selection belongs upstream of IsoNetAnalyzeR. The package preserves the supplied table after validation and does not implement a particular statistical workflow.

Isoform fractions are compositional. When all isoforms of a gene use a common denominator in the same samples, their changes sum approximately to zero, so opposing dIF signs are coupled rather than independent confirmation. A gene-level differential-usage screen followed by a pre-specified within-gene or stagewise driver rule is one defensible design. Any minimum absolute effect floor should be chosen for the study and declared before network analysis; it is not a universal package threshold. Do not add a non-significant opposite-direction transcript merely to force two seed sets.

One-direction genes are valid inputs. Unmapped genes are retained as typed empty results rather than silently removed.

## Interpretation

The package expands each seed set through `interact == 1` edges and applies the stored or supplied internal-edge cutoff. Comparisons describe two seed-derived neighbourhoods in the same static predicted network:

- `n_edges_gained` is retained for compatibility and means edges unique to the `up`-seed neighbourhood;
- `n_edges_lost` means edges unique to the `down`-seed neighbourhood;
- `mixed` means a component contains seeds from both labels;
- topology classes such as `hub-spoke` are descriptive heuristics, not evidence of a master regulator or disease mechanism.

Use `run_component_analysis()` to summarize merged components and `lookup_virtual_pulldown()` for an ad hoc seed neighbourhood. `run_network_diagnostics()` describes coverage, isolation, components, and high-degree nodes in the loaded network.

## Reproducible enrichment

For a reproducible enrichment analysis, save the exact named pathway list with the analysis, use stable gene identifiers, set the random seed, and record package versions and parameters. Fetching a live pathway collection by name can change results over time.

```r
pathways <- readRDS("frozen_pathways.rds")
set.seed(2026)
obj <- run_fgsea(
  obj, "seed_neighbourhoods", "frozen_gsea",
  pathways = pathways, rank_by = "delta_nodes"
)
```

The example workflow is in `inst/examples/full_workflow.R`, and the vignette provides a small fully synthetic run.
