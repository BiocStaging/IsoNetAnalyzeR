# IsoNetAnalyzeR

IsoNetAnalyzeR builds and compares predicted interaction neighbourhoods around caller-supplied isoform seed sets. The package is analysis-agnostic: it does not test differential transcript usage, select seeds, infer biological states, or claim condition-specific protein interactions.

## Installation

Install the 0.6.5 source archive into an R environment with the package
dependencies available:

```r
install.packages("IsoNetAnalyzeR_0.6.5.tar.gz", repos = NULL, type = "source")
packageVersion("IsoNetAnalyzeR")
```

After the `v0.6.5` tag has been published to the official repository, the
matching pinned GitHub installation is:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "IsoformAnalysisGroup/IsoNetAnalyzeR@v0.6.5",
  build_vignettes = TRUE,
  upgrade = "never"
)
```

The unfiltered bundled static network contains 58,164 transcript vertices and 2,853,266 total edges, of which 158,458 have `interact == 1`. Versions 0.5.5 and later exclude annotated NMD and NSD transcripts by default when creating a session. A custom undirected `igraph` can be supplied instead; vertices require unique `name` values and a non-missing `gene` attribute, and edges require a binary `interact` attribute. Custom vertices not covered by the pinned annotation also require a `transcript_biotype` attribute or an explicit annotation table when either exclusion is enabled.

## Protein-network eligibility: NMD and NSD

The default introduced in 0.5.5 is conservative for protein-interaction interpretation:

```r
obj <- IsoNetAnalyzeR("protein_network", include_nmd = FALSE, include_nsd = FALSE)
# Explicit opt-in for exploratory or legacy comparisons:
all_transcripts <- IsoNetAnalyzeR("unfiltered", include_nmd = TRUE, include_nsd = TRUE)
# The two flags are independent.
annotation <- isonet_transcript_biotypes()
```

`nonsense_mediated_decay` (NMD) and `non_stop_decay` (NSD) are **transcript** biotypes. A gene labelled `protein_coding`, an annotated ORF, or a predicted interaction does not establish a stable functional protein product for every transcript. GENCODE associates these decay classes with premature termination or absence of a stop codon, respectively. Excluding them by default reduces the risk of interpreting decay-associated RNA changes as changes between stable protein interaction partners. This is an annotation-based eligibility policy, not proof that these RNAs are never translated. NMD is translation-dependent, and protein products can be made and subsequently destabilized. [GENCODE biotype definitions](https://www.gencodegenes.org/pages/biotypes.html), [experimental evidence](https://pmc.ncbi.nlm.nih.gov/articles/PMC8694075/).

The session constructor removes excluded vertices and their incident edges **before** neighbourhood expansion and internal-edge fraction filtering. Thus the same policy covers seeds and neighbours. Submitted decay-associated aliases are blocked before protein-identical remapping; they cannot bypass exclusion through a coding representative. The vertex audit is in `obj@metadata$network_biotype_filter`; excluded submitted seed rows are in `obj@metadata$analysis_biotype_exclusions[[label]]`. Unknown transcript biotypes produce an actionable error if either exclusion is enabled. An explicit table overrides vertex annotations, which override the packaged GENCODE v41 mapping; gene-level biotypes are never substituted.

These options do not filter RNA count matrices, recompute isoform fractions, or change pairedGSEA tests. Keep decay-associated transcripts in the RNA analysis when the question concerns total transcript usage, then apply protein eligibility before network seed selection. Removing them from the RNA denominator would change that question. Recheck bidirectionality after eligibility and mapping; do not silently fill a missing direction with an ineligible transcript.

All remaining vertices are only eligible predictions, not experimentally validated proteins. Other biotypes require their own interpretation. Removing vertices also changes the denominator of the network's internal-edge fraction, so a previously excluded coding neighbour may cross the same cutoff. The inherited network cutoff was calibrated on the original network; reusing it in a filtered-network sensitivity run does not establish that it remains optimal.

Existing saved sessions and their cached results are not retroactively filtered. Construct a new session from the original network and rerun analyses when changing policy. Do not try to restore excluded vertices merely by changing metadata on an already filtered graph. Historical versions through 0.5.0 retain their unfiltered defaults; the 0.6.5 release includes the exclusions introduced in 0.5.5.

## Input contract

### Retained significant support

Version 0.5.5 introduced the NMD/NSD exclusions. The post-filter support check
introduced in 0.6.0 is included in 0.6.5 for statistically selected seed sets:

```r
# seeds already passed the study's gene screen and transcript effect floor.
obj <- run_differential_pulldown(
  obj, "T2D_minus_Lean", seeds, order = 1,
  seed_padj_col = "within_gene_padj", seed_fdr = 0.05,
  confirm = FALSE
)
```

This mode requires **at least one significant retained mapped seed, at least
one up seed, and at least one down seed per gene**. It is checked after decay
exclusions and identifier resolution. Unmapped rows are removed, and a
significant transcript that was excluded or collapsed away cannot support
the surviving pair. A non-significant opposite-direction partner is allowed;
all qualifying seeds are used, with no maximum. Mapped isolated seeds are
allowed: statistical support does not require an observed or predicted edge.
Gene/row audits are recorded in `obj@metadata$analysis_seed_support[[label]]`.

The package does not calculate or readjust these p-values. The caller must
declare the statistical method and apply its gene/effect filters before this
check. The GSE141432 0.6.0 pipeline enables both its retained-seed selector
and this package guard. The generic three-column interface remains available
with `seed_padj_col = NULL` and makes no statistical-support claim.

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

In differential overlay plots, focal seeds use a reddish-purple fill
(`#CC79A7`) and directional seed shapes. Non-seed nodes and edges use green
for the up-seed neighbourhood, orange for the down-seed neighbourhood, and
blue when present in both. This keeps seed role visually distinct from shared
neighbourhood membership.

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

Isoform fractions are compositional. When all isoforms of a gene use a common denominator in the same samples, their changes sum approximately to zero, so opposing dIF signs are coupled rather than independent confirmation. A gene-level differential-usage screen followed by a pre-specified within-gene or stagewise driver rule is one defensible design. Any minimum absolute effect floor should be chosen for the study and declared before network analysis; it is not a universal package threshold. A non-significant complementary transcript may be included under a declared gene-screen/driver/effect-size pairing rule, but its selection is not independent statistical confirmation. Do not describe such a partner as individually significant. All qualifying transcripts can be used jointly within each direction; the package imposes no maximum seed count. A study comparing opposing usage changes should require at least one eligible up and one eligible down seed after annotation filtering and network mapping.

This joint-seed design is also described in Nikolas A. Mumm, *An R Package for Isoform-Specific Interaction Rewiring and Disease Mechanism Discovery* (DTU master's thesis, July 2026), sections 2.3.3–2.3.4, pp. 12–13. The thesis used a significant mapped transcript plus an opposite-sign log2FC proxy where needed, with input selection implemented in project analysis code. A study using signed dIF and a gene-level driver rule is a different selection protocol and should not claim exact reproduction of that method. In particular, a significant driver elsewhere in the gene is not equivalent to a significant retained network seed.

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
