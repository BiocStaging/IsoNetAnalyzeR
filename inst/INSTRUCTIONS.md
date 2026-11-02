# IsoNetAnalyzeR user instructions

These instructions cover a complete analysis of your own selected transcript isoforms. The [example script](examples/full_workflow.R) supplies synthetic data for trying the package first.

## Installation

IsoNetAnalyzeR requires R 4.6 or later. Until it is available from Bioconductor, install it from GitHub; `build_vignettes = TRUE` also builds the tutorial:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("IsoformAnalysisGroup/IsoNetAnalyzeR",
                     dependencies = TRUE, build_vignettes = TRUE)
library(IsoNetAnalyzeR)
```

Once it is available from Bioconductor, use `BiocManager::install("IsoNetAnalyzeR")`.

## Try the example

The installed tutorial walks through a complete analysis of the bundled GSE141432 example seeds:

```r
vignette("IsoNetAnalyzeR", package = "IsoNetAnalyzeR")
```

The installed example script creates a small synthetic network, analyses its seed sets, generates plots, and demonstrates component analysis and optional seed-support validation:

```r
library(IsoNetAnalyzeR)
source(system.file("examples", "full_workflow.R", package = "IsoNetAnalyzeR"))
```

## Prepare your input

Supply a data frame, CSV file, or Excel file with these columns:

| Column | Required values |
| --- | --- |
| `gene_id` | Gene identifier corresponding to each transcript; use Ensembl gene IDs with the bundled human network. |
| `isoform_id` | Transcript identifier; use Ensembl transcript IDs with the bundled network. |
| `direction` | Exactly `up` or `down`, identifying your two seed sets. |

Each row supplies a seed transcript. Multiple transcripts per gene and direction are allowed. Choose the seeds using your study's statistical and effect-size criteria before running the network analysis. The package does not fit those statistical models.

Version suffixes and compatible protein-identical transcript identifiers are resolved automatically. Exact graph vertices take precedence. Same-direction seeds that resolve to one representative are collapsed; a gene is excluded if its opposite directions resolve to the same representative. Set `resolve_duplicates = FALSE` when exact-ID matching is required.

## Run your own analysis

The following uses `seeds.csv`, containing your prepared table:

```r
library(IsoNetAnalyzeR)

seeds <- read.csv("seeds.csv", stringsAsFactors = FALSE)
obj <- IsoNetAnalyzeR("my_analysis")
obj <- run_differential_pulldown(
  obj,
  label = "comparison",
  input = seeds,
  order = 1,
  confirm = FALSE
)

results <- get_cross_gene_summary(obj, "comparison")
results
write.csv(results, "network_summary.csv", row.names = FALSE)
```

`order = 1` expands first-order seed neighbourhoods. If `cutoff` is omitted, the session's stored internal-edge cutoff is used. See `?run_differential_pulldown` for the full argument reference.

### Check coverage

Review mapping and filtering before interpreting a result:

```r
get_metadata(obj, "analysis_network_coverage")$comparison
get_metadata(obj, "analysis_duplicate_resolution")$comparison
get_metadata(obj, "analysis_biotype_exclusions")$comparison
```

An unmapped seed is absent from the active reference network. A mapped isolated seed is present but has no retained predicted interaction; it remains a valid result. The generic interface accepts one-direction genes. Neither outcome establishes absence of biological interactions.

### Optional validation of retained seed support

If you have an adjusted p-value column, supply its name explicitly. For example, if your input includes `adjusted_pvalue`:

```r
obj <- run_differential_pulldown(
  obj,
  label = "supported_comparison",
  input = seeds,
  order = 1,
  seed_padj_col = "adjusted_pvalue",
  seed_fdr = 0.05,
  confirm = FALSE
)
get_metadata(obj, "analysis_seed_support")$supported_comparison$gene_audit
```

This optional mode requires both mapped seed directions and at least one significant retained seed per gene, after eligibility filtering and identifier resolution. Nonsignificant complementary seeds are allowed. The supplied p-values are used as given. With `seed_padj_col = NULL` (the default), this statistical-support check is disabled.

## Plot and save results

Choose a gene identifier from your results. This example uses the first gene with mapped up and down seeds:

```r
mapped <- results$gene_id[which(results$pct_up_found > 0 & results$pct_down_found > 0)]
if (length(mapped) > 0) {
  gene_result <- get_gene_result(obj, "comparison", mapped[[1]])
  p <- plot_network(gene_result, direction = "overlay")
  print(p)
  ggplot2::ggsave("gene_network.pdf", plot = p, width = 7, height = 5)
}

saveRDS(obj, "isonet_session.rds")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

Focal seeds use purple fill and directional shapes. Neighbours and edges use green for up-only, orange for down-only (with dashed edges), and blue for shared membership. Add `level = "gene"` to merge each gene's transcripts into one node. These categories compare predicted neighbourhoods in the same reference network.

Nodes are labelled with gene symbols from `isonet_gene_names()`; each node is a transcript, so a gene can appear more than once. Use `node_labels = "id"` for transcript identifiers or `node_labels = "none"` to hide labels. `plot_component()` takes the same argument.

Reload a saved session with `obj <- readRDS("isonet_session.rds")`.

## Network eligibility and custom networks

The [network eligibility explanation](NETWORK_ELIGIBILITY.md) gives the biological rationale, annotation details and references for these defaults.

New sessions exclude transcripts annotated as nonsense-mediated decay (NMD) or non-stop decay (NSD) from both seeds and neighbours by default. To include either class, set its constructor flag explicitly:

```r
obj <- IsoNetAnalyzeR("include_nmd", include_nmd = TRUE)
obj <- IsoNetAnalyzeR("include_both", include_nmd = TRUE, include_nsd = TRUE)
```

The flags are independent. The exclusion audit is returned by `get_metadata(obj, "network_biotype_filter")`, and `isonet_transcript_biotypes()` returns the packaged transcript annotation. These options affect the network, not upstream RNA count matrices or statistical tests. Construct a new session and rerun analyses when changing this policy; saved graphs are not retroactively rebuilt.

A custom network must be an undirected `igraph` with unique vertex `name` values, non-missing vertex `gene` values, and a binary edge `interact` attribute. When either exclusion is enabled, provide known transcript biotypes through the vertex `transcript_biotype` attribute or the constructor's annotation table for IDs not covered by the packaged annotation. See `?IsoNetAnalyzeR` and the synthetic example for a complete custom network.

## Further analyses

- Use `run_component_analysis()` and `get_component_overview()` to summarize connected components across genes.
- Use `lookup_virtual_pulldown()` for a selected transcript neighbourhood.
- Use `run_network_diagnostics()` for coverage, isolation, and degree summaries.
- Use `run_fgsea()` or `run_component_ora()` for enrichment; set a random seed first, because `run_fgsea()` uses random permutations.
- Use `aggregate_duplicate_counts()` before statistical fitting if your study aggregates protein-identical transcripts. Refit the statistical models after aggregation; do not reuse adjusted p-values from the original fit.

Function details and additional examples are available through `help(package = "IsoNetAnalyzeR")` and the installed vignette. In the cross-gene summary, `n_edges_gained` and `n_edges_lost` count the edges unique to the up-seed and down-seed neighbourhoods.
