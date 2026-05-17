# IsoNetAnalyzeR

An R package for isoform-resolved network analysis of differential splicing data. IsoNetAnalyzeR maps differentially expressed isoforms onto a protein interaction network, runs virtual pulldowns to extract condition-specific subnetworks, and compares their topology across genes and conditions.

## Overview

Alternative splicing can rewire protein interaction networks without changing gene expression. IsoNetAnalyzeR makes this rewiring visible and quantifiable. Given a table of up- and down-regulated isoforms, it:

-   Expands each isoform seed set through a predicted interaction network via BFS (virtual pulldown)
-   Computes per-node and subgraph-level topology statistics for each condition
-   Summarises rewiring across all genes in a single 35-column table
-   Decomposes the merged cross-gene network into connected components and classifies their topology
-   Runs gene set enrichment (GSEA) on topology delta scores and over-representation analysis (ORA) per component
-   Provides diagnostic plots for the underlying network and all result layers

## Installation

``` r
# Install from GitHub
remotes::install_github("your-org/IsoNetAnalyzeR")
```

Dependencies: `igraph`, `cli`, `ggplot2`, `ggraph`, `ggrepel`, `patchwork`, `fgsea`, `org.Hs.eg.db`, `AnnotationDbi`.\
Optional (parallel execution): `future`, `future.apply`.\
Optional (Excel input): `readxl`.

## Quick start

``` r
library(IsoNetAnalyzeR)

# 1. Initialise session — loads the built-in isoform interaction network
obj <- IsoNetAnalyzeR("my_experiment")

# 2. Run differential pulldown
#    Input: CSV / Excel / data.frame with columns gene_id, isoform_id, direction
obj <- run_differential_pulldown(obj, "cond_A", "de_isoforms.csv", order = 1)

# 3. Inspect cross-gene summary
cgs <- get_cross_gene_summary(obj, "cond_A")
head(cgs[order(-abs(cgs$n_edges_gained - cgs$n_edges_lost)), ])

# 4. Visualise a single gene
gr <- get_gene_result(obj, "cond_A", "ENSG00000120333")
plot_network(gr)

# 5. Run component analysis on the merged cross-gene network
obj <- run_component_analysis(obj, "cond_A", label = "comp_both")

# 6. Explore components
overview <- get_component_overview(obj, "cond_A", "comp_both")
overview[overview$mixed, ]   # rewiring candidates

# 7. GSEA on topology deltas
obj <- run_fgsea(obj, "cond_A", label = "gsea_nodes", rank_by = "delta_nodes")
get_gsea_table(obj, "cond_A", "gsea_nodes", padj_cutoff = 0.05)
```

## Core workflow

```         
IsoNetAnalyzeR()                    # create session, load network
  └── run_differential_pulldown()   # per-gene virtual pulldown (up vs down)
        ├── get_cross_gene_summary()    # 35-col table: topology, Jaccard, edge gain/loss
        ├── get_isoform_stats()         # per-isoform: degree, clustering, pi confidence
        ├── plot_network()              # subgraph visualisation (gene, analysis, or lookup)
        │
        ├── run_component_analysis()    # connected components of merged cross-gene graph
        │     ├── get_component_overview()   # Tier-1 + Tier-2 stats per component
        │     ├── get_component()            # per-node table for one component
        │     ├── plot_component()           # direction territory / dominance plot
        │     └── run_component_ora()        # GO over-representation per component
        │           └── get_ora_table()
        │
        └── run_fgsea()              # GSEA on delta-topology gene ranking
              └── get_gsea_table()
              └── get_gsea_stats()

lookup_virtual_pulldown()           # quick single-seed pulldown (independent)
run_network_diagnostics()           # degree distribution, sticky nodes, components
plot_network_diagnostics()
```

## Input format

A long-format table with exactly three required columns:

| Column       | Description                                                  |
|----------------------------|--------------------------------------------|
| `gene_id`    | Gene identifier (e.g. Ensembl gene ID)                       |
| `isoform_id` | Isoform / transcript ID matching vertex names in the network |
| `direction`  | `"up"` or `"down"` (case-insensitive)                        |

Accepted formats: CSV, Excel (`.xlsx`/`.xls`), or an R `data.frame`.

## Key statistics produced

### Cross-gene summary (`get_cross_gene_summary`)

One row per gene. Key columns include:

| Column group | Columns |
|------------------------------------|------------------------------------|
| Input coverage | `n_up_input`, `n_down_input`, `pct_up_found`, `pct_down_found` |
| Subgraph topology | `up_nodes`, `up_edges`, `up_density`, `up_avg_degree`, `up_clustering`, `up_starness` (+ down equivalents) |
| Node-set overlap | `n_shared_nodes`, `n_up_unique`, `n_down_unique`, `jaccard` |
| Edge rewiring | `n_edges_gained`, `n_edges_lost`, `n_edges_shared` |
| Interaction confidence | `up_seed_mean_pi`, `up_seed_pi_entropy` (+ down equivalents) |

### Component overview (`get_component_overview`)

One row per connected component. Tier-1 (structural) + Tier-2 (topological) columns joined automatically:

| Tier | Key columns |
|------------------------------------|------------------------------------|
| Tier-1 | `n_nodes`, `n_edges`, `density`, `n_seeds_up`, `n_seeds_down`, `mixed`, `shared_fraction` |
| Tier-2 | `diameter`, `radius`, `degree_centralization`, `topology_class` (`"hub-spoke"` / `"chain-like"` / `"clique-like"`) |

### Network diagnostics (`get_network_diagnostics`)

19 global statistics plus degree tables, component table, clustering table, isoform-per-gene connectivity, and sticky node tables.

## Parallel execution

``` r
library(future)
plan(multisession, workers = 4)

obj <- run_differential_pulldown(obj, "cond_A", input, order = 1, parallel = TRUE)

plan(sequential)
```

## Session metadata

``` r
obj <- set_metadata(obj, "author", "Jane Doe")
obj <- set_metadata(obj, "notes",  "Pilot run — March 2025")
print_metadata(obj)
```
