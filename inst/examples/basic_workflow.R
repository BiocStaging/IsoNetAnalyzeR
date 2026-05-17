# =============================================================================
# IsoNetAnalyzeR — Basic Workflow Example
# =============================================================================
#
# This script demonstrates the core functionality of the package.
# Run interactively line-by-line, not as a batch script.

library(IsoNetAnalyzeR)
library(future)

# ── 1. Create a session object ───────────────────────────────────────────────
main_obj <- IsoNetAnalyzeR("example")
print(main_obj)
print_metadata(main_obj)

# ── 2. Network summary statistics ────────────────────────────────────────────
print_network_summary(main_obj)

# Sticky genes / isoforms
sticky <- find_sticky_genes(obj, method = "sd", cutoff = 2)
sticky$threshold
head(sort(sticky$sticky_genes, decreasing = TRUE))

sticky <- find_sticky_isoforms(obj, method = "sd", cutoff = 2)
sticky$threshold
head(sort(sticky$sticky_isoforms, decreasing = TRUE))

# ── 3. Quick lookup — single seed ────────────────────────────────────────────

main_obj <- lookup_virtual_pulldown(
  main_obj,
  label = "example_virtual_pulldown",
  seed  = "ENST00000647435",
  level = "isoform",
  order = 1,
  cutoff = NULL
)
print(main_obj)

example_vp <- get_lookup(main_obj, "example_virtual_pulldown")
print(example_vp)

plot_network(example_vp,
             group_by_gene = TRUE)

# ── 4. Quick lookup — multiple seeds, gene level ─────────────────────────────
obj <- lookup_virtual_pulldown(
  obj,
  seed  = c("ENSG00000143549",
            "ENSG00000118194"),
  order = 1,
  level = "gene",
  label = "TPM3_TNNT2"
)

# ── 5. Differential pulldown analysis ────────────────────────────────────────
de_table <- data.frame(
  gene_id    = c(
    # ENSG00000004478: the hub gene — one isoform up, others down
    "ENSG00000004478", "ENSG00000004478", "ENSG00000004478",
    # ENSG00000003056: both isoforms present, split up/down
    "ENSG00000003056", "ENSG00000003056",
    # ENSG00000100083: only down-regulated hits
    "ENSG00000100083", "ENSG00000100083",
    # ENSG00000105355: no interact == 1 edges — both directions will produce NULL graphs
    "ENSG00000105355", "ENSG00000105355"
  ),
  isoform_id = c(
    "ENST00000001008", "ENST00000538622", "ENST00000539181",
    "ENST00000000412", "ENST00000536844",
    "ENST00000343632", "ENST00000325180",
    "ENST00000221957", "ENST00000585479"
  ),
  direction  = c(
    "up",   "down", "down",
    "up",   "down",
    "down", "down",
    "up",   "down"
  ),
  stringsAsFactors = FALSE
)

plan(multisession, workers = 4)

obj <- run_differential_pulldown(
  obj,
  input = de_table,
  order = 1,
  label = "condA_vs_B",
  parallel = TRUE
)

plan(sequential)

list_analyses(obj)
get_analysis(obj, "condA_vs_B")

# ── 6. Inspect per-gene results ──────────────────────────────────────────────
result <- get_analysis(obj, "condA_vs_B")
cross_gene_summary(result)

plot_network(result, direction = "both")

gene_results(result, "ENSG00000004478")

# ── 7. Plotting ──────────────────────────────────────────────────────────────
# Full network (only run when the network is small enough)
# plot_network(obj, level = "gene")

# Lookup subgraph
plot_network(get_lookup(obj, "TPM2_lookup"), group_by_gene = TRUE)

# Differential pulldown — both directions side by side
plot_network(gene_results(result, "ENSG00000100083"),
             direction = "overlay",
             group_by_gene = TRUE,
             level = "gene")


# ── 8. Session overview ──────────────────────────────────────────────────────
print(obj)
