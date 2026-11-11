# Data shipped with IsoNetAnalyzeR

This directory holds the scripts that produced each data file in the package.
Full provenance, including source checksums, is in `inst/PROVENANCE.md`.

| File | Contents | Produced by |
| --- | --- | --- |
| `R/sysdata.rda` | Bundled isoform interaction network (58,164 transcripts, 158,447 predicted interaction edges), neighbourhood cutoff and pi threshold | `data-raw/update_network.R` (full scored network, 0.6.5), then `reduce_network.R` |
| `inst/extdata/protein_identical_isoform_map.csv` | 16,630 redundant transcripts mapped to 9,760 protein-identical representatives | `data-raw/build_duplicate_map.R` in the source repository |
| `inst/extdata/protein_identical_isoform_map_provenance.csv` | Build checksums and counts for the map | `data-raw/build_duplicate_map.R` |
| `inst/extdata/gse141432_t2d_vs_lean_seeds.csv` | Example isoform seeds: 137 transcripts of 65 genes, Obese-T2D versus Lean, GEO GSE141432 | `make_example_data.R` |

`data-raw/` is excluded from the built package but is part of the source
repository.

## Scripts

- `reduce_network.R` keeps the predicted interaction edges (`interact == 1`)
  of the 0.6.5 network together with every vertex, merges the six
  transcript pairs that the source table lists twice and removes five
  self-pairs. All analyses use only predicted interaction edges; the scored
  non-interacting pairs are not shipped.
- `make_example_data.R` writes the GSE141432 seed table from the frozen
  annotated seed table.

Transcript biotypes and gene symbols come from Ensembl 107 (AnnotationHub
record AH104864), and the gene sets used in the examples and the vignette
from the msigdb package; neither is shipped.
