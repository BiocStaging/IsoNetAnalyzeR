# Data shipped with IsoNetAnalyzeR

This directory holds the scripts that produced each data file in the package.
Full provenance, including source checksums, is in `inst/PROVENANCE.md`.

| File | Contents | Produced by |
| --- | --- | --- |
| `R/sysdata.rda` | Bundled isoform interaction network (58,164 transcripts, 158,458 predicted interaction edges), neighbourhood cutoff and pi threshold | `data-raw/update_network.R` (full scored network, 0.6.5), then `reduce_network.R` |
| `inst/extdata/protein_identical_isoform_map.csv` | 16,630 redundant transcripts mapped to 9,760 protein-identical representatives | `data-raw/build_duplicate_map.R` in the source repository |
| `inst/extdata/protein_identical_isoform_map_provenance.csv` | Build checksums and counts for the map | `data-raw/build_duplicate_map.R` |
| `inst/extdata/transcript_biotypes_gencode_v41.csv.gz` | GENCODE v41 transcript biotypes for 250,161 transcripts | Exported from the project's GENCODE v41 isoform annotation (see `inst/PROVENANCE.md`) |
| `inst/extdata/gene_names_gencode_v41.csv.gz` | GENCODE v41 gene symbols for the 10,554 genes of the bundled network, used to label plots | `make_gene_names.R` |
| `inst/extdata/gse141432_t2d_vs_lean_seeds.csv` | Example isoform seeds: 137 transcripts of 65 genes, Obese-T2D versus Lean, GEO GSE141432 | `make_example_data.R` |
| `inst/extdata/msigdb_hallmark_ensembl.rds` | MSigDB Hallmark gene sets (50 sets, Ensembl gene IDs, CC BY 4.0) for offline examples | `make_example_data.R` |

`data-raw/` is excluded from the built package but is part of the source
repository.

## Scripts

- `reduce_network.R` keeps the predicted interaction edges (`interact == 1`)
  of the 0.6.5 network together with every vertex. All analyses use only these
  edges, so results are unchanged; the scored non-interacting pairs are not
  shipped.
- `make_gene_names.R` extracts gene symbols from the project's GENCODE v41
  isoform annotation, the source of the biotype table, for the genes of the
  bundled network.
- `make_example_data.R` writes the two example files from the frozen
  GSE141432 seed table and `msigdbr`.
