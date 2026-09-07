# IsoNetAnalyzeR

IsoNetAnalyzeR builds and compares predicted protein-interaction neighbourhoods around user-selected transcript isoforms. It provides directional network comparisons, gene and component summaries, plots, and enrichment analysis using a shared reference network.

## Documentation

**Start with the [user instructions](inst/INSTRUCTIONS.md)** for input preparation, running an analysis, plotting, saving results, and available options.

- [Runnable example script](inst/examples/full_workflow.R)
- [Worked tutorial (R vignette)](vignettes/vignette_introduction.Rmd)
- [Version history](NEWS.md)

## Installation

Install the tagged **0.6.5** release from the official repository. The [release page](https://github.com/IsoformAnalysisGroup/IsoNetAnalyzeR/releases/tag/v0.6.5) provides the source package and its SHA-256 checksum:

```r
install.packages(c("remotes", "BiocManager", "knitr", "rmarkdown"))
BiocManager::install("fgsea", ask = FALSE, update = FALSE)

remotes::install_github(
  "IsoformAnalysisGroup/IsoNetAnalyzeR@v0.6.5",
  build_vignettes = TRUE,
  upgrade = "never"
)

library(IsoNetAnalyzeR)
packageVersion("IsoNetAnalyzeR")
```

See the [installation instructions](inst/INSTRUCTIONS.md#installation) for source-archive installation and help with dependencies or vignette builds.

## Run the example

This small synthetic example requires no external dataset:

```r
source(system.file("examples", "full_workflow.R", package = "IsoNetAnalyzeR"))
```

For your own data, prepare a table containing `gene_id`, `isoform_id`, and `direction` (`"up"` or `"down"`), then follow [Run your own analysis](inst/INSTRUCTIONS.md#run-your-own-analysis).

The results compare seed-derived neighbourhoods in a static predicted network. Statistical seed selection is supplied by the user; network differences are predictions to interpret alongside independent biological evidence.

## Package reference and provenance

Use `help(package = "IsoNetAnalyzeR")` for function documentation and `vignette("vignette_introduction", package = "IsoNetAnalyzeR")` for the installed tutorial.

[Reference-data provenance](inst/PROVENANCE.md) and the [maintainer validation protocol](inst/validation/VALIDATION_PROTOCOL.md) document the packaged resources and release checks. The package is distributed under the [MIT licence](LICENSE).
