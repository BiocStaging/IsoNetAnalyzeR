# IsoNetAnalyzeR

<!-- badges: start -->
[![R CMD check on Linux](https://img.shields.io/github/check-runs/IsoformAnalysisGroup/IsoNetAnalyzeR/main?nameFilter=R-CMD-check%20%28ubuntu-latest%29&label=R%20CMD%20check%20Linux)](https://github.com/IsoformAnalysisGroup/IsoNetAnalyzeR/actions/workflows/R-CMD-check.yaml)
[![R CMD check on Windows](https://img.shields.io/github/check-runs/IsoformAnalysisGroup/IsoNetAnalyzeR/main?nameFilter=R-CMD-check%20%28windows-latest%29&label=R%20CMD%20check%20Windows)](https://github.com/IsoformAnalysisGroup/IsoNetAnalyzeR/actions/workflows/R-CMD-check.yaml)
[![R CMD check on macOS](https://img.shields.io/github/check-runs/IsoformAnalysisGroup/IsoNetAnalyzeR/main?nameFilter=R-CMD-check%20%28macos-latest%29&label=R%20CMD%20check%20macOS)](https://github.com/IsoformAnalysisGroup/IsoNetAnalyzeR/actions/workflows/R-CMD-check.yaml)
[![Bioconductor devel check](https://github.com/IsoformAnalysisGroup/IsoNetAnalyzeR/actions/workflows/bioc-check.yaml/badge.svg?branch=main)](https://github.com/IsoformAnalysisGroup/IsoNetAnalyzeR/actions/workflows/bioc-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
<!-- badges: end -->

IsoNetAnalyzeR compares the predicted protein-interaction neighbourhoods of transcript isoforms. Given isoforms whose usage goes up and isoforms whose usage goes down, for example from a differential transcript usage analysis, it expands each set in an isoform-resolved predicted interaction network and compares the resulting neighbourhoods per gene and across genes. It also provides connected-component analysis, network diagnostics, plots and enrichment analysis.

## Installation

IsoNetAnalyzeR requires R 4.6 or later. Until it is available from Bioconductor, install it from GitHub; `build_vignettes = TRUE` also builds the tutorial:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("IsoformAnalysisGroup/IsoNetAnalyzeR",
                     dependencies = TRUE, build_vignettes = TRUE)
```

Once it is available from Bioconductor:

```r
BiocManager::install("IsoNetAnalyzeR")
```

## Getting started

The vignette walks through a complete analysis of isoform switches between obese type 2 diabetic and lean adipose tissue (GEO GSE141432), using example data shipped with the package:

```r
library(IsoNetAnalyzeR)
vignette("IsoNetAnalyzeR", package = "IsoNetAnalyzeR")
```

For your own data, prepare a table with `gene_id`, `isoform_id` and `direction` (`"up"` or `"down"`) and follow the [user instructions](inst/INSTRUCTIONS.md). Statistical seed selection happens upstream; network differences are predictions to interpret alongside independent biological evidence.

## Documentation

- [Vignette source](vignettes/IsoNetAnalyzeR.Rmd)
- [User instructions](inst/INSTRUCTIONS.md)
- [Runnable synthetic example](inst/examples/full_workflow.R)
- [Network eligibility rationale](inst/NETWORK_ELIGIBILITY.md)
- [Data provenance](inst/PROVENANCE.md) and [data scripts](inst/scripts/README.md)
- [Version history](NEWS.md)

Function reference: `help(package = "IsoNetAnalyzeR")`.

## Licence

MIT; see [LICENSE.md](LICENSE.md). The bundled interaction network is derived from STRING physical interactions, which STRING distributes under CC BY 4.0; see [inst/PROVENANCE.md](inst/PROVENANCE.md) for the attribution.
