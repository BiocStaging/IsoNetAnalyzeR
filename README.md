# IsoNetAnalyzeR

IsoNetAnalyzeR compares the predicted protein-interaction neighbourhoods of transcript isoforms. Given isoforms whose usage goes up and isoforms whose usage goes down, for example from a differential transcript usage analysis, it expands each set in an isoform-resolved predicted interaction network and compares the resulting neighbourhoods per gene and across genes. It also provides connected-component analysis, network diagnostics, plots and enrichment analysis.

## Installation

From Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IsoNetAnalyzeR")
```

Development version from GitHub:

```r
BiocManager::install("IsoformAnalysisGroup/IsoNetAnalyzeR")
```

Earlier releases, such as [0.6.5](https://github.com/IsoformAnalysisGroup/IsoNetAnalyzeR/releases/tag/v0.6.5), remain available as tagged GitHub releases.

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
- [Maintainer validation protocol](inst/validation/VALIDATION_PROTOCOL.md)
- [Version history](NEWS.md)

Function reference: `help(package = "IsoNetAnalyzeR")`.

## Licence

MIT; see [LICENSE.md](LICENSE.md). The bundled interaction network is derived from STRING physical interactions, which STRING distributes under CC BY 4.0; see [inst/PROVENANCE.md](inst/PROVENANCE.md) for the attribution.
