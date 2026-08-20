# IsoNetAnalyzeR 0.5.0 release verification

- Canonical tarball: `IsoNetAnalyzeR_0.5.0.tar.gz`
- Tarball SHA-256: `52be4416788ae621730ba043484a11a26b9cf3bd8ed37ae1e3c253d016a828d8`
- Source map SHA-256: `cfe6fb57ed06c7231b9c77d7246b148be3cfdb1ac6aad4f9895cc422efaa36f5`
- Enriched packaged map SHA-256: `95ff0ac877397f8c0e99d2789f2b833dec8a01009b6d2686170950da3f56c85a`
- Embedded `R/sysdata.rda` SHA-256: `3130209f090b9b83bbf2ebf328f29a0ad9a8fcd89922a277d0966325bc163ebb`
- Pinned GSE141432 evidence SE SHA-256: `8ab52523b1564b061c9424e5189dcf6014f0f4e0161d1eb1e6f87c855dbe7b35`

Verification results:

- complete testthat suite: passed;
- clean tarball install and installed-resource gate: passed;
- `R CMD check --no-manual`: `Status: OK`;
- `R CMD check --as-cran`: code, tests, examples, HTML manual, and vignette
  rebuild passed. The run could not build the PDF manual because `pdflatex` is
  not installed. Its other notes are the expected unavailable-network and
  unverifiable-clock checks in the restricted environment.

The full check logs are archived beside this manifest.
