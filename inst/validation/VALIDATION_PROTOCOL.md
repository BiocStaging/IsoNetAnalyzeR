# IsoNetAnalyzeR 0.3.2 validation protocol

## Automated package checks

1. Run `testthat::test_local("source_fork")` and archive the complete result.
2. Run `R CMD build source_fork` with vignettes enabled.
3. Run `R CMD check --no-manual <tarball>` and require `Status: OK`.
4. With LaTeX, qpdf, current HTML Tidy, and V8 available, run
   `R CMD check --as-cran <tarball>`. Archive `00check.log`, the generated
   manual, and every remaining finding.
5. Verify the built package version is 0.3.2 and includes this validation
   directory and `NEWS.md`.

## Required regression cases

- ENSG bookkeeping at gene level and ENST bookkeeping at isoform level.
- Interacting mapped seed, mapped isolated seed, and unmapped seed.
- Same mapped seed in both directions.
- Connected and singleton components at isoform and gene levels.
- Explicit `min_edges` removal of singleton components.
- Single-direction, overlay, lookup, and component plotting for a singleton.
- Filtered and singleton plotting emits no igraph deprecation warnings while
  preserving isolated seeds and the retained graph inventory.
- Save/read RDS round trip followed by component access and plotting.
- Explicit seed flags for disconnected multi-seed graphs.
- Legacy node-stat tables without an explicit `is_seed` column.
- Invalid labels, orders, cutoffs, required columns, directions, lookup seeds,
  component references, and diagnostic thresholds.
- Fully unmapped inputs retain typed empty statistics and component tables.
- Single-direction inputs remain reportable.
- Plot tests use an explicit temporary device and leave no `Rplots.pdf`.

## Real-data repair-only validation

Use the saved GSE141432 and GSE231509 pulldown sessions as frozen graph inputs
with the unchanged pinned network. Regenerate component analyses under a new
label without overwriting the historical component labels. This isolates the
component repair from any change in differential inputs or pulldown expansion.
Compare the historical and repaired component inventories.

Acceptance criteria:

- 24 previously omitted gene-level singleton components are restored in total;
- the multiset of node and edge counts for previously valid non-singleton
  components is unchanged;
- former one-node overlay failures render successfully;
- genuinely absent directions are recorded as coverage limitations and are
  never counted as plotting failures;
- every mapped/unmapped input, per-direction loss, singleton, and mixed
  component is reported.

## Interpretation guardrail

All components are inventoried. Rewiring claims are restricted to nontrivial
components containing mapped seed support in both directions. A singleton is a
valid loss-of-interaction result, not sufficient evidence of comparative
rewiring by itself.

## Release artifacts

- source tarball with SHA-256;
- Git patch against the `source_fork` base commit with SHA-256;
- issue register and audited change log;
- detailed repair report in Markdown and standalone HTML;
- automated-test output, package build log, and `R CMD check` log;
- real-data before/after component comparison and checksums.
