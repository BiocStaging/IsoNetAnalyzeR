# Network eligibility: NMD and NSD

This document explains the transcript eligibility policy used by IsoNetAnalyzeR. For the analysis workflow, see the [user instructions](INSTRUCTIONS.md#network-eligibility-and-custom-networks).

The default introduced in 0.5.5 is conservative for protein-interaction interpretation:

```r
obj <- IsoNetAnalyzeR("protein_network", include_nmd = FALSE, include_nsd = FALSE)
# Explicit opt-in for exploratory or legacy comparisons:
all_transcripts <- IsoNetAnalyzeR("unfiltered", include_nmd = TRUE, include_nsd = TRUE)
# The two flags are independent.
annotation <- isonet_transcript_biotypes()
```

`nonsense_mediated_decay` (NMD) and `non_stop_decay` (NSD) are **transcript** biotypes. A gene labelled `protein_coding`, an annotated ORF, or a predicted interaction does not establish a stable functional protein product for every transcript. GENCODE associates these decay classes with premature termination or absence of a stop codon, respectively. Excluding them by default reduces the risk of interpreting decay-associated RNA changes as changes between stable protein interaction partners. This is an annotation-based eligibility policy, not proof that these RNAs are never translated. NMD is translation-dependent, and protein products can be made and subsequently destabilized. [GENCODE biotype definitions](https://www.gencodegenes.org/pages/biotypes.html), [experimental evidence](https://pmc.ncbi.nlm.nih.gov/articles/PMC8694075/).

The session constructor removes excluded vertices and their incident edges **before** neighbourhood expansion and internal-edge fraction filtering. Thus the same policy covers seeds and neighbours. Submitted decay-associated aliases are blocked before protein-identical remapping; they cannot bypass exclusion through a coding representative. The vertex audit is returned by `get_metadata(obj, "network_biotype_filter")`; excluded submitted seed rows are in `get_metadata(obj, "analysis_biotype_exclusions")[[label]]`. Unknown transcript biotypes produce an actionable error if either exclusion is enabled. An explicit table overrides vertex annotations, which override the packaged GENCODE v41 mapping; gene-level biotypes are never substituted.

These options do not filter RNA count matrices, recompute isoform fractions, or change pairedGSEA tests. Keep decay-associated transcripts in the RNA analysis when the question concerns total transcript usage, then apply protein eligibility before network seed selection. Removing them from the RNA denominator would change that question. Recheck bidirectionality after eligibility and mapping; do not silently fill a missing direction with an ineligible transcript.

All remaining vertices are only eligible predictions, not experimentally validated proteins. Other biotypes require their own interpretation. Removing vertices also changes the denominator of the network's internal-edge fraction, so a previously excluded coding neighbour may cross the same cutoff. The inherited network cutoff was calibrated on the original network; reusing it in a filtered-network sensitivity run does not establish that it remains optimal.

Existing saved sessions and their cached results are not retroactively filtered. Construct a new session from the original network and rerun analyses when changing policy. Do not try to restore excluded vertices merely by changing metadata on an already filtered graph. Historical versions through 0.5.0 retain their unfiltered defaults; releases from 0.5.5 onwards include the exclusions.
