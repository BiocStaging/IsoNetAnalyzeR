# Network reduction equivalence report

Generated 2026-10-30 22:57. Full: 0.6.5 network (SHA-256 3130209f...3ebb). 
Reduced: interact == 1 edges + all vertices (SHA-256 41df4ce0...bf9b).
Both runs use the identical 0.6.5 code; only R/sysdata.rda differs.

**43 checks, 43 PASS, 0 FAIL.**

| Area | Check | Result | Detail |
| --- | --- | --- | --- |
| session | vertex count | PASS | 46433 vs 46433 |
| session | vertex table, in order | PASS | identical |
| session | eligibility audit (vertex fields) | PASS | identical |
| session | eligibility audit edge counts (expected change) | PASS | edges_before 2853266 -> 158458; edges_after 1793213 -> 146267 |
| session | session edge count (expected change) | PASS | 1793213 -> 146267 |
| session | pi_threshold | PASS | identical |
| session | network_cutoff | PASS | identical |
| session | network_source | PASS | identical |
| interact graph | vertices, in order | PASS | identical |
| interact graph | edges with pi/interact, in order | PASS | identical |
| interact graph | size | PASS | 46433 vertices, 146267 edges |
| 65-gene grid | o1_c0 (65 genes: graphs, stats, summaries, audits) | PASS | identical |
| 65-gene grid | o1_c0.2609 (65 genes: graphs, stats, summaries, audits) | PASS | identical |
| 65-gene grid | o1_c0.5 (65 genes: graphs, stats, summaries, audits) | PASS | identical |
| 65-gene grid | o2_c0 (65 genes: graphs, stats, summaries, audits) | PASS | identical |
| 65-gene grid | o2_c0.2609 (65 genes: graphs, stats, summaries, audits) | PASS | identical |
| 65-gene grid | o2_c0.5 (65 genes: graphs, stats, summaries, audits) | PASS | identical |
| components | all_gene_components | PASS | identical |
| components | iso_components | PASS | identical |
| enrichment | component ORA (GO:BP 2026.1) | PASS | identical |
| enrichment | fgsea (GO:BP 2026.1, seed 141432) | PASS | identical |
| golden (frozen study) | full network: all 65 gene graphs | PASS | 0 mismatching graphs |
| golden (frozen study) | full network: CLTCL1 8 up / 33 down | PASS | 8 up / 33 down |
| golden (frozen study) | full network: gene-level components | PASS |  |
| golden (frozen study) | reduced network: all 65 gene graphs | PASS | 0 mismatching graphs |
| golden (frozen study) | reduced network: CLTCL1 8 up / 33 down | PASS | 8 up / 33 down |
| golden (frozen study) | reduced network: gene-level components | PASS |  |
| random seed sets | same 500 specifications | PASS | identical |
| random seed sets | virtual_pulldown outputs identical | PASS | 0 mismatches; 471 graphs with edges, 29 edgeless, 0 errors (same message on both) |
| lookups | 10 isoform + 10 gene lookups | PASS | identical |
| diagnostics | network_stats except full_edges | PASS | identical |
| diagnostics | full_edges (expected change) | PASS | 1793213 -> 146267 |
| diagnostics | isoform_degree_table | PASS | identical |
| diagnostics | gene_degree_table | PASS | identical |
| diagnostics | component_table | PASS | identical |
| diagnostics | clustering_table | PASS | identical |
| diagnostics | isoform_per_gene_table | PASS | identical |
| diagnostics | sticky_isoform_table | PASS | identical |
| diagnostics | sticky_gene_table | PASS | identical |
| diagnostics | sticky_method | PASS | identical |
| diagnostics | sticky_cutoff | PASS | identical |
| diagnostics | sticky_isoform_threshold | PASS | identical |
| diagnostics | sticky_gene_threshold | PASS | identical |

## Timings (seconds)

| Step | Full | Reduced |
| --- | --- | --- |
| session | 0.7 | 0.2 |
| filter_interact | 1.2 | 0.1 |
| pulldown_o1_c0 | 158.2 | 17.2 |
| pulldown_o1_c0.2609 | 156.6 | 17 |
| pulldown_o1_c0.5 | 141.7 | 17.4 |
| pulldown_o2_c0 | 152.2 | 17.3 |
| pulldown_o2_c0.2609 | 156.1 | 17.1 |
| pulldown_o2_c0.5 | 144.5 | 16.2 |
| components_gene | 0.3 | 0.3 |
| components_isoform | 0.2 | 0.2 |
| component_ora | 7.3 | 6.9 |
| fgsea | 0.7 | 0.7 |
| random_500 | 502 | 66.5 |
| diagnostics | 38.6 | 40.5 |
