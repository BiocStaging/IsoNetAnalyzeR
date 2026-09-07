make_isolated_seed_fixture <- function() {
  vertices <- data.frame(
    name = c("tx_up", "tx_down", "tx_neighbour", "tx_isolated_up",
             "tx_isolated_down", "tx_other"),
    gene = c("G1", "G1", "G2", "G3", "G3", "G4"),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    from = c("tx_up", "tx_down"),
    to = c("tx_neighbour", "tx_neighbour"),
    interact = c(1L, 1L),
    pi = c(0.9, 0.8),
    pi_score = c(0.9, 0.8),
    stringsAsFactors = FALSE
  )
  vertices$transcript_biotype <- "protein_coding"
  network <- igraph::graph_from_data_frame(
    edges, directed = FALSE, vertices = vertices
  )
  obj <- IsoNetAnalyzeR(
    "isolated-seed-test",
    network = network,
    metadata = list(network_cutoff = 0, pi_threshold = 0.5)
  )
  input <- data.frame(
    gene_id = c("G1", "G1", "G3", "G3", "G_missing", "G_missing"),
    isoform_id = c(
      "tx_up", "tx_down", "tx_isolated_up", "tx_isolated_down",
      "tx_absent_up", "tx_absent_down"
    ),
    direction = c("up", "down", "up", "down", "up", "down"),
    stringsAsFactors = FALSE
  )
  list(obj = obj, input = input)
}
