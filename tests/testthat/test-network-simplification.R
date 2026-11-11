# Seed S has one neighbour, N; N also interacts with X1 and X2. The S-N pair
# is listed twice, N-X2 twice (one copy not predicted) and N with itself.
make_repeated_edge_network <- function() {
  edges <- data.frame(
    from     = c("S", "S", "N",  "N",  "N",  "N"),
    to       = c("N", "N", "X1", "X2", "X2", "N"),
    interact = c(1L,  1L,  1L,   0L,   1L,   1L),
    pi       = c(0.2, 0.6, 0.5,  0.1,  0.3,  0.9)
  )
  vertices <- data.frame(name = c("S", "N", "X1", "X2"),
                         gene = c("G1", "G2", "G3", "G4"),
                         transcript_biotype = "protein_coding")
  igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
}

test_that("supplied networks are simplified when a session is created", {
  expect_message(
    obj <- IsoNetAnalyzeR("repeated", network = make_repeated_edge_network(),
                          metadata = list(network_cutoff = 0.4)),
    "Merged 2 repeated edge(s) and removed 1 self-loop(s)", fixed = TRUE
  )
  net <- get_network(obj)
  expect_true(igraph::is_simple(net))
  expect_equal(igraph::ecount(net), 3)
  sn <- igraph::get_edge_ids(net, c("S", "N"))
  nx2 <- igraph::get_edge_ids(net, c("N", "X2"))
  # S-N: both copies predicted, pi 0.2 and 0.6
  expect_equal(igraph::E(net)$interact[sn], 1)
  expect_equal(igraph::E(net)$pi[sn], 0.4)
  # N-X2: one copy predicted, pi 0.1 and 0.3
  expect_equal(igraph::E(net)$interact[nx2], 1)
  expect_equal(igraph::E(net)$pi[nx2], 0.2)
})

test_that("repeated edges and self-loops do not move the neighbourhood cutoff", {
  obj <- suppressMessages(IsoNetAnalyzeR(
    "repeated", network = make_repeated_edge_network(),
    metadata = list(network_cutoff = 0.4)
  ))
  # N has 3 predicted edges, 1 of them inside the order-1 neighbourhood
  # {S, N}: frac_internal = 1/3. Counting the repeated S-N edge and the N-N
  # loop (degree 2) would give (2 + 2) / (2 + 1 + 1 + 2) = 4/6 and keep N at
  # cutoff 0.4.
  kept <- suppressMessages(virtual_pulldown(obj, "S", order = 1,
                                            cutoff = 0.4))
  expect_identical(igraph::V(kept)$name, "S")
  kept <- suppressMessages(virtual_pulldown(obj, "S", order = 1,
                                            cutoff = 0.3))
  expect_setequal(igraph::V(kept)$name, c("S", "N"))
  expect_equal(igraph::V(kept)$frac_internal[igraph::V(kept)$name == "N"],
               1 / 3)
})

test_that("simple networks are used unchanged", {
  fx <- make_isolated_seed_fixture()
  expect_true(igraph::is_simple(get_network(fx$obj)))
  expect_equal(igraph::ecount(get_network(fx$obj)), 2)
})
