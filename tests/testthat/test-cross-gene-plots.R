# The multi-gene fixture (helper-fixtures.R) has six genes.
test_that("cross-gene scatter plots draw one point per gene", {
  fx <- make_multi_gene_fixture()
  p <- plot_cross_gene_scatter(fx$obj, "a", "up_nodes", "down_nodes")
  expect_true(ggplot2::is_ggplot(p))
  expect_identical(nrow(p$data), 6L)

  p_log <- plot_cross_gene_scatter(
    fx$obj, "a", "up_edges", "down_edges", x_axis_log = TRUE,
    y_axis_log = TRUE, label = c("delta", "x"), label_n = 3L, title = "Edges"
  )
  expect_true(ggplot2::is_ggplot(p_log))
  with_temp_plot_device({
    expect_no_warning(print(p))
    expect_no_warning(print(p_log))
  })
})

test_that("cross-gene boxplots and overviews draw", {
  fx <- make_multi_gene_fixture()
  plots <- list(
    plot_cross_gene_boxplot(fx$obj, "a", "up_nodes", "down_nodes"),
    plot_cross_gene_boxplot(fx$obj, "a", "up_edges", "down_edges",
                            log_ratio = TRUE, label = c("delta", "up"),
                            label_n = 2L),
    plot_cross_gene_overview(fx$obj, "a"),
    plot_cross_gene_boxplot_overview(fx$obj, "a", log_ratio = TRUE)
  )
  for (p in plots)
    expect_true(ggplot2::is_ggplot(p) || inherits(p, "patchwork"))
  with_temp_plot_device(for (p in plots) expect_no_warning(print(p)))
})

test_that("cross-gene plots reject unknown or non-numeric columns", {
  fx <- make_multi_gene_fixture()
  expect_error(plot_cross_gene_scatter(fx$obj, "a", "up_nodes", "nonexistent"),
               "not found")
  expect_error(plot_cross_gene_scatter(fx$obj, "a", "up_nodes", "gene_id"),
               "numeric")
  expect_error(plot_cross_gene_boxplot(fx$obj, "a", "up_nodes", "nonexistent"),
               "not found")
  expect_error(plot_cross_gene_scatter(fx$obj, "a", "up_nodes", "down_nodes",
                                       label = "bogus"))
})

test_that("each label strategy labels its own top genes", {
  df <- data.frame(gene_id = paste0("g", 1:5), x__ = c(50, 40, 3, 2, 1),
                   y__ = c(49, 39, 30, 2, 1))
  df$delta__ <- df$x__ - df$y__
  # x - y is 1, 1, -27, 0, 0: g3 is furthest from the diagonal, g1 has the
  # highest x and (first of two) the largest positive difference
  expect_identical(.select_labels(df, "x", 1L), "g1")
  expect_identical(.select_labels(df, "delta", 1L), "g3")
  expect_identical(.select_labels(df, "up", 1L), "g1")
  expect_identical(.select_labels(df, "down", 1L), "g3")
  # Strategies take turns in the order given: x's best (g1), then delta's
  # best (g3), and the other way round
  expect_identical(.select_labels(df, c("x", "delta"), 2L), c("g1", "g3"))
  expect_identical(.select_labels(df, c("delta", "x"), 2L), c("g3", "g1"))
  # "up" takes only positive differences (g1, g2), "down" only negative (g3)
  expect_identical(.select_labels(df, "up", 3L), c("g1", "g2"))
  expect_identical(.select_labels(df, "down", 3L), "g3")

  # Differences -5..-1 and 1..5 for g1..g10: "up" ranks g10, g9, ... and
  # "down" g1, g2, ...; taking turns gives g10, g1, g9, g2, g8, g3
  ten <- data.frame(gene_id = paste0("g", 1:10), x__ = 0, y__ = 0,
                    delta__ = c(-5:-1, 1:5))
  expect_identical(.select_labels(ten, c("up", "down"), 6L),
                   c("g10", "g1", "g9", "g2", "g8", "g3"))
})

test_that("up and down have the same colours in every plot", {
  # Up is green and down orange, as in the overlay plot
  expect_identical(unname(.membership_palette[c("up_only", "down_only")]),
                   c("#009E73", "#D55E00"))
  fx <- make_multi_gene_fixture()
  obj <- suppressMessages(run_component_analysis(fx$obj, "a", "comp"))
  p <- suppressMessages(plot_component(obj, "a", "comp", 1,
                                       colour_by = "direction",
                                       node_labels = "none"))
  nodes <- ggplot2::ggplot_build(p)$data[[2L]]
  role <- as.character(p$data$colour_role)
  expect_identical(unique(nodes$colour[role == "up"]), "#009E73")
  expect_identical(unique(nodes$colour[role == "down"]), "#D55E00")
  # Box plots: up_nodes - down_nodes is GA 2, GB -2, GC 1, GD 0, GE 2,
  # GF -1, so GA, GC and GE are green, GB and GF orange, and GD grey; the
  # keys name all three
  box <- plot_cross_gene_boxplot(fx$obj, "a", "up_nodes", "down_nodes")
  points <- ggplot2::ggplot_build(box)$data[[3L]]
  expect_identical(sort(table(points$colour)),
                   sort(table(c(rep("#009E73", 3), rep("#D55E00", 2),
                                "#999999"))))
  expect_match(box$labels$subtitle, paste0(
    "^Green = up_nodes > down_nodes.  Orange = up_nodes < down_nodes.  ",
    "Grey = equal."))
  overview <- plot_cross_gene_boxplot_overview(fx$obj, "a")
  expect_match(overview$patches$annotation$subtitle,
               "Green = up > down.  Orange = up < down.  Grey = equal.",
               fixed = TRUE)
})

test_that("plot_component names the available component analyses", {
  fx <- make_multi_gene_fixture()
  # The fixture's analysis has no component analysis yet
  expect_error(plot_component(fx$obj, "a", "missing", 1), "Available: none")
})

test_that("a missing component analysis label gives an argument error", {
  fx <- make_multi_gene_fixture()
  expect_error(plot_component(fx$obj, "a", NA_character_, 1),
               "`component_analysis` must be a single non-empty")
})

test_that("a missing analysis label gives an argument error", {
  fx <- make_multi_gene_fixture()
  expect_error(get_component_overview(fx$obj, NA_character_, "comp"),
               "`analysis` must be a single non-empty")
})

test_that("group_by_gene = NA gives an argument error", {
  fx <- make_multi_gene_fixture()
  gr <- get_gene_result(fx$obj, "a", "GA")
  expect_error(plot_network(gr, group_by_gene = NA),
               "`group_by_gene` must be a single logical value")
})

test_that("cross-gene plot_network() checks its arguments like the others", {
  fx <- make_multi_gene_fixture()
  analysis <- get_analysis(fx$obj, "a")
  expect_error(plot_network(analysis, group_by_gene = NA),
               "`group_by_gene` must be a single logical value")
})

test_that("an empty cross-gene summary warning names the analysis", {
  fx <- make_multi_gene_fixture()
  obj <- fx$obj
  obj@analyses$a@cross_gene_summary <- obj@analyses$a@cross_gene_summary[0, ]
  expect_warning(get_cross_gene_summary(obj, "a"),
                 "Analysis 'a' has an empty cross-gene summary", fixed = TRUE)
})
