test_that("cross-gene scatter plots draw one point per gene", {
  fx <- make_multi_gene_fixture()
  p <- plot_cross_gene_scatter(fx$obj, "a", "up_nodes", "down_nodes")
  expect_s3_class(p, "ggplot")
  expect_identical(nrow(p$data), 6L)

  p_log <- plot_cross_gene_scatter(
    fx$obj, "a", "up_edges", "down_edges", x_axis_log = TRUE,
    y_axis_log = TRUE, label = c("delta", "x"), label_n = 3L, title = "Edges"
  )
  expect_s3_class(p_log, "ggplot")
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
    expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
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
