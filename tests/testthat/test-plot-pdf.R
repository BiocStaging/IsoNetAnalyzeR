test_that("every plot type draws to a PDF device without warnings", {
  fx <- make_multi_gene_fixture()
  obj <- suppressMessages(run_component_analysis(fx$obj, "a", label = "comp"))
  obj <- suppressMessages(lookup_virtual_pulldown(obj, "look", "GA_up",
                                                  order = 1))
  analysis <- get_analysis(obj, "a")
  gene <- get_gene_result(obj, "a", "GA")

  plots <- list(
    plot_network(gene, direction = "side_by_side"),
    plot_network(gene, direction = "overlay"),
    plot_network(gene, direction = "overlay", level = "gene"),
    plot_network(gene, direction = "up", node_labels = "id"),
    plot_network(analysis, direction = "both"),
    plot_network(analysis, direction = "down", level = "gene"),
    plot_network(obj),
    plot_network(get_lookup(obj, "look"), node_labels = "none"),
    plot_component(obj, "a", "comp", 1),
    plot_component(obj, "a", "comp", 1, colour_by = "hub_degree"),
    plot_component(obj, "a", "comp", 1, colour_by = "direction"),
    plot_component(obj, "a", "comp", 1, colour_by = "direction",
                   colour_neighbours = TRUE),
    plot_cross_gene_scatter(obj, "a", "up_nodes", "down_nodes"),
    plot_cross_gene_boxplot(obj, "a", "up_nodes", "down_nodes"),
    plot_cross_gene_overview(obj, "a"),
    plot_cross_gene_boxplot_overview(obj, "a")
  )
  # The standard PDF device cannot draw characters outside its encoding,
  # such as Greek letters, and warns when a plot contains them.
  with_temp_plot_device(for (p in plots) expect_no_warning(print(p)))
})
