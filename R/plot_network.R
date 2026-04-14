## plot_network.R
## Generic and methods for plotting isoform interaction networks.
## Dispatches on IsoNetAnalyzeR and VirtualPulldownLookup; shared rendering
## logic lives in internal helpers that accept plain igraph objects.
##
## Exported: plot_network

#' @include isonetanalyzer_class.R
NULL
#' @include lookup_class.R
NULL
#' @include diff_pulldown_gene_class.R
NULL

# ==============================================================================
#  Generic
# ==============================================================================

#' Plot an isoform interaction network
#'
#' Renders a network using \pkg{ggraph} and \pkg{ggplot2}.  Dispatches on
#' \linkS4class{IsoNetAnalyzeR} (full session network) and
#' \linkS4class{VirtualPulldownLookup} (single pulldown subgraph).
#'
#' @param object An \linkS4class{IsoNetAnalyzeR} or
#'   \linkS4class{VirtualPulldownLookup} object.
#' @param level Character scalar.  \code{"isoform"} (default) plots every
#'   isoform vertex.  \code{"gene"} collapses isoforms to a single node per
#'   gene, retaining only unique gene-gene edges (unweighted).
#' @param all_edges Logical.  \code{FALSE} (default) respects the
#'   \code{pi_threshold} stored in the object.  \code{TRUE} shows all edges
#'   present in the graph regardless of their \code{pi_score} attribute.
#'   See Details.
#' @param group_by_gene Logical.  Only used when \code{level = "isoform"}.
#'   Colours nodes by their parent gene when \code{TRUE}.  Default
#'   \code{FALSE}.
#' @param layout Character scalar.  Any layout supported by
#'   \code{\link[ggraph]{create_layout}}.  Defaults to \code{"stress"}, which
#'   is faster and generally better than \code{"fr"} for biological networks.
#' @param node_size Numeric scalar.  Base node size.  Default \code{3}.
#' @param label_size Numeric scalar.  Node label text size.  Default \code{2.5}.
#' @param edge_alpha Numeric in \[0, 1\].  Edge opacity.  Default \code{0.4}.
#' @param title Character scalar.  Plot title.  Auto-generated when \code{NULL}
#'   (default).
#' @param ... Additional arguments forwarded to \code{\link[ggraph]{ggraph}}.
#'
#' @details
#' **\code{all_edges}**: For \linkS4class{IsoNetAnalyzeR} the stored network
#' is already filtered at construction time, so \code{all_edges = TRUE} only
#' has an effect when a \code{pi_score} edge attribute is still present on the
#' graph.  For \linkS4class{VirtualPulldownLookup} \code{all_edges = TRUE}
#' re-displays edges whose \code{pi_score} falls below \code{@@params$cutoff},
#' which is meaningful for inspecting nodes that were filtered out.
#'
#' **Gene-level plot**: Isoforms are contracted by the \code{gene} vertex
#' attribute.  Isoforms with a missing \code{gene} attribute are grouped as
#' \code{"NA_gene"}.  Isolated gene nodes (no edges after contraction) are
#' removed.
#'
#' @return A \code{ggplot} object, invisibly.  Can be further modified with
#'   \pkg{ggplot2} layers or saved with \code{\link[ggplot2]{ggsave}}.
#'
#' @examples
#' \dontrun{
#' obj <- IsoNetAnalyzeR("demo")
#'
#' # Full network - isoform level, coloured by gene
#' plot_network(obj, group_by_gene = TRUE)
#'
#' # Full network - gene level
#' plot_network(obj, level = "gene")
#'
#' # Lookup subgraph - default (isoform, interaction edges only)
#' lkp <- lookup_virtual_pulldown(obj, seed = "ENST00000379236", order = 1)
#' plot_network(lkp)
#'
#' # Lookup subgraph - show all edges including below-cutoff ones
#' plot_network(lkp, all_edges = TRUE)
#'
#' # Lookup subgraph - gene level, save to file
#' p <- plot_network(lkp, level = "gene")
#' ggplot2::ggsave("lookup_gene.pdf", p, width = 8, height = 6)
#' }
#'
#' @export
setGeneric(
  "plot_network",
  function(object,
           level         = "isoform",
           all_edges     = FALSE,
           group_by_gene = FALSE,
           layout        = "stress",
           node_size     = 3,
           label_size    = 2.5,
           edge_alpha    = 0.4,
           title         = NULL,
           ...)
    standardGeneric("plot_network")
)


# ==============================================================================
#  Method - IsoNetAnalyzeR
# ==============================================================================

#' @describeIn plot_network Plot the full network of an
#'   \linkS4class{IsoNetAnalyzeR} session object.
setMethod(
  "plot_network",
  "IsoNetAnalyzeR",
  function(object,
           level         = "isoform",
           all_edges     = FALSE,
           group_by_gene = FALSE,
           layout        = "stress",
           node_size     = 3,
           label_size    = 2.5,
           edge_alpha    = 0.4,
           title         = NULL,
           ...) {

    level <- match.arg(level, c("isoform", "gene"))
    .check_plot_args(all_edges, group_by_gene, level)
    .check_plot_pkgs()

    # Pull the graph and threshold out of the session object, then delegate
    g <- .build_plot_graph(
      graph     = object@network,
      level     = level,
      all_edges = all_edges,
      threshold = object@metadata[["pi_threshold"]]
    )

    if (is.null(title)) {
      edge_label <- if (all_edges) "all edges" else "interaction edges"
      title <- sprintf("%s  -  %s level  (%s)", object@name, level, edge_label)
    }

    .render_plot(
      g             = g,
      level         = level,
      group_by_gene = group_by_gene,
      layout        = layout,
      node_size     = node_size,
      label_size    = label_size,
      edge_alpha    = edge_alpha,
      title         = title,
      seed_vertices = character(0),   # no seeds to highlight on the full network
      ...
    )
  }
)


# ==============================================================================
#  Method - VirtualPulldownLookup
# ==============================================================================

#' @describeIn plot_network Plot the subgraph of a
#'   \linkS4class{VirtualPulldownLookup} result.  The seed node is highlighted
#'   in a contrasting colour.
setMethod(
  "plot_network",
  "VirtualPulldownLookup",
  function(object,
           level         = "isoform",
           all_edges     = FALSE,
           group_by_gene = FALSE,
           layout        = "stress",
           node_size     = 3,
           label_size    = 2.5,
           edge_alpha    = 0.4,
           title         = NULL,
           ...) {

    level <- match.arg(level, c("isoform", "gene"))
    .check_plot_args(all_edges, group_by_gene, level)
    .check_plot_pkgs()

    # For a lookup the threshold comes from params, not session metadata
    g <- .build_plot_graph(
      graph     = object@subgraph,
      level     = level,
      all_edges = all_edges,
      threshold = object@params[["cutoff"]]
    )

    if (is.null(title)) {
      edge_label <- if (all_edges) "all edges" else "interaction edges"
      title <- sprintf(
        "Lookup: %s  -  %s level  (%s)",
        object@seed, level, edge_label
      )
    }

    # Resolve seed vertex name(s) after possible gene contraction
    seed_vertices <- if (level == "gene") {
      # After contraction the node name is the gene; look up the seed's gene
      seed_gene <- igraph::V(object@subgraph)$gene[
        igraph::V(object@subgraph)$name == object@seed
      ]
      if (length(seed_gene) && !is.na(seed_gene)) seed_gene else object@seed
    } else {
      object@seed
    }

    .render_plot(
      g             = g,
      level         = level,
      group_by_gene = group_by_gene,
      layout        = layout,
      node_size     = node_size,
      label_size    = label_size,
      edge_alpha    = edge_alpha,
      title         = title,
      seed_vertices = seed_vertices,
      ...
    )
  }
)


# ==============================================================================
#  Shared rendering - internal
# ==============================================================================

#' Core ggraph rendering shared by all plot_network methods
#'
#' @param g             igraph object already prepared by \code{.build_plot_graph}.
#' @param level         "isoform" or "gene".
#' @param group_by_gene Logical.
#' @param layout        ggraph layout string.
#' @param node_size     Numeric.
#' @param label_size    Numeric.
#' @param edge_alpha    Numeric.
#' @param title         Character scalar.
#' @param seed_vertices Character vector of seed node names to highlight.
#'   Empty vector = no highlighting.
#' @param ...           Forwarded to \code{ggraph::ggraph()}.
#'
#' @return A ggplot object, invisibly.
#' @noRd
.render_plot <- function(g, level, group_by_gene, layout,
                         node_size, label_size, edge_alpha,
                         title, seed_vertices, ...) {

  subtitle <- sprintf(
    "%d nodes  \u00b7  %d edges",
    igraph::vcount(g), igraph::ecount(g)
  )

  # Mark seed nodes as a vertex attribute so ggraph can access it via aes()
  igraph::V(g)$is_seed <- igraph::V(g)$name %in% seed_vertices

  p <- ggraph::ggraph(g, layout = layout, ...) +
    ggraph::geom_edge_link(alpha = edge_alpha, colour = "grey60") +
    .node_geom(level         = level,
               group_by_gene = group_by_gene,
               node_size     = node_size,
               has_seeds     = any(igraph::V(g)$is_seed)) +
    .label_geom(label_size = label_size) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle   = ggplot2::element_text(size = 9,  colour = "grey40"),
      legend.position = if (group_by_gene && level == "isoform") "right" else "none"
    )

  if (group_by_gene && level == "isoform") {
    p <- p + ggplot2::guides(
      colour = ggplot2::guide_legend(
        title        = "Gene",
        override.aes = list(size = 3),
        ncol         = 1
      )
    )
  }

  print(p)
  invisible(p)
}


# ==============================================================================
#  Graph preparation - internal
# ==============================================================================

#' Prepare an igraph for plotting
#'
#' Applies optional pi_score edge filtering and gene-level contraction.
#' Accepts a plain igraph and a threshold scalar, no dependency on any S4
#' class, so it can be called from any method.
#'
#' @param graph     An igraph object.
#' @param level     "isoform" or "gene".
#' @param all_edges Logical; if TRUE skip pi_score filtering.
#' @param threshold Numeric scalar or NULL; the pi_score cutoff to apply when
#'   \code{all_edges = FALSE}.
#'
#' @return A prepared igraph object with isolated vertices removed.
#' @noRd
.build_plot_graph <- function(graph, level, all_edges, threshold = NULL) {

  g <- graph

  # -- Edge filtering ----------------------------------------------------------
  # Only acts when the pi_score attribute exists AND all_edges is FALSE AND a
  # threshold is available.  If the graph is already fully filtered (no
  # pi_score attribute) this block is a no-op.
  if (!all_edges &&
      !is.null(igraph::edge_attr(g, "pi_score")) &&
      !is.null(threshold) && !is.na(threshold)) {
    keep <- igraph::E(g)$pi_score >= threshold
    g    <- igraph::subgraph.edges(g, igraph::E(g)[keep],
                                   delete.vertices = FALSE)
  }

  # -- Gene-level contraction --------------------------------------------------
  if (level == "gene") {
    g <- .contract_to_gene(g)
  }

  # -- Sanitise vertex attributes ----------------------------------------------
  # ggraph coerces vertex attributes to a tibble; NULL columns (e.g. `order`
  # and `frac_internal` absent after gene contraction) cause an error.
  # Delete any attribute whose value is NULL so ggraph only sees clean vectors.
  for (attr in igraph::vertex_attr_names(g)) {
    if (is.null(igraph::vertex_attr(g, attr)))
      g <- igraph::delete_vertex_attr(g, attr)
  }

  # Remove isolated vertices for a cleaner plot
  igraph::delete.vertices(g, igraph::degree(g) == 0L)
}


#' Contract an isoform-level igraph to gene level
#'
#' Groups vertices by the \code{gene} vertex attribute, contracts each group
#' to a single node named after the gene, and simplifies to unique unweighted
#' gene-gene edges.
#'
#' @param g igraph with a \code{gene} vertex attribute.
#' @return Simplified gene-level igraph.
#' @noRd
.contract_to_gene <- function(g) {

  gene_attr <- igraph::V(g)$gene

  if (is.null(gene_attr)) {
    warning(
      "No 'gene' vertex attribute found; plotting at isoform level instead.",
      call. = FALSE
    )
    return(g)
  }

  gene_attr[is.na(gene_attr)] <- "NA_gene"

  gene_levels <- unique(gene_attr)
  mapping     <- match(gene_attr, gene_levels)

  g_contracted <- igraph::contract(
    g,
    mapping          = mapping,
    vertex.attr.comb = list(
      name = function(x) x[1L],
      gene = function(x) x[1L],
      .    = "ignore"
    )
  )

  igraph::V(g_contracted)$name <- gene_levels

  # Drop every vertex attribute except name and gene - contraction can produce
  # NULL or list-valued attributes from the combination functions that break
  # ggraph's tibble conversion downstream.
  keep_attrs <- c("name", "gene")
  for (attr in igraph::vertex_attr_names(g_contracted)) {
    if (!attr %in% keep_attrs)
      g_contracted <- igraph::delete_vertex_attr(g_contracted, attr)
  }

  igraph::simplify(
    g_contracted,
    remove.multiple = TRUE,
    remove.loops    = TRUE,
    edge.attr.comb  = "ignore"
  )
}


# ==============================================================================
#  ggraph layer helpers - internal
# ==============================================================================

#' Node geom: uniform colour, gene-coloured, or seed-highlighted
#' @noRd
.node_geom <- function(level, group_by_gene, node_size, has_seeds) {

  # Seed highlighting takes precedence over gene colouring - if there are seed
  # nodes we use shape + fill to distinguish them; gene colouring is applied on
  # top via colour when group_by_gene is also TRUE.
  if (has_seeds) {
    if (group_by_gene && level == "isoform") {
      return(list(
        ggraph::geom_node_point(
          ggplot2::aes(
            colour = .data[["gene"]],
            size   = .data[["is_seed"]],
            shape  = .data[["is_seed"]]
          )
        ),
        ggplot2::scale_size_manual( values = c("FALSE" = node_size,
                                               "TRUE"  = node_size * 1.8)),
        ggplot2::scale_shape_manual(values = c("FALSE" = 16L, "TRUE" = 18L))
      ))
    }

    # No gene colouring - use a two-colour scheme: seed vs neighbour
    return(list(
      ggraph::geom_node_point(
        ggplot2::aes(
          fill  = .data[["is_seed"]],
          size  = .data[["is_seed"]],
          shape = .data[["is_seed"]]
        )
      ),
      ggplot2::scale_fill_manual(
        values = c("FALSE" = "#4B6FA8", "TRUE" = "#E05C3A"),
        guide  = "none"
      ),
      ggplot2::scale_size_manual(
        values = c("FALSE" = node_size, "TRUE" = node_size * 1.8),
        guide  = "none"
      ),
      ggplot2::scale_shape_manual(
        values = c("FALSE" = 21L, "TRUE" = 23L),
        guide  = "none"
      )
    ))
  }

  # No seeds to highlight
  if (group_by_gene && level == "isoform") {
    return(
      ggraph::geom_node_point(
        ggplot2::aes(colour = .data[["gene"]]),
        size = node_size
      )
    )
  }

  ggraph::geom_node_point(size = node_size, colour = "#4B6FA8")
}


#' Label geom: repelled text for all nodes
#' @noRd
.label_geom <- function(label_size) {
  ggraph::geom_node_text(
    ggplot2::aes(label = .data[["name"]]),
    size  = label_size,
    repel = TRUE
  )
}


# ==============================================================================
#  Argument / package guards - internal
# ==============================================================================

#' @noRd
.check_plot_args <- function(all_edges, group_by_gene, level) {
  if (!is.logical(all_edges)     || length(all_edges)     != 1L)
    stop("`all_edges` must be a single logical value.",     call. = FALSE)
  if (!is.logical(group_by_gene) || length(group_by_gene) != 1L)
    stop("`group_by_gene` must be a single logical value.", call. = FALSE)
  if (level == "gene" && group_by_gene)
    warning("`group_by_gene` is ignored when `level = 'gene'`.",
            call. = FALSE)
}

#' @noRd
.check_plot_pkgs <- function() {
  if (!requireNamespace("ggraph",  quietly = TRUE))
    stop("Package 'ggraph' is required: install.packages('ggraph').",
         call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required: install.packages('ggplot2').",
         call. = FALSE)
}

# ==============================================================================
#  plot_network methods
# ==============================================================================

#' Plot the subgraph(s) of a GeneDiffPulldown result
#'
#' Dispatches on \linkS4class{GeneDiffPulldown}.  Two display modes are
#' available via the \code{direction} argument:
#' \itemize{
#'   \item \code{"up"} or \code{"down"} - renders a single subgraph,
#'     highlighting the seed isoforms.
#'   \item \code{"both"} (default) - renders the two subgraphs side-by-side
#'     as a two-panel \pkg{ggplot2} / \pkg{patchwork} figure.
#' }
#'
#' All other arguments are forwarded to the underlying \code{.render_plot()}
#' helper and have the same meaning as in the
#' \linkS4class{IsoNetAnalyzeR} and \linkS4class{VirtualPulldownLookup}
#' methods.
#'
#' @param object    A \linkS4class{GeneDiffPulldown} object.
#' @param direction Character scalar: \code{"up"}, \code{"down"}, or
#'   \code{"both"} (default).
#' @param level     \code{"isoform"} (default) or \code{"gene"}.
#' @param layout    ggraph layout string.  Default \code{"stress"}.
#' @param node_size Numeric. Base node size.  Default \code{3}.
#' @param label_size Numeric. Label text size.  Default \code{2.5}.
#' @param edge_alpha Numeric in \[0, 1\].  Default \code{0.4}.
#' @param title     Character scalar or \code{NULL} (auto-generated).
#' @param ...       Forwarded to \code{ggraph::ggraph()}.
#'
#' @return A \code{ggplot} / \code{patchwork} object, invisibly.
#'
#' @examples
#' \dontrun{
#' # Plot both directions side-by-side
#' plot_network(gene_result)
#'
#' # Plot only the up-regulated subgraph
#' plot_network(gene_result, direction = "up")
#' }
#'
#' @describeIn plot_network Plot a \linkS4class{GeneDiffPulldown} result.
setMethod(
  "plot_network",
  "GeneDiffPulldown",
  function(object,
           direction  = c("both", "up", "down"),
           level      = "isoform",
           layout     = "stress",
           node_size  = 3,
           label_size = 2.5,
           edge_alpha = 0.4,
           title      = NULL,
           ...) {

    direction <- match.arg(direction)
    level     <- match.arg(level, c("isoform", "gene"))
    .check_plot_pkgs()

    if (direction == "both") {
      return(.plot_gene_both(object, level, layout,
                             node_size, label_size, edge_alpha, title, ...))
    }

    # Single-direction plot
    g_slot   <- if (direction == "up") object@up_graph else object@down_graph
    iso_seed <- if (direction == "up") object@up_isoforms else object@down_isoforms

    if (is.null(g_slot))
      stop(sprintf(
        "No %s-regulated subgraph available for gene '%s'.",
        direction, object@gene_id
      ), call. = FALSE)

    g <- .build_plot_graph(
      graph     = g_slot,
      level     = level,
      all_edges = FALSE,
      threshold = NULL        # graph is already interaction-filtered
    )

    if (is.null(title))
      title <- sprintf("%s  -  %s  (%s level)",
                       object@gene_id, direction, level)

    seed_vertices <- if (level == "gene") {
      unique(stats::na.omit(
        igraph::V(g_slot)$gene[igraph::V(g_slot)$name %in% iso_seed]
      ))
    } else {
      iso_seed
    }

    p <- .render_plot(
      g             = g,
      level         = level,
      group_by_gene = FALSE,
      layout        = layout,
      node_size     = node_size,
      label_size    = label_size,
      edge_alpha    = edge_alpha,
      title         = title,
      seed_vertices = seed_vertices,
      ...
    )

    invisible(p)
  }
)


#' Render the two-panel combined plot for a GeneDiffPulldown
#'
#' @param object     GeneDiffPulldown.
#' @param level      "isoform" or "gene".
#' @param layout     ggraph layout.
#' @param node_size  Numeric.
#' @param label_size Numeric.
#' @param edge_alpha Numeric.
#' @param title      Character or NULL.
#' @param ...        Forwarded to ggraph.
#' @return patchwork object, invisibly.
#' @noRd
.plot_gene_both <- function(object, level, layout,
                            node_size, label_size, edge_alpha, title, ...) {

  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required for direction = 'both': ",
         "install.packages('patchwork').", call. = FALSE)

  make_panel <- function(g_slot, iso_seed, dir_label) {
    if (is.null(g_slot)) {
      # Return an empty annotated ggplot when the direction has no subgraph
      return(
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = paste("No", dir_label, "subgraph"),
                            colour = "grey60", size = 4) +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste(object@gene_id, "-", dir_label))
      )
    }

    g <- .build_plot_graph(
      graph     = g_slot,
      level     = level,
      all_edges = FALSE,
      threshold = NULL
    )

    seed_vertices <- if (level == "gene") {
      unique(stats::na.omit(
        igraph::V(g_slot)$gene[igraph::V(g_slot)$name %in% iso_seed]
      ))
    } else {
      iso_seed
    }

    .render_plot(
      g             = g,
      level         = level,
      group_by_gene = FALSE,
      layout        = layout,
      node_size     = node_size,
      label_size    = label_size,
      edge_alpha    = edge_alpha,
      title         = paste(object@gene_id, "-", dir_label),
      seed_vertices = seed_vertices,
      ...
    )
  }

  p_up   <- make_panel(object@up_graph,   object@up_isoforms,   "up-regulated")
  p_down <- make_panel(object@down_graph, object@down_isoforms, "down-regulated")

  combined <- p_up + p_down +
    patchwork::plot_annotation(
      title   = title %||% sprintf("Differential pulldown  -  %s", object@gene_id),
      caption = sprintf("level: %s", level)
    )

  print(combined)
  invisible(combined)
}
