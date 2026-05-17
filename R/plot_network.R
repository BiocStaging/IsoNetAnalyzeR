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
#' @include diff_pulldown_analysis_class.R
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
    # Apply CVD-safe gene colour palette
    gene_ids <- igraph::V(g)$gene
    if (!is.null(gene_ids)) {
      gene_ids[is.na(gene_ids)] <- "NA_gene"
      unique_genes <- unique(gene_ids)
      gene_cols    <- .gene_palette(length(unique_genes))
      names(gene_cols) <- unique_genes
      p <- p + ggplot2::scale_colour_manual(values = gene_cols)
    }
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
#  Colour-blind-safe palettes - internal
# ==============================================================================

#' Okabe-Ito-inspired palette for membership categories
#'
#' Three colours chosen for deuteranopia/protanopia/tritanopia safety.
#' @noRd
.membership_palette <- c(
  up_only   = "#E69F00",   # orange
  shared    = "#999999",   # grey
  down_only = "#56B4E9"    # sky blue
)

#' Edge colours matching the membership palette
#' @noRd
.edge_membership_palette <- .membership_palette

#' Seed/neighbour two-colour palette (CVD-safe)
#' @noRd
.seed_colour   <- "#D55E00"   # vermillion (Okabe-Ito)
.default_colour <- "#0072B2"  # blue (Okabe-Ito)

#' Qualitative palette for gene colouring (CVD-safe, Okabe-Ito extended)
#'
#' Returns a character vector of n colours, recycling if needed.
#' @param n Integer. Number of colours needed.
#' @return Character vector of hex colours.
#' @noRd
.gene_palette <- function(n) {
  if (n == 0L) return(character(0))
  # For small counts use the Okabe-Ito palette: 8 hand-tuned CVD-safe colours.
  if (n <= 8L) {
    oi <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
            "#0072B2", "#D55E00", "#CC79A7", "#000000")
    return(oi[seq_len(n)])
  }
  # For larger counts use HCL qualitative colours: evenly spaced around the hue
  # wheel at fixed chroma and lightness, so no gene visually dominates and
  # complementary hues are naturally produced at any N.
  if (requireNamespace("colorspace", quietly = TRUE)) {
    return(colorspace::qualitative_hcl(n, palette = "Dynamic"))
  }
  # Fallback if colorspace is not installed: recycle Okabe-Ito with a warning.
  warning(
    "Package 'colorspace' is needed for gene palettes with more than 8 colours.\n",
    "Install it with: install.packages('colorspace')\n",
    "Recycling the 8-colour Okabe-Ito palette in the meantime.",
    call. = FALSE
  )
  rep_len(c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
            "#0072B2", "#D55E00", "#CC79A7", "#000000"), n)
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

    # No gene colouring - use a two-colour scheme: seed vs neighbour (CVD-safe)
    return(list(
      ggraph::geom_node_point(
        ggplot2::aes(
          fill  = .data[["is_seed"]],
          size  = .data[["is_seed"]],
          shape = .data[["is_seed"]]
        )
      ),
      ggplot2::scale_fill_manual(
        values = c("FALSE" = .default_colour, "TRUE" = .seed_colour),
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

  ggraph::geom_node_point(size = node_size, colour = .default_colour)
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
#' Dispatches on \linkS4class{GeneDiffPulldown}.  Three display modes are
#' available via the \code{direction} argument:
#' \itemize{
#'   \item \code{"up"} or \code{"down"} - renders a single subgraph,
#'     highlighting the seed isoforms.
#'   \item \code{"side_by_side"} (default) - renders the two subgraphs
#'     side-by-side as a two-panel \pkg{ggplot2} / \pkg{patchwork} figure.
#'   \item \code{"overlay"} - merges the up- and down-regulated subgraphs
#'     into a single network.  Nodes and edges are coloured by membership:
#'     \emph{up-only} (orange), \emph{down-only} (sky blue), or
#'     \emph{shared} (grey).  Seed isoforms are shown with a distinct shape.
#' }
#'
#' All colour palettes are colour-blind-safe (Okabe-Ito).
#'
#' @section Gene colouring:
#' When \code{group_by_gene = TRUE} and \code{level = "isoform"}, nodes are
#' coloured by their parent gene.  In \code{"up"}, \code{"down"}, and
#' \code{"side_by_side"} modes this replaces the default uniform node colour.
#' In \code{"overlay"} mode, gene identity is mapped to node \emph{fill}
#' colour while membership (up-only / shared / down-only) is shown via the
#' node \emph{border} (stroke) colour and edge colour, so both dimensions
#' remain visible simultaneously.
#'
#' @param object    A \linkS4class{GeneDiffPulldown} object.
#' @param level     \code{"isoform"} (default) or \code{"gene"}.
#' @param group_by_gene Logical.  Only used when \code{level = "isoform"}.
#'   Colours nodes by their parent gene when \code{TRUE}.  Default
#'   \code{FALSE}.
#' @param layout    ggraph layout string.  Default \code{"stress"}.
#' @param node_size Numeric. Base node size.  Default \code{3}.
#' @param label_size Numeric. Label text size.  Default \code{2.5}.
#' @param edge_alpha Numeric in \[0, 1\].  Default \code{0.4}.
#' @param title     Character scalar or \code{NULL} (auto-generated).
#' @param ...       Additional arguments forwarded to
#'   \code{\link[ggraph]{ggraph}}.
#' @param direction Character scalar: \code{"up"}, \code{"down"},
#'   \code{"side_by_side"} (default), or \code{"overlay"}.  Must be passed
#'   as a named argument.
#'
#' @return A \code{ggplot} / \code{patchwork} object, invisibly.
#'
#' @examples
#' \dontrun{
#' # Plot both directions side-by-side (default)
#' plot_network(gene_result)
#'
#' # Side-by-side with gene colouring
#' plot_network(gene_result, group_by_gene = TRUE)
#'
#' # Overlay up and down subgraphs in a single plot
#' plot_network(gene_result, direction = "overlay")
#'
#' # Overlay with gene colouring (membership shown via border + edge colour)
#' plot_network(gene_result, direction = "overlay", group_by_gene = TRUE)
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
           level         = "isoform",
           all_edges     = FALSE,
           group_by_gene = FALSE,
           layout        = "stress",
           node_size     = 3,
           label_size    = 2.5,
           edge_alpha    = 0.4,
           title         = NULL,
           ...,
           direction     = c("side_by_side", "up", "down", "overlay")) {

    direction <- match.arg(direction)
    level     <- match.arg(level, c("isoform", "gene"))
    .check_plot_args(all_edges = FALSE, group_by_gene = group_by_gene, level = level)
    .check_plot_pkgs()

    if (direction == "side_by_side") {
      return(.plot_gene_side_by_side(object, level, group_by_gene, layout,
                                     node_size, label_size, edge_alpha,
                                     title, ...))
    }

    if (direction == "overlay") {
      return(.plot_gene_overlay(object, level, group_by_gene, layout,
                                node_size, label_size, edge_alpha,
                                title, ...))
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
      group_by_gene = group_by_gene,
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
#' @param object        GeneDiffPulldown.
#' @param level         "isoform" or "gene".
#' @param group_by_gene Logical.
#' @param layout        ggraph layout.
#' @param node_size     Numeric.
#' @param label_size    Numeric.
#' @param edge_alpha    Numeric.
#' @param title         Character or NULL.
#' @param ...           Forwarded to ggraph.
#' @return patchwork object, invisibly.
#' @noRd
.plot_gene_side_by_side <- function(object, level, group_by_gene, layout,
                                    node_size, label_size, edge_alpha,
                                    title, ...) {

  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required for direction = 'side_by_side': ",
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
      group_by_gene = group_by_gene,
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


# ==============================================================================
#  Overlay rendering - internal
# ==============================================================================

#' Merge up and down subgraphs into a single overlay graph
#'
#' Performs a graph union of the up- and down-regulated subgraphs and annotates
#' every vertex and edge with a \code{membership} attribute:
#' \code{"up_only"}, \code{"down_only"}, or \code{"shared"}.
#'
#' @param up_graph   igraph or NULL.
#' @param down_graph igraph or NULL.
#' @param up_isoforms   Character vector of up-regulated seed isoform IDs.
#' @param down_isoforms Character vector of down-regulated seed isoform IDs.
#' @param level      "isoform" or "gene".
#'
#' @return A list with components:
#'   \item{graph}{The merged igraph with \code{membership} and
#'     \code{seed_direction} vertex attributes, and a \code{membership}
#'     edge attribute.}
#'   \item{up_seeds}{Character vector of seed vertex names (post-contraction
#'     when \code{level = "gene"}).}
#'   \item{down_seeds}{Character vector of seed vertex names.}
#' @noRd
.build_overlay_graph <- function(up_graph, down_graph,
                                 up_isoforms, down_isoforms, level) {

  # Collect vertex and edge sets from each direction --------------------------
  up_vnames   <- if (!is.null(up_graph))   igraph::V(up_graph)$name   else character(0)
  down_vnames <- if (!is.null(down_graph)) igraph::V(down_graph)$name else character(0)

  up_edges <- if (!is.null(up_graph)) {
    edf <- igraph::as_data_frame(up_graph, what = "edges")
    paste(edf$from, edf$to, sep = ":::")
  } else character(0)

  down_edges <- if (!is.null(down_graph)) {
    edf <- igraph::as_data_frame(down_graph, what = "edges")
    paste(edf$from, edf$to, sep = ":::")
  } else character(0)

  # Build vertex data frame with membership -----------------------------------
  all_vnames <- union(up_vnames, down_vnames)
  if (length(all_vnames) == 0L)
    stop("Both up and down subgraphs are empty; nothing to overlay.",
         call. = FALSE)

  in_up   <- all_vnames %in% up_vnames
  in_down <- all_vnames %in% down_vnames

  node_membership <- ifelse(
    in_up & in_down, "shared",
    ifelse(in_up, "up_only", "down_only")
  )

  # Seed direction: which direction(s) seeded this vertex --------------------
  is_up_seed   <- all_vnames %in% up_isoforms
  is_down_seed <- all_vnames %in% down_isoforms
  seed_direction <- ifelse(
    is_up_seed & is_down_seed, "both_seeds",
    ifelse(is_up_seed, "up_seed",
           ifelse(is_down_seed, "down_seed", "non_seed"))
  )

  # Build edge data frame with membership ------------------------------------
  # Start from both edge sets (as_data_frame preserves from/to)
  up_edf <- if (!is.null(up_graph)) {
    igraph::as_data_frame(up_graph, what = "edges")
  } else {
    data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
  }
  down_edf <- if (!is.null(down_graph)) {
    igraph::as_data_frame(down_graph, what = "edges")
  } else {
    data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
  }

  # Canonical edge key (sort endpoints so undirected comparison works)
  .edge_key <- function(edf) {
    lo <- pmin(edf$from, edf$to)
    hi <- pmax(edf$from, edf$to)
    paste(lo, hi, sep = ":::")
  }

  up_keys   <- .edge_key(up_edf)
  down_keys <- .edge_key(down_edf)

  # Union of edges (deduplicate)
  all_edf <- rbind(up_edf[, c("from", "to")], down_edf[, c("from", "to")])
  all_keys <- c(up_keys, down_keys)
  dup      <- duplicated(all_keys)
  all_edf  <- all_edf[!dup, , drop = FALSE]
  all_keys <- all_keys[!dup]

  edge_in_up   <- all_keys %in% up_keys
  edge_in_down <- all_keys %in% down_keys
  edge_membership <- ifelse(
    edge_in_up & edge_in_down, "shared",
    ifelse(edge_in_up, "up_only", "down_only")
  )

  # Carry the gene vertex attribute through if available ---------------------
  gene_attr <- rep(NA_character_, length(all_vnames))
  if (!is.null(up_graph) && !is.null(igraph::V(up_graph)$gene)) {
    idx <- match(igraph::V(up_graph)$name, all_vnames)
    gene_attr[idx] <- igraph::V(up_graph)$gene
  }
  if (!is.null(down_graph) && !is.null(igraph::V(down_graph)$gene)) {
    idx <- match(igraph::V(down_graph)$name, all_vnames)
    # Fill only NAs (up already set above)
    still_na <- is.na(gene_attr[idx])
    gene_attr[idx[still_na]] <- igraph::V(down_graph)$gene[still_na]
  }

  # Assemble igraph ----------------------------------------------------------
  g <- igraph::graph_from_data_frame(all_edf, directed = FALSE,
                                     vertices = data.frame(
                                       name = all_vnames,
                                       stringsAsFactors = FALSE
                                     ))
  igraph::V(g)$membership     <- node_membership[match(igraph::V(g)$name, all_vnames)]
  igraph::V(g)$seed_direction <- seed_direction[match(igraph::V(g)$name, all_vnames)]
  igraph::V(g)$gene           <- gene_attr[match(igraph::V(g)$name, all_vnames)]
  igraph::E(g)$membership     <- edge_membership

  # Resolve seed names for highlighting --------------------------------------
  up_seeds   <- up_isoforms
  down_seeds <- down_isoforms

  if (level == "gene") {
    # Map isoform seeds to their gene names
    if (!is.null(up_graph) && !is.null(igraph::V(up_graph)$gene)) {
      up_seeds <- unique(stats::na.omit(
        igraph::V(up_graph)$gene[igraph::V(up_graph)$name %in% up_isoforms]
      ))
    }
    if (!is.null(down_graph) && !is.null(igraph::V(down_graph)$gene)) {
      down_seeds <- unique(stats::na.omit(
        igraph::V(down_graph)$gene[igraph::V(down_graph)$name %in% down_isoforms]
      ))
    }
  }

  list(graph = g, up_seeds = up_seeds, down_seeds = down_seeds)
}


#' Render the overlay plot for a GeneDiffPulldown
#'
#' Merges up- and down-regulated subgraphs and renders them as a single
#' network with colour-coded membership.  When \code{group_by_gene = TRUE}
#' and \code{level = "isoform"}, node fill encodes gene identity while
#' the node border (stroke) and edge colour encode membership.
#'
#' @param object        GeneDiffPulldown.
#' @param level         "isoform" or "gene".
#' @param group_by_gene Logical.
#' @param layout        ggraph layout.
#' @param node_size     Numeric.
#' @param label_size    Numeric.
#' @param edge_alpha    Numeric.
#' @param title         Character or NULL.
#' @param ...           Forwarded to ggraph.
#' @return A ggplot object, invisibly.
#' @noRd
.plot_gene_overlay <- function(object, level, group_by_gene, layout,
                               node_size, label_size, edge_alpha,
                               title, ...) {

  # Build the merged graph with membership annotations
  overlay <- .build_overlay_graph(
    up_graph      = object@up_graph,
    down_graph    = object@down_graph,
    up_isoforms   = object@up_isoforms,
    down_isoforms = object@down_isoforms,
    level         = level
  )

  g <- overlay$graph

  # Apply gene-level contraction if requested
  if (level == "gene") {
    # Store membership before contraction so we can re-derive it
    pre_membership <- stats::setNames(igraph::V(g)$membership, igraph::V(g)$name)
    pre_seed_dir   <- stats::setNames(igraph::V(g)$seed_direction, igraph::V(g)$name)
    pre_edge_memb  <- igraph::E(g)$membership

    g <- .contract_to_gene(g)

    # Re-apply membership: a gene node is "shared" if any of its collapsed
    # isoforms carried different memberships, otherwise keeps the uniform value
    # For simplicity, use the overlay seeds to re-derive
    gv_names <- igraph::V(g)$name
    is_up_seed   <- gv_names %in% overlay$up_seeds
    is_down_seed <- gv_names %in% overlay$down_seeds
    igraph::V(g)$seed_direction <- ifelse(
      is_up_seed & is_down_seed, "both_seeds",
      ifelse(is_up_seed, "up_seed",
             ifelse(is_down_seed, "down_seed", "non_seed"))
    )

    # For gene-level membership, re-derive from the original node memberships
    # by checking which genes appeared in up and/or down graphs
    up_genes <- if (!is.null(object@up_graph) && !is.null(igraph::V(object@up_graph)$gene)) {
      unique(stats::na.omit(igraph::V(object@up_graph)$gene))
    } else character(0)
    down_genes <- if (!is.null(object@down_graph) && !is.null(igraph::V(object@down_graph)$gene)) {
      unique(stats::na.omit(igraph::V(object@down_graph)$gene))
    } else character(0)

    in_up   <- gv_names %in% up_genes
    in_down <- gv_names %in% down_genes
    igraph::V(g)$membership <- ifelse(
      in_up & in_down, "shared",
      ifelse(in_up, "up_only", "down_only")
    )

    # Re-derive edge membership at gene level
    edf <- igraph::as_data_frame(g, what = "edges")
    if (nrow(edf) > 0) {
      # An edge is "shared" if it existed in both direction graphs at gene level
      up_gene_edges <- if (!is.null(object@up_graph)) {
        ue <- igraph::as_data_frame(object@up_graph, what = "edges")
        if (nrow(ue) > 0 && !is.null(igraph::V(object@up_graph)$gene)) {
          iso2gene <- stats::setNames(igraph::V(object@up_graph)$gene,
                                      igraph::V(object@up_graph)$name)
          lo <- pmin(iso2gene[ue$from], iso2gene[ue$to])
          hi <- pmax(iso2gene[ue$from], iso2gene[ue$to])
          unique(paste(lo, hi, sep = ":::"))
        } else character(0)
      } else character(0)

      down_gene_edges <- if (!is.null(object@down_graph)) {
        de <- igraph::as_data_frame(object@down_graph, what = "edges")
        if (nrow(de) > 0 && !is.null(igraph::V(object@down_graph)$gene)) {
          iso2gene <- stats::setNames(igraph::V(object@down_graph)$gene,
                                      igraph::V(object@down_graph)$name)
          lo <- pmin(iso2gene[de$from], iso2gene[de$to])
          hi <- pmax(iso2gene[de$from], iso2gene[de$to])
          unique(paste(lo, hi, sep = ":::"))
        } else character(0)
      } else character(0)

      merged_keys <- paste(pmin(edf$from, edf$to), pmax(edf$from, edf$to),
                           sep = ":::")
      e_in_up   <- merged_keys %in% up_gene_edges
      e_in_down <- merged_keys %in% down_gene_edges
      igraph::E(g)$membership <- ifelse(
        e_in_up & e_in_down, "shared",
        ifelse(e_in_up, "up_only", "down_only")
      )
    }
  }

  # Remove isolated vertices
  g <- igraph::delete.vertices(g, igraph::degree(g) == 0L)

  if (igraph::vcount(g) == 0L)
    stop("Overlay graph has no vertices after cleanup.", call. = FALSE)

  # -- Colour palettes (CVD-safe, Okabe-Ito) ----------------------------------
  membership_colours <- .membership_palette
  edge_colours       <- .edge_membership_palette

  seed_shapes <- c(
    non_seed   = 16L,   # filled circle
    up_seed    = 17L,   # triangle (pointing up)
    down_seed  = 25L,   # inverted triangle (pointing down)
    both_seeds = 18L    # diamond
  )

  # Shapes that support fill + stroke (for group_by_gene overlay)
  seed_shapes_filled <- c(
    non_seed   = 21L,   # filled circle
    up_seed    = 24L,   # filled triangle up
    down_seed  = 25L,   # filled triangle down
    both_seeds = 23L    # filled diamond
  )

  # -- Build ggraph plot -------------------------------------------------------
  if (is.null(title))
    title <- sprintf("Differential pulldown overlay  -  %s  (%s level)",
                     object@gene_id, level)

  subtitle <- sprintf(
    "%d nodes  \u00b7  %d edges    (up-only: %d  |  shared: %d  |  down-only: %d)",
    igraph::vcount(g), igraph::ecount(g),
    sum(igraph::V(g)$membership == "up_only"),
    sum(igraph::V(g)$membership == "shared"),
    sum(igraph::V(g)$membership == "down_only")
  )

  # Ensure factor levels for consistent legend ordering
  igraph::V(g)$membership     <- factor(igraph::V(g)$membership,
                                        levels = c("up_only", "shared", "down_only"))
  igraph::V(g)$seed_direction <- factor(igraph::V(g)$seed_direction,
                                        levels = c("non_seed", "up_seed",
                                                   "down_seed", "both_seeds"))
  igraph::E(g)$membership     <- factor(igraph::E(g)$membership,
                                        levels = c("up_only", "shared", "down_only"))

  # Scale seed nodes slightly larger
  size_vals <- stats::setNames(
    c(node_size, node_size * 1.5, node_size * 1.5, node_size * 1.5),
    c("non_seed", "up_seed", "down_seed", "both_seeds")
  )

  membership_labels <- c(up_only = "Up only", shared = "Shared",
                         down_only = "Down only")
  seed_labels       <- c(non_seed = "Neighbour", up_seed = "Up seed",
                         down_seed = "Down seed", both_seeds = "Both seeds")

  # ---- Branch: group_by_gene on overlay ------------------------------------
  # When group_by_gene is TRUE at isoform level, we use:
  #   fill  = gene   (colour identity)
  #   stroke/colour = membership  (border colour encodes up/down/shared)
  #   shape = seed_direction
  #   edge colour = membership
  use_gene_fill <- group_by_gene && level == "isoform" &&
    !is.null(igraph::V(g)$gene)

  if (use_gene_fill) {
    gene_ids <- igraph::V(g)$gene
    gene_ids[is.na(gene_ids)] <- "NA_gene"
    igraph::V(g)$gene <- gene_ids
    unique_genes <- unique(gene_ids)
    gene_cols    <- .gene_palette(length(unique_genes))
    names(gene_cols) <- unique_genes

    p <- ggraph::ggraph(g, layout = layout, ...) +
      ggraph::geom_edge_link(
        ggplot2::aes(colour = .data[["membership"]]),
        alpha = edge_alpha
      ) +
      ggraph::scale_edge_colour_manual(
        name   = "Edge membership",
        values = edge_colours,
        labels = membership_labels,
        drop   = FALSE
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(
          fill   = .data[["gene"]],
          colour = .data[["membership"]],
          shape  = .data[["seed_direction"]],
          size   = .data[["seed_direction"]]
        ),
        stroke = 1.5
      ) +
      ggplot2::scale_fill_manual(
        name   = "Gene",
        values = gene_cols
      ) +
      ggplot2::scale_colour_manual(
        name   = "Membership (border)",
        values = membership_colours,
        labels = membership_labels,
        drop   = FALSE
      ) +
      ggplot2::scale_shape_manual(
        name   = "Seed status",
        values = seed_shapes_filled,
        labels = seed_labels,
        drop   = FALSE
      ) +
      ggplot2::scale_size_manual(
        name   = "Seed status",
        values = size_vals,
        labels = seed_labels,
        drop   = FALSE
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label = .data[["name"]]),
        size  = label_size,
        repel = TRUE
      ) +
      ggplot2::labs(title = title, subtitle = subtitle) +
      ggraph::theme_graph(base_family = "sans") +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 13, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40"),
        legend.position = "right"
      ) +
      ggplot2::guides(
        fill        = ggplot2::guide_legend(title = "Gene",
                                            override.aes = list(shape = 21,
                                                                size = 3,
                                                                stroke = 0),
                                            order = 1,
                                            ncol  = 1),
        colour      = ggplot2::guide_legend(title = "Membership (border)",
                                            override.aes = list(shape = 21,
                                                                size = 3,
                                                                fill  = "white"),
                                            order = 2),
        edge_colour = ggplot2::guide_legend(title = "Edge", order = 3),
        shape       = ggplot2::guide_legend(title = "Seed status", order = 4),
        size        = ggplot2::guide_legend(title = "Seed status", order = 4)
      )

  } else {

    # ---- Default overlay: colour = membership -------------------------------
    p <- ggraph::ggraph(g, layout = layout, ...) +
      ggraph::geom_edge_link(
        ggplot2::aes(colour = .data[["membership"]]),
        alpha = edge_alpha
      ) +
      ggraph::scale_edge_colour_manual(
        name   = "Edge membership",
        values = edge_colours,
        labels = membership_labels,
        drop   = FALSE
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(
          colour = .data[["membership"]],
          shape  = .data[["seed_direction"]],
          size   = .data[["seed_direction"]]
        )
      ) +
      ggplot2::scale_colour_manual(
        name   = "Node membership",
        values = membership_colours,
        labels = membership_labels,
        drop   = FALSE
      ) +
      ggplot2::scale_shape_manual(
        name   = "Seed status",
        values = seed_shapes,
        labels = seed_labels,
        drop   = FALSE
      ) +
      ggplot2::scale_size_manual(
        name   = "Seed status",
        values = size_vals,
        labels = seed_labels,
        drop   = FALSE
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label = .data[["name"]]),
        size  = label_size,
        repel = TRUE
      ) +
      ggplot2::labs(title = title, subtitle = subtitle) +
      ggraph::theme_graph(base_family = "sans") +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 13, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40"),
        legend.position = "right"
      ) +
      ggplot2::guides(
        colour      = ggplot2::guide_legend(title = "Membership", order = 1),
        edge_colour = ggplot2::guide_legend(title = "Edge", order = 2),
        shape       = ggplot2::guide_legend(title = "Seed status", order = 3),
        size        = ggplot2::guide_legend(title = "Seed status", order = 3)
      )
  }

  print(p)
  invisible(p)
}


# ==============================================================================
#  plot_network - DifferentialPulldownAnalysis (cross-gene view)
# ==============================================================================

#' Plot the combined network across all (or selected) genes
#'
#' Dispatches on \linkS4class{DifferentialPulldownAnalysis}.  Merges the
#' virtual pulldown subgraphs from multiple genes into a single network to
#' reveal cross-gene clustering and shared interaction partners.
#'
#' Nodes are coloured by the gene whose pulldown they originated from.
#' Nodes that appear in the subgraphs of more than one gene receive a
#' distinct \emph{multi-gene} colour (dark grey) so that shared interaction
#' hubs are immediately visible.  Seed isoforms are highlighted via shape.
#'
#' All colour palettes are colour-blind-safe (Okabe-Ito).
#'
#' @param object    A \linkS4class{DifferentialPulldownAnalysis} object.
#' @param level     \code{"isoform"} (default) or \code{"gene"}.
#' @param layout    ggraph layout string.  Default \code{"stress"}.
#' @param node_size Numeric. Base node size.  Default \code{3}.
#' @param label_size Numeric. Label text size.  Default \code{2.5}.
#' @param edge_alpha Numeric in \[0, 1\].  Default \code{0.4}.
#' @param title     Character scalar or \code{NULL} (auto-generated).
#' @param ...       Additional arguments forwarded to
#'   \code{\link[ggraph]{ggraph}}.
#' @param direction Character scalar: \code{"up"} (merge all up-regulated
#'   subgraphs), \code{"down"} (merge all down-regulated subgraphs), or
#'   \code{"both"} (default; union of up and down per gene).  Must be
#'   passed as a named argument.
#' @param genes     Character vector of gene IDs to include.  Defaults to
#'   \code{NULL}, which includes all genes in the analysis.  Must be
#'   passed as a named argument.
#' @param min_nodes Integer or \code{NULL} (default).  Pre-merge per-gene
#'   filter: genes whose merged subgraph (union of up + down vertices for the
#'   requested \code{direction}) contains fewer than \code{min_nodes} nodes are
#'   dropped before the cross-gene graph is built.  Useful for removing genes
#'   with very small or absent pulldown neighborhoods.  Set to \code{NULL} to
#'   disable.  Must be passed as a named argument.
#' @param min_shared_nodes Integer or \code{NULL} (default).  Post-merge
#'   component filter: connected components of the merged graph that contain
#'   fewer than \code{min_shared_nodes} multi-gene nodes (nodes attributed to
#'   more than one gene's pulldown) are removed.  This is the strongest lever
#'   for keeping only biologically interesting cross-gene hubs.  Set to
#'   \code{NULL} to disable.  Must be passed as a named argument.
#' @param min_component_nodes Integer or \code{NULL} (default).  Post-merge
#'   component filter: connected components with fewer than
#'   \code{min_component_nodes} total nodes are removed after merging.
#'   Complements \code{min_shared_nodes} by also dropping small single-gene
#'   islands.  Set to \code{NULL} to disable.  Must be passed as a named
#'   argument.
#' @param min_edges Integer or \code{NULL} (default).  Post-merge component
#'   filter: connected components with fewer than \code{min_edges} edges are
#'   removed.  Helps discard loosely connected star-shaped or chain components
#'   that pass a node-count threshold.  Set to \code{NULL} to disable.  Must
#'   be passed as a named argument.
#' @param colour_by Character scalar: \code{"gene"} (default) or
#'   \code{"hub_degree"}.  Controls the node colour encoding.
#'   \describe{
#'     \item{\code{"gene"}}{Each node is coloured by the gene whose pulldown it
#'       originated from (Okabe-Ito for N <= 8, HCL Dynamic for larger N).
#'       Multi-gene shared nodes are shown in dark grey; their size is scaled
#'       by degree so high-connectivity hubs stand out.}
#'     \item{\code{"hub_degree"}}{Non-shared nodes are collapsed to light grey.
#'       Shared multi-gene nodes are coloured by degree on a sequential Mako
#'       scale (dark = low degree, cream = high degree), making cross-gene
#'       interaction hubs the visual focus of the plot.}
#'   }
#'   Must be passed as a named argument.
#'
#' @return A \code{ggplot} object, invisibly.
#'
#' @examples
#' \dontrun{
#' # After running a differential pulldown:
#' obj <- run_differential_pulldown(obj, "my_run", de_table, order = 1)
#' result <- obj@@analyses$my_run
#'
#' # Cross-gene view of up-regulated subgraphs
#' plot_network(result, direction = "up")
#'
#' # Both directions, subset to specific genes
#' plot_network(result, genes = c("BRCA1", "TP53"))
#'
#' # Gene-level view
#' plot_network(result, level = "gene")
#'
#' # Keep only components with >= 2 multi-gene nodes (shared hubs)
#' plot_network(result, min_shared_nodes = 2)
#'
#' # Drop small components and loosely connected ones
#' plot_network(result, min_component_nodes = 5, min_edges = 4)
#'
#' # Pre-filter genes with tiny neighborhoods, then keep shared-hub components
#' plot_network(result, min_nodes = 10, min_shared_nodes = 2)
#'
#' # Hub-degree colour mode: shared nodes coloured by connectivity
#' plot_network(result, colour_by = "hub_degree", min_shared_nodes = 2)
#' }
#'
#' @describeIn plot_network Plot cross-gene network from a
#'   \linkS4class{DifferentialPulldownAnalysis} result.
setMethod(
  "plot_network",
  "DifferentialPulldownAnalysis",
  function(object,
           level      = "isoform",
           all_edges  = FALSE,
           group_by_gene = FALSE,
           layout     = "stress",
           node_size  = 3,
           label_size = 2.5,
           edge_alpha = 0.4,
           title      = NULL,
           ...,
           direction          = c("both", "up", "down"),
           genes              = NULL,
           min_nodes          = NULL,
           min_shared_nodes   = NULL,
           min_component_nodes = NULL,
           min_edges          = NULL,
           colour_by          = c("gene", "hub_degree")) {

    direction <- match.arg(direction)
    level     <- match.arg(level, c("isoform", "gene"))
    colour_by <- match.arg(colour_by)
    .check_plot_pkgs()

    # -- Validate filter arguments ---------------------------------------------
    .check_filter_arg <- function(x, name) {
      if (!is.null(x) && (!is.numeric(x) || length(x) != 1L || x < 1L))
        stop(sprintf("`%s` must be a single positive integer or NULL.", name),
             call. = FALSE)
    }
    .check_filter_arg(min_nodes,           "min_nodes")
    .check_filter_arg(min_shared_nodes,    "min_shared_nodes")
    .check_filter_arg(min_component_nodes, "min_component_nodes")
    .check_filter_arg(min_edges,           "min_edges")

    # -- Resolve gene subset ---------------------------------------------------
    all_gene_ids <- names(object@gene_results)

    if (is.null(genes)) {
      genes <- all_gene_ids
    } else {
      bad <- setdiff(genes, all_gene_ids)
      if (length(bad))
        stop(sprintf(
          "Gene(s) not found in this analysis: %s. Available: %s",
          paste(bad, collapse = ", "),
          paste(all_gene_ids, collapse = ", ")
        ), call. = FALSE)
    }

    # -- Stage 1: per-gene pre-filter (min_nodes) ------------------------------
    # Count how many unique nodes each gene contributes for the requested
    # direction before merging, then drop genes that fall below the threshold.
    if (!is.null(min_nodes)) {
      .get_gene_node_count <- function(gr, direction) {
        vnames <- character(0)
        if (direction %in% c("up",   "both") && !is.null(gr@up_graph))
          vnames <- union(vnames, igraph::V(gr@up_graph)$name)
        if (direction %in% c("down", "both") && !is.null(gr@down_graph))
          vnames <- union(vnames, igraph::V(gr@down_graph)$name)
        length(vnames)
      }
      n_nodes_per_gene <- vapply(
        object@gene_results[genes],
        .get_gene_node_count,
        direction = direction,
        FUN.VALUE = integer(1L)
      )
      dropped_pre <- genes[n_nodes_per_gene < min_nodes]
      genes       <- genes[n_nodes_per_gene >= min_nodes]
      if (length(dropped_pre) > 0L)
        message(sprintf(
          "min_nodes filter: removed %d gene(s) with fewer than %d nodes: %s",
          length(dropped_pre), min_nodes,
          paste(dropped_pre, collapse = ", ")
        ))
      if (length(genes) == 0L)
        stop("All genes were removed by the min_nodes filter. ",
             "Lower the threshold or set min_nodes = NULL.", call. = FALSE)
    }

    if (length(genes) < 2L)
      warning(
        "Cross-gene plot with fewer than 2 genes \u2014 ",
        "consider using the GeneDiffPulldown plot_network() method instead.",
        call. = FALSE
      )

    # -- Build the merged cross-gene graph ------------------------------------
    built <- .build_cross_gene_graph(
      gene_results = object@gene_results[genes],
      direction    = direction,
      level        = level
    )

    g          <- built$graph
    seed_names <- built$seed_names

    if (igraph::vcount(g) == 0L)
      stop("No vertices in the merged graph. ",
           "Check that the selected genes and direction have non-empty subgraphs.",
           call. = FALSE)

    # -- Stage 2: post-merge component filters ---------------------------------
    if (!is.null(min_shared_nodes) || !is.null(min_component_nodes) ||
        !is.null(min_edges)) {
      g <- .filter_cross_gene_components(
        g                   = g,
        min_shared_nodes    = min_shared_nodes,
        min_component_nodes = min_component_nodes,
        min_edges           = min_edges
      )
      if (igraph::vcount(g) == 0L)
        stop(
          "No components survived the post-merge filters ",
          "(min_shared_nodes / min_component_nodes / min_edges). ",
          "Try relaxing the thresholds.",
          call. = FALSE
        )
      # Restrict seed_names to vertices that survived filtering
      seed_names <- seed_names[seed_names %in% igraph::V(g)$name]
    }

    # -- Render ---------------------------------------------------------------
    dir_label <- switch(direction,
                        up   = "up-regulated",
                        down = "down-regulated",
                        both = "up + down")

    if (is.null(title))
      title <- sprintf(
        "Cross-gene network  \u2014  %s  (%d gene(s), %s, %s level)",
        object@label, length(genes), dir_label, level
      )

    .render_cross_gene_plot(
      g          = g,
      seed_names = seed_names,
      level      = level,
      layout     = layout,
      node_size  = node_size,
      label_size = label_size,
      edge_alpha = edge_alpha,
      title      = title,
      colour_by  = colour_by,
      ...
    )
  }
)


# ==============================================================================
#  Internal: post-merge component filter
# ==============================================================================

#' Filter connected components of the merged cross-gene graph
#'
#' Decomposes the merged graph into connected components and removes any
#' component that fails to meet the requested thresholds.  All three criteria
#' are applied in combination: a component must pass \emph{all} active
#' (non-NULL) filters to be retained.
#'
#' @param g                   igraph.  The merged cross-gene graph with a
#'   \code{source_gene} vertex attribute.
#' @param min_shared_nodes    Integer or NULL.  Minimum number of multi-gene
#'   nodes (vertices with \code{source_gene == "multi-gene"}) a component must
#'   contain.
#' @param min_component_nodes Integer or NULL.  Minimum total node count per
#'   component.
#' @param min_edges           Integer or NULL.  Minimum edge count per
#'   component.
#'
#' @return The filtered igraph (isolated vertices already removed upstream, so
#'   every retained vertex has degree >= 1).
#' @noRd
.filter_cross_gene_components <- function(g,
                                          min_shared_nodes    = NULL,
                                          min_component_nodes = NULL,
                                          min_edges           = NULL) {

  comps      <- igraph::decompose(g)
  n_before   <- length(comps)
  keep_flags <- logical(n_before)

  for (i in seq_along(comps)) {
    cg <- comps[[i]]

    # -- min_component_nodes ---------------------------------------------------
    if (!is.null(min_component_nodes) &&
        igraph::vcount(cg) < min_component_nodes) {
      next
    }

    # -- min_edges -------------------------------------------------------------
    if (!is.null(min_edges) && igraph::ecount(cg) < min_edges) {
      next
    }

    # -- min_shared_nodes ------------------------------------------------------
    if (!is.null(min_shared_nodes)) {
      n_multi <- sum(igraph::V(cg)$source_gene == "multi-gene", na.rm = TRUE)
      if (n_multi < min_shared_nodes) next
    }

    keep_flags[i] <- TRUE
  }

  n_kept    <- sum(keep_flags)
  n_dropped <- n_before - n_kept

  if (n_dropped > 0L)
    message(sprintf(
      "Component filters removed %d of %d component(s); %d retained.",
      n_dropped, n_before, n_kept
    ))

  if (n_kept == 0L)
    return(igraph::make_empty_graph(directed = FALSE))

  # Collect the vertex names that belong to surviving components and induce
  # a subgraph from the *original* graph.  This preserves all vertex
  # attributes (source_gene, is_seed, gene, ...) exactly as they were.
  # Using igraph::union() on sub-graphs would drop or mangle attributes.
  kept_vnames <- unlist(
    lapply(comps[keep_flags], function(cg) igraph::V(cg)$name),
    use.names = FALSE
  )
  igraph::induced_subgraph(g, vids = kept_vnames)
}


# ==============================================================================
#  Internal: build the merged cross-gene graph
# ==============================================================================

#' Merge subgraphs across genes into a single annotated graph
#'
#' For each gene, extracts the requested direction subgraph(s), collects
#' vertices and edges, and annotates every vertex with the gene(s) whose
#' pulldown it belongs to.  Vertices appearing in more than one gene's
#' subgraph are labelled \code{"multi-gene"}.
#'
#' @param gene_results Named list of GeneDiffPulldown objects.
#' @param direction    \code{"up"}, \code{"down"}, or \code{"both"}.
#' @param level        \code{"isoform"} or \code{"gene"}.
#'
#' @return A list with:
#'   \item{graph}{igraph with vertex attributes \code{source_gene}
#'     (gene ID or \code{"multi-gene"}) and \code{is_seed} (logical).}
#'   \item{seed_names}{Character vector of all seed vertex names.}
#' @noRd
.build_cross_gene_graph <- function(gene_results, direction, level) {

  # Collect per-gene subgraphs -------------------------------------------------
  .get_graphs <- function(gr, direction) {
    graphs <- list()
    if (direction %in% c("up", "both") && !is.null(gr@up_graph))
      graphs <- c(graphs, list(gr@up_graph))
    if (direction %in% c("down", "both") && !is.null(gr@down_graph))
      graphs <- c(graphs, list(gr@down_graph))
    graphs
  }

  # For each gene, gather: vertex names, edges, and seed isoform names
  per_gene_vnames <- list()   # gene_id -> character vector of vertex names
  all_edge_dfs    <- list()
  all_seed_names  <- character(0)

  for (gid in names(gene_results)) {
    gr      <- gene_results[[gid]]
    graphs  <- .get_graphs(gr, direction)

    gene_vnames <- character(0)
    for (g in graphs) {
      gene_vnames <- union(gene_vnames, igraph::V(g)$name)

      edf <- igraph::as_data_frame(g, what = "edges")[, c("from", "to"),
                                                      drop = FALSE]
      if (nrow(edf) > 0)
        all_edge_dfs <- c(all_edge_dfs, list(edf))
    }

    per_gene_vnames[[gid]] <- gene_vnames

    # Collect seed isoforms for this gene (depending on direction)
    if (direction %in% c("up", "both"))
      all_seed_names <- union(all_seed_names, gr@up_isoforms)
    if (direction %in% c("down", "both"))
      all_seed_names <- union(all_seed_names, gr@down_isoforms)
  }

  # Determine source_gene per vertex ------------------------------------------
  all_vnames <- unique(unlist(per_gene_vnames, use.names = FALSE))

  if (length(all_vnames) == 0L)
    return(list(
      graph      = igraph::make_empty_graph(directed = FALSE),
      seed_names = character(0)
    ))

  # Map each vertex to the gene(s) it belongs to
  vertex_gene_list <- stats::setNames(
    vector("list", length(all_vnames)),
    all_vnames
  )

  for (gid in names(per_gene_vnames)) {
    for (v in per_gene_vnames[[gid]]) {
      vertex_gene_list[[v]] <- c(vertex_gene_list[[v]], gid)
    }
  }

  source_gene <- vapply(vertex_gene_list, function(gids) {
    gids <- unique(gids)
    if (length(gids) == 1L) gids else "multi-gene"
  }, character(1))

  # Carry the gene (parent gene) vertex attribute through if available --------
  gene_attr <- rep(NA_character_, length(all_vnames))
  names(gene_attr) <- all_vnames

  for (gid in names(gene_results)) {
    gr     <- gene_results[[gid]]
    graphs <- .get_graphs(gr, direction)
    for (g in graphs) {
      g_gene <- igraph::V(g)$gene
      if (!is.null(g_gene)) {
        g_names <- igraph::V(g)$name
        still_na <- is.na(gene_attr[g_names])
        gene_attr[g_names[still_na]] <- g_gene[still_na]
      }
    }
  }

  # Build union edge data frame (deduplicated) --------------------------------
  if (length(all_edge_dfs) > 0) {
    combined_edf <- do.call(rbind, all_edge_dfs)
    # Canonical key for deduplication (undirected)
    ekey <- paste(pmin(combined_edf$from, combined_edf$to),
                  pmax(combined_edf$from, combined_edf$to),
                  sep = ":::")
    combined_edf <- combined_edf[!duplicated(ekey), , drop = FALSE]
  } else {
    combined_edf <- data.frame(from = character(0), to = character(0),
                               stringsAsFactors = FALSE)
  }

  # Assemble igraph -----------------------------------------------------------
  g <- igraph::graph_from_data_frame(
    combined_edf,
    directed = FALSE,
    vertices = data.frame(name = all_vnames, stringsAsFactors = FALSE)
  )

  igraph::V(g)$source_gene <- source_gene[igraph::V(g)$name]
  igraph::V(g)$gene        <- gene_attr[igraph::V(g)$name]
  igraph::V(g)$is_seed     <- igraph::V(g)$name %in% all_seed_names

  # Gene-level contraction if requested ----------------------------------------
  if (level == "gene") {
    pre_source <- stats::setNames(igraph::V(g)$source_gene, igraph::V(g)$name)

    g <- .contract_to_gene(g)

    # Re-derive source_gene at gene level: a gene node is "multi-gene" if
    # isoforms from different pulldown genes contributed to it
    gv_names <- igraph::V(g)$name

    gene_source <- vapply(gv_names, function(gn) {
      # Which pulldown genes contributed isoforms for this parent gene?
      contributing <- names(per_gene_vnames)[
        vapply(per_gene_vnames, function(vn) {
          any(gene_attr[vn] == gn, na.rm = TRUE)
        }, logical(1))
      ]
      contributing <- unique(contributing)
      if (length(contributing) == 1L) contributing else "multi-gene"
    }, character(1))

    igraph::V(g)$source_gene <- gene_source

    # Re-derive is_seed at gene level
    seed_genes <- unique(stats::na.omit(
      gene_attr[all_seed_names[all_seed_names %in% names(gene_attr)]]
    ))
    igraph::V(g)$is_seed <- gv_names %in% seed_genes
  }

  # Remove isolated vertices
  g <- igraph::delete.vertices(g, igraph::degree(g) == 0L)

  list(
    graph      = g,
    seed_names = if (level == "gene") {
      igraph::V(g)$name[igraph::V(g)$is_seed]
    } else {
      all_seed_names
    }
  )
}


# ==============================================================================
#  Internal: render the cross-gene plot
# ==============================================================================

#' Render the cross-gene overlay plot
#'
#' @param g          igraph with \code{source_gene} and \code{is_seed}
#'   vertex attributes.
#' @param seed_names Character vector of seed vertex names.
#' @param level      "isoform" or "gene".
#' @param layout     ggraph layout string.
#' @param node_size  Numeric.
#' @param label_size Numeric.
#' @param edge_alpha Numeric.
#' @param title      Character scalar.
#' @param ...        Forwarded to ggraph.
#' @return A ggplot object, invisibly.
#' @noRd
.render_cross_gene_plot <- function(g, seed_names, level, layout,
                                    node_size, label_size, edge_alpha,
                                    title, colour_by = "gene", ...) {

  source_genes <- igraph::V(g)$source_gene
  is_multi     <- source_genes == "multi-gene"
  n_multi      <- sum(is_multi)
  deg          <- igraph::degree(g)

  # -- Subtitle ----------------------------------------------------------------
  subtitle <- sprintf(
    "%d nodes  \u00b7  %d edges  \u00b7  %d multi-gene node(s)",
    igraph::vcount(g), igraph::ecount(g), n_multi
  )

  # -- Seed shapes (shared by both modes) --------------------------------------
  igraph::V(g)$is_seed <- factor(
    as.character(igraph::V(g)$name %in% seed_names),
    levels = c("FALSE", "TRUE")
  )
  shape_vals <- c("FALSE" = 16L, "TRUE" = 18L)   # circle vs diamond
  size_vals  <- c("FALSE" = node_size, "TRUE" = node_size * 1.6)

  # ============================================================================
  #  MODE A: colour_by = "gene"
  #  Gene colours for single-gene nodes; grey for shared nodes.
  #  Fixed node sizes — seeds slightly larger.
  # ============================================================================
  if (colour_by == "gene") {

    unique_sources <- unique(source_genes)
    real_genes     <- sort(setdiff(unique_sources, "multi-gene"))
    has_multi      <- "multi-gene" %in% unique_sources

    gene_cols <- .gene_palette(length(real_genes))
    names(gene_cols) <- real_genes
    if (has_multi)
      gene_cols <- c(gene_cols, "multi-gene" = "#444444")

    source_levels <- c(real_genes, if (has_multi) "multi-gene")
    igraph::V(g)$source_gene <- factor(source_genes, levels = source_levels)

    p <- ggraph::ggraph(g, layout = layout, ...) +
      ggraph::geom_edge_link(alpha = edge_alpha, colour = "grey60") +
      ggraph::geom_node_point(
        ggplot2::aes(
          colour = .data[["source_gene"]],
          shape  = .data[["is_seed"]],
          size   = .data[["is_seed"]]
        )
      ) +
      ggplot2::scale_colour_manual(
        name   = "Source gene",
        values = gene_cols,
        drop   = FALSE
      ) +
      ggplot2::scale_shape_manual(
        name   = "Role",
        values = shape_vals,
        labels = c("FALSE" = "Neighbour", "TRUE" = "Seed"),
        drop   = FALSE
      ) +
      ggplot2::scale_size_manual(
        name   = "Role",
        values = size_vals,
        labels = c("FALSE" = "Neighbour", "TRUE" = "Seed"),
        drop   = FALSE
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label = .data[["name"]]),
        size  = label_size,
        repel = TRUE
      ) +
      ggplot2::labs(title = title, subtitle = subtitle) +
      ggraph::theme_graph(base_family = "sans") +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 13, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 9,  colour = "grey40"),
        legend.position = "right"
      ) +
      ggplot2::guides(
        colour = "none",
        shape  = ggplot2::guide_legend(title = "Role", order = 1),
        size   = ggplot2::guide_legend(title = "Role", order = 1)
      )

    print(p)
    return(invisible(p))
  }

  # ============================================================================
  #  MODE B: colour_by = "hub_degree"
  #  Non-shared nodes → light grey.  Shared nodes → Mako sequential scale
  #  by degree (dark = low connectivity, cream = high connectivity).
  #  Fixed node sizes — seeds slightly larger.
  # ============================================================================

  # Set degree to NA for non-shared nodes so they map to the na.value grey
  igraph::V(g)$hub_degree <- ifelse(is_multi, deg, NA_real_)

  # Mako palette via colorspace (falls back to viridis if unavailable)
  if (requireNamespace("colorspace", quietly = TRUE)) {
    scale_col <- colorspace::scale_colour_continuous_sequential(
      palette  = "Mako",
      name     = "Degree\n(shared nodes)",
      na.value = "#CCCCCC",
      guide    = ggplot2::guide_colorbar(order = 2)
    )
  } else {
    warning("colorspace not installed; falling back to viridis for hub_degree mode.",
            call. = FALSE)
    scale_col <- ggplot2::scale_colour_viridis_c(
      name     = "Degree\n(shared nodes)",
      na.value = "#CCCCCC",
      option   = "D"
    )
  }

  p <- ggraph::ggraph(g, layout = layout, ...) +
    ggraph::geom_edge_link(alpha = edge_alpha, colour = "grey75") +
    ggraph::geom_node_point(
      ggplot2::aes(
        colour = .data[["hub_degree"]],
        shape  = .data[["is_seed"]],
        size   = .data[["is_seed"]]
      )
    ) +
    scale_col +
    ggplot2::scale_shape_manual(
      name   = "Role",
      values = shape_vals,
      labels = c("FALSE" = "Neighbour", "TRUE" = "Seed"),
      drop   = FALSE
    ) +
    ggplot2::scale_size_manual(
      name   = "Role",
      values = size_vals,
      labels = c("FALSE" = "Neighbour", "TRUE" = "Seed"),
      drop   = FALSE
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data[["name"]]),
      size  = label_size,
      repel = TRUE
    ) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9,  colour = "grey40"),
      legend.position = "right"
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(title = "Role", order = 1),
      size  = ggplot2::guide_legend(title = "Role", order = 1)
    )

  print(p)
  invisible(p)
}
