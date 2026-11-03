## plot_network.R
## Generic and methods for plotting isoform interaction networks.
## Dispatches on IsoNetAnalyzeR and VirtualPulldownLookup; shared rendering
## logic lives in internal helpers that accept plain igraph objects.
##
## Exported: plot_network

#' @include isonetanalyzer_class.R
#' @include lookup_class.R
#' @include diff_pulldown_gene_class.R
#' @include diff_pulldown_analysis_class.R
#' @include graph_utils.R
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
#'   \code{ggraph::create_layout()}.  Defaults to \code{"stress"}, which
#'   is faster and generally better than \code{"fr"} for biological networks.
#' @param node_size Numeric scalar.  Base node size.  Default \code{3}.
#' @param label_size Numeric scalar.  Node label text size.  Default \code{2.5}.
#' @param edge_alpha Numeric in \[0, 1\].  Edge opacity.  Default \code{0.4}.
#' @param title Character scalar.  Plot title.  Auto-generated when \code{NULL}
#'   (default).
#' @param node_labels Character scalar.  \code{"symbol"} (default) labels
#'   nodes with gene symbols, \code{"id"} with transcript (or, at gene level,
#'   gene) identifiers, and \code{"none"} draws no labels.  Symbols come from a
#'   \code{gene_name} vertex attribute when the network has one, otherwise
#'   from \code{\link{isonet_gene_names}()}; nodes without a known symbol keep
#'   their identifier.
#' @param ... Additional arguments forwarded to \code{ggraph::ggraph()}.
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
#' \code{"NA_gene"}. Isolated seed nodes are retained and rendered as valid
#' one-node/zero-edge results; isolated non-seed nodes are omitted.
#'
#' @return A \code{ggplot} object, invisibly.  Can be further modified with
#'   \pkg{ggplot2} layers or saved with \code{ggplot2::ggsave()}.
#'
#' @examples
#' obj <- IsoNetAnalyzeR("example")
#'
#' # Neighbourhood of one isoform
#' obj <- lookup_virtual_pulldown(obj, "cltcl1_up", seed = "ENST00000621271",
#'                                order = 1)
#' p <- plot_network(get_lookup(obj, "cltcl1_up"))
#'
#' # The same neighbourhood at gene level
#' p <- plot_network(get_lookup(obj, "cltcl1_up"), level = "gene")
#'
#' # Transcript identifiers instead of gene symbols, or no labels
#' p <- plot_network(get_lookup(obj, "cltcl1_up"), node_labels = "id")
#' p <- plot_network(get_lookup(obj, "cltcl1_up"), node_labels = "none")
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
           node_labels   = c("symbol", "id", "none"),
           ...)
    standardGeneric("plot_network")
)


# ==============================================================================
#  Method - IsoNetAnalyzeR
# ==============================================================================
#' @describeIn plot_network Plot the full network of an
#'   \linkS4class{IsoNetAnalyzeR} session object.
#' @param force Logical. If \code{TRUE}, bypass the size guard and attempt to
#'   plot the full network regardless of its size. Default \code{FALSE}.
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
           node_labels   = c("symbol", "id", "none"),
           force         = FALSE,
           ...) {

    # -- Size guard -------------------------------------------------------------
    n_vertices <- igraph::vcount(object@network)
    n_edges    <- igraph::ecount(object@network)

    if (!force && (n_vertices > 500L || n_edges > 2000L))
      stop(
        sprintf("The full session network has %d vertices and %d edges, ",
                n_vertices, n_edges),
        "which is too large to plot meaningfully.\n",
        "Consider using:\n",
        "  plot_network(get_analysis(obj, label))              ",
        "# cross-gene subnetwork\n",
        "  plot_network(get_gene_result(obj, label, gene_id))  ",
        "# single gene\n",
        "  plot_network(get_lookup(obj, label))                ",
        "# single lookup\n",
        "To plot the full network anyway: plot_network(obj, force = TRUE)",
        call. = FALSE
      )

    level <- match.arg(level, c("isoform", "gene"))
    node_labels <- match.arg(node_labels)
    .check_plot_args(all_edges, group_by_gene, level)
    .check_plot_pkgs()

    # Pull the graph and threshold out of the session object, then delegate
    g <- .build_plot_graph(
      graph     = object@network,
      level     = level,
      all_edges = all_edges,
      threshold = object@metadata[["pi_threshold"]],
      seed_vertices = character(0)
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
      node_labels   = node_labels,
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
           node_labels   = c("symbol", "id", "none"),
           ...) {

    level <- match.arg(level, c("isoform", "gene"))
    node_labels <- match.arg(node_labels)
    .check_plot_args(all_edges, group_by_gene, level)
    .check_plot_pkgs()

    # Resolve the actual order-0 seed vertices before filtering/contraction so
    # isolated seeds can be protected by .build_plot_graph().
    raw_seed_vertices <- igraph::V(object@subgraph)$name[
      !is.na(igraph::V(object@subgraph)$order) &
        igraph::V(object@subgraph)$order == 0
    ]

    # For a lookup the threshold comes from params, not session metadata
    g <- .build_plot_graph(
      graph     = object@subgraph,
      level     = level,
      all_edges = all_edges,
      threshold = object@params[["cutoff"]],
      seed_vertices = raw_seed_vertices
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
      seed_gene <- igraph::V(object@subgraph)$gene[
        igraph::V(object@subgraph)$name %in% raw_seed_vertices
      ]
      unique(stats::na.omit(seed_gene))
    } else {
      raw_seed_vertices
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
      node_labels   = node_labels,
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
#' @param node_labels   "symbol", "id" or "none".
#' @param ...           Forwarded to \code{ggraph::ggraph()}.
#'
#' @return A ggplot object, invisibly.
#' @noRd
.render_plot <- function(g, level, group_by_gene, layout,
                         node_size, label_size, edge_alpha,
                         title, seed_vertices, node_labels = "symbol", ...) {

  subtitle <- sprintf(
    "%d nodes  \u00b7  %d edges",
    igraph::vcount(g), igraph::ecount(g)
  )

  # Mark seed nodes as a vertex attribute so ggraph can access it via aes()
  igraph::V(g)$is_seed <- igraph::V(g)$name %in% seed_vertices
  g <- .add_node_labels(g, node_labels)

  p <- ggraph::ggraph(g, layout = .plot_layout(g, layout), ...) +
    ggraph::geom_edge_link(alpha = edge_alpha, colour = "grey60") +
    .node_geom(level         = level,
               group_by_gene = group_by_gene,
               node_size     = node_size,
               has_seeds     = any(igraph::V(g)$is_seed)) +
    .label_geom(label_size = label_size, node_labels = node_labels) +
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

  p
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
#' @param seed_vertices Character vector of seed vertex names in the input
#'   graph. Isolated vertices in this set are protected from cleanup.
#'
#' @return A prepared igraph object. Isolated seed vertices are retained;
#'   isolated non-seed vertices are removed.
#' @noRd
.build_plot_graph <- function(graph, level, all_edges, threshold = NULL,
                              seed_vertices = character(0)) {

  g <- graph
  protected_seeds <- intersect(seed_vertices, igraph::V(g)$name)

  # -- Edge filtering ----------------------------------------------------------
  # Only acts when the pi_score attribute exists AND all_edges is FALSE AND a
  # threshold is available.  If the graph is already fully filtered (no
  # pi_score attribute) this block is a no-op.
  if (!all_edges &&
      !is.null(igraph::edge_attr(g, "pi_score")) &&
      !is.null(threshold) && !is.na(threshold)) {
    keep <- igraph::E(g)$pi_score >= threshold
    g    <- igraph::subgraph_from_edges(g, igraph::E(g)[keep],
                                        delete.vertices = FALSE)
  }

  # -- Gene-level contraction --------------------------------------------------
  if (level == "gene") {
    gene_attr <- igraph::V(g)$gene
    if (!is.null(gene_attr)) {
      gene_attr[is.na(gene_attr)] <- "NA_gene"
      protected_seeds <- unique(stats::na.omit(
        gene_attr[igraph::V(g)$name %in% protected_seeds]
      ))
    }
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

  # Keep biologically meaningful isolated seeds, while retaining the original
  # cleanup for isolated non-seed vertices (especially in full-network plots).
  drop <- igraph::degree(g) == 0L &
    !(igraph::V(g)$name %in% protected_seeds)
  igraph::delete_vertices(g, drop)
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

#' Direction palette for neighbourhood membership
#'
#' Matches the analysis figures: green for up, orange for down, and blue for
#' membership in both directional neighbourhoods.
#' @noRd
.membership_palette <- c(
  up_only   = "#009E73",   # green
  shared    = "#0072B2",   # blue
  down_only = "#D55E00"    # orange
)

#' Edge colours matching the membership palette
#' @noRd
.edge_membership_palette <- .membership_palette

#' Seed/neighbour two-colour palette
#' @noRd
.seed_colour   <- "#CC79A7"   # approved reddish purple (Okabe-Ito)
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


#' Layout for ggraph(), avoiding the stress layout on two-node graphs
#'
#' graphlayouts' stress layout warns on graphs with exactly two vertices,
#' although the layout it returns is a straight line; "linear" gives the same
#' picture without the warning.
#' @noRd
.plot_layout <- function(g, layout) {
  if (identical(layout, "stress") && igraph::vcount(g) == 2L) "linear" else layout
}


#' Label geom: repelled text for all nodes, or none
#'
#' Reads the \code{node_label} vertex attribute set by
#' \code{.add_node_labels()}.
#' @noRd
.label_geom <- function(label_size, node_labels = "symbol") {
  if (node_labels == "none") return(NULL)
  ggraph::geom_node_text(
    ggplot2::aes(label = .data[["node_label"]]),
    size  = label_size,
    repel = TRUE
  )
}

#' Add the text of each node's label as a \code{node_label} vertex attribute
#'
#' \code{"symbol"} uses a \code{gene_name} vertex attribute when the network
#' has one, otherwise the bundled gene symbols looked up through the
#' \code{gene} attribute (or the vertex name of gene-level graphs). Nodes
#' without a symbol keep their ID. \code{"id"} and \code{"none"} use the
#' vertex name; \code{"none"} labels are not drawn.
#' @noRd
.add_node_labels <- function(g, node_labels) {
  ids <- igraph::V(g)$name
  labels <- ids
  if (node_labels == "symbol" && length(ids) > 0L) {
    symbols <- igraph::vertex_attr(g, "gene_name")
    if (is.null(symbols)) {
      genes <- igraph::vertex_attr(g, "gene")
      symbols <- .gene_symbols(if (is.null(genes)) ids else genes)
    }
    known <- !is.na(symbols) & nzchar(symbols)
    labels[known] <- symbols[known]
  }
  igraph::V(g)$node_label <- labels
  g
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
#'   \item \code{"overlay"} - merges the up- and down-labelled subgraphs
#'     into a single network. Non-seed nodes and edges use green for the
#'     up-seed neighbourhood, orange for the down-seed neighbourhood, and blue
#'     when present in both. Focal seeds use a reddish-purple fill and directional
#'     seed shapes, separating node role from neighbourhood membership.
#'     Down-only edges are dashed and shared edges thicker, and the legend
#'     sits below the plot.
#' }
#'
#' Node shapes and black outlines provide redundant role encoding alongside
#' colour.
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
#'   \code{ggraph::ggraph()}.
#' @param direction Character scalar: \code{"up"}, \code{"down"},
#'   \code{"side_by_side"} (default), or \code{"overlay"}.  Must be passed
#'   as a named argument.
#'
#' @return A \code{ggplot} / \code{patchwork} object, invisibly.
#'
#' @examples
#' seeds <- read.csv(system.file("extdata", "gse141432_t2d_vs_lean_seeds.csv",
#'                               package = "IsoNetAnalyzeR"))
#' obj <- IsoNetAnalyzeR("example")
#' obj <- run_differential_pulldown(
#'   obj, "t2d", seeds[seeds$gene_name == "CLTCL1", ],
#'   order = 1, confirm = FALSE, verbose = FALSE
#' )
#' gr <- get_gene_result(obj, "t2d", "ENSG00000070371")
#'
#' # Up and down neighbourhoods side by side
#' p <- plot_network(gr)
#'
#' # Both neighbourhoods in one overlay
#' p <- plot_network(gr, direction = "overlay")
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
           node_labels   = c("symbol", "id", "none"),
           ...,
           direction     = c("side_by_side", "up", "down", "overlay")) {

    direction   <- match.arg(direction)
    level       <- match.arg(level, c("isoform", "gene"))
    node_labels <- match.arg(node_labels)
    .check_plot_args(all_edges = FALSE, group_by_gene = group_by_gene, level = level)
    .check_plot_pkgs()

    if (direction == "side_by_side") {
      return(.plot_gene_side_by_side(object, level, group_by_gene, layout,
                                     node_size, label_size, edge_alpha,
                                     title, node_labels, ...))
    }

    if (direction == "overlay") {
      return(.plot_gene_overlay(object, level, group_by_gene, layout,
                                node_size, label_size, edge_alpha,
                                title, node_labels, ...))
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
      threshold = NULL,       # graph is already interaction-filtered
      seed_vertices = iso_seed
    )

    if (is.null(title))
      title <- sprintf("%s  -  %s  (%s level)",
                       .gene_title(object@gene_id), direction, level)

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
      node_labels   = node_labels,
      ...
    )

    p
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
#' @param node_labels   "symbol", "id" or "none".
#' @param ...           Forwarded to ggraph.
#' @return patchwork object, invisibly.
#' @noRd
.plot_gene_side_by_side <- function(object, level, group_by_gene, layout,
                                    node_size, label_size, edge_alpha,
                                    title, node_labels, ...) {

  gene_title <- .gene_title(object@gene_id)

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
          ggplot2::labs(title = paste(gene_title, "-", dir_label))
      )
    }

    g <- .build_plot_graph(
      graph     = g_slot,
      level     = level,
      all_edges = FALSE,
      threshold = NULL,
      seed_vertices = iso_seed
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
      title         = paste(gene_title, "-", dir_label),
      seed_vertices = seed_vertices,
      node_labels   = node_labels,
      ...
    )
  }

  p_up   <- make_panel(object@up_graph,   object@up_isoforms,   "up seeds")
  p_down <- make_panel(object@down_graph, object@down_isoforms, "down seeds")

  combined <- p_up + p_down +
    patchwork::plot_annotation(
      title   = title %||% sprintf("Differential pulldown  -  %s", gene_title),
      caption = sprintf("level: %s", level)
    )

  combined
}


# ==============================================================================
#  Overlay rendering - internal
# ==============================================================================

#' Merge up and down subgraphs into a single overlay graph
#'
#' Performs a graph union of the up- and down-labelled subgraphs and annotates
#' every vertex and edge with a \code{membership} attribute:
#' \code{"up_only"}, \code{"down_only"}, or \code{"shared"}.
#'
#' @param up_graph   igraph or NULL.
#' @param down_graph igraph or NULL.
#' @param up_isoforms   Character vector of up-labelled seed isoform IDs.
#' @param down_isoforms Character vector of down-labelled seed isoform IDs.
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
#' Merges up- and down-labelled subgraphs and renders them as a single
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
#' @param node_labels   "symbol", "id" or "none".
#' @param ...           Forwarded to ggraph.
#' @return A ggplot object, invisibly.
#' @noRd
.plot_gene_overlay <- function(object, level, group_by_gene, layout,
                               node_size, label_size, edge_alpha,
                               title, node_labels, ...) {

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

  # Preserve isolated seed nodes. Isolated non-seeds add no information to the
  # overlay and can still be removed for readability.
  drop <- igraph::degree(g) == 0L &
    igraph::V(g)$seed_direction == "non_seed"
  g <- igraph::delete_vertices(g, drop)

  if (igraph::vcount(g) == 0L)
    stop("Overlay graph has no vertices to plot.", call. = FALSE)

  g <- .add_node_labels(g, node_labels)

  # -- Colour palettes (CVD-safe, Okabe-Ito) ----------------------------------
  membership_colours <- .membership_palette
  edge_colours       <- .edge_membership_palette

  # All shapes support fill + stroke so both branches render consistently
  seed_shapes <- c(
    non_seed   = 21L,   # filled circle
    up_seed    = 24L,   # filled triangle up
    down_seed  = 25L,   # filled triangle down
    both_seeds = 23L    # filled diamond
  )

  # -- Build ggraph plot -------------------------------------------------------
  if (is.null(title))
    title <- .gene_title(object@gene_id)

  subtitle <- sprintf(
    paste("Predicted neighbourhoods, %s level  \u00b7  %d nodes, %d edges",
          " \u00b7  up-only %d | shared %d | down-only %d"),
    level, igraph::vcount(g), igraph::ecount(g),
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

  membership_labels <- c(up_only = "Up-only", shared = "Shared",
                         down_only = "Down-only")

  # Down-only edges are dashed and shared edges thicker, so membership stays
  # readable without colour
  edge_linetypes <- c(up_only = "solid", shared = "solid",
                      down_only = "dashed")
  edge_widths    <- c(up_only = 0.5, shared = 0.9, down_only = 0.5)
  edge_scales <- if (igraph::ecount(g) > 0L) {
    list(
      ggraph::scale_edge_colour_manual(
        name = "Edge", values = edge_colours, labels = membership_labels,
        drop = TRUE
      ),
      ggraph::scale_edge_linetype_manual(
        name = "Edge", values = edge_linetypes, labels = membership_labels,
        drop = TRUE
      ),
      ggraph::scale_edge_width_manual(
        name = "Edge", values = edge_widths, labels = membership_labels,
        drop = TRUE
      )
    )
  }
  edge_aes <- ggplot2::aes(
    edge_colour   = .data[["membership"]],
    edge_linetype = .data[["membership"]],
    edge_width    = .data[["membership"]]
  )
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

    p <- ggraph::ggraph(g, layout = .plot_layout(g, layout), ...) +
      ggraph::geom_edge_link(edge_aes, alpha = edge_alpha) +
      edge_scales +
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
        values = seed_shapes,
        labels = seed_labels,
        drop   = TRUE
      ) +
      ggplot2::scale_size_manual(
        name   = "Seed status",
        values = size_vals,
        labels = seed_labels,
        drop   = TRUE
      ) +
      .label_geom(label_size, node_labels) +
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
        edge_linetype = ggplot2::guide_legend(title = "Edge", order = 3),
        edge_width  = ggplot2::guide_legend(title = "Edge", order = 3),
        shape       = ggplot2::guide_legend(title = "Seed status", order = 4),
        size        = ggplot2::guide_legend(title = "Seed status", order = 4)
      )

  } else {

    # ---- Default overlay: fill = role / neighbourhood -----------------------
    # Direction colours describe non-seed neighbourhood membership. Seed fill
    # is deliberately independent so a shared neighbour cannot be mistaken for
    # the focal seed. Shape continues to retain up/down/both seed direction.
    node_fill <- ifelse(
      as.character(igraph::V(g)$seed_direction) == "non_seed",
      as.character(igraph::V(g)$membership),
      "focal_seed"
    )
    igraph::V(g)$node_fill <- factor(
      node_fill,
      levels = c("up_only", "shared", "down_only", "focal_seed")
    )
    node_fill_colours <- c(membership_colours, focal_seed = .seed_colour)
    node_fill_labels <- c(
      up_only = "Up-only",
      shared = "Shared",
      down_only = "Down-only",
      focal_seed = "Seed"
    )

    p <- ggraph::ggraph(g, layout = .plot_layout(g, layout), ...) +
      ggraph::geom_edge_link(edge_aes, alpha = edge_alpha) +
      edge_scales +
      ggraph::geom_node_point(
        ggplot2::aes(
          fill   = .data[["node_fill"]],
          shape  = .data[["seed_direction"]],
          size   = .data[["seed_direction"]]
        ),
        colour = "black",
        stroke = 0.65
      ) +
      ggplot2::scale_fill_manual(
        name   = "Predicted neighbourhood",
        values = node_fill_colours,
        labels = node_fill_labels,
        drop   = TRUE
      ) +
      ggplot2::scale_shape_manual(
        name   = "Node",
        values = seed_shapes,
        labels = seed_labels,
        drop   = TRUE
      ) +
      ggplot2::scale_size_manual(
        name   = "Node",
        values = size_vals,
        labels = seed_labels,
        drop   = TRUE
      ) +
      .label_geom(label_size, node_labels) +
      ggplot2::labs(title = title, subtitle = subtitle) +
      ggraph::theme_graph(base_family = "sans") +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(size = 13, face = "bold",
                                                hjust = 0.5),
        plot.subtitle   = ggplot2::element_text(size = 9, colour = "grey40",
                                                hjust = 0.5),
        legend.position = "bottom",
        legend.box      = "vertical",
        legend.title    = ggplot2::element_text(face = "bold", size = 10),
        legend.margin   = ggplot2::margin(0, 0, 0, 0)
      ) +
      ggplot2::guides(
        fill          = ggplot2::guide_legend(
          order = 1,
          override.aes = list(shape = 21L, size = 3.5, colour = "black")
        ),
        edge_colour   = ggplot2::guide_legend(order = 2),
        edge_linetype = ggplot2::guide_legend(order = 2),
        edge_width    = ggplot2::guide_legend(order = 2),
        shape         = ggplot2::guide_legend(order = 3),
        size          = ggplot2::guide_legend(order = 3)
      )
  }

  p
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
#'   \code{ggraph::ggraph()}.
#' @param direction Character scalar: \code{"up"} (merge all up-labelled
#'   subgraphs), \code{"down"} (merge all down-labelled subgraphs), or
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
#' seeds <- read.csv(system.file("extdata", "gse141432_t2d_vs_lean_seeds.csv",
#'                               package = "IsoNetAnalyzeR"))
#' genes <- c("CLTCL1", "WAS", "VAV1", "HLA-DMA", "ICOS", "PRKACB")
#' obj <- IsoNetAnalyzeR("example")
#' obj <- run_differential_pulldown(
#'   obj, "t2d", seeds[seeds$gene_name %in% genes, ],
#'   order = 1, confirm = FALSE, verbose = FALSE
#' )
#' result <- get_analysis(obj, "t2d")
#'
#' # Cross-gene view of the up-labelled seed neighbourhoods
#' p <- plot_network(result, direction = "up")
#'
#' # Keep components that share at least one node between genes
#' p <- plot_network(result, min_shared_nodes = 1)
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
           node_labels = c("symbol", "id", "none"),
           ...,
           direction          = c("both", "up", "down"),
           genes              = NULL,
           min_nodes          = NULL,
           min_shared_nodes   = NULL,
           min_component_nodes = NULL,
           min_edges          = NULL,
           colour_by          = c("gene", "hub_degree")) {

    direction   <- match.arg(direction)
    level       <- match.arg(level, c("isoform", "gene"))
    colour_by   <- match.arg(colour_by)
    node_labels <- match.arg(node_labels)
    .check_plot_pkgs()

    # -- Validate filter arguments ---------------------------------------------
    .check_filter_arg <- function(x, name) {
      if (!is.null(x)) .validate_positive_integer(x, name)
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
                        up   = "up seeds",
                        down = "down seeds",
                        both = "up + down")

    if (is.null(title))
      title <- sprintf(
        "Cross-gene network: %s  (%d gene(s), %s, %s level)",
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
      node_labels = node_labels,
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
#' @return The filtered igraph. Explicit filters may remove singleton
#'   components, but unfiltered isolated seed components are retained.
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
#' @param node_labels "symbol", "id" or "none".
#' @param ...        Forwarded to ggraph.
#' @return A ggplot object, invisibly.
#' @noRd
.render_cross_gene_plot <- function(g, seed_names, level, layout,
                                    node_size, label_size, edge_alpha,
                                    title, colour_by = "gene",
                                    node_labels = "symbol", ...) {

  g <- .add_node_labels(g, node_labels)

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

    p <- ggraph::ggraph(g, layout = .plot_layout(g, layout), ...) +
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
      .label_geom(label_size, node_labels) +
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

    return(p)
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

  p <- ggraph::ggraph(g, layout = .plot_layout(g, layout), ...) +
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
    .label_geom(label_size, node_labels) +
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

  p
}
