## plot_component.R
## User-facing function to plot a single connected component from a
## ComponentAnalysis result.
##
## Rebuilds the merged cross-gene graph (via .build_cross_gene_graph, the same
## helper used by run_component_analysis and plot_network), sorts components
## by descending node count / edge count to match the IDs in the
## ComponentAnalysis tables, extracts the requested component, and delegates
## rendering to .render_cross_gene_plot (gene / hub_degree modes) or the new
## direction-aware renderers.
##
## colour_by = "direction" behaviour:
##
##   ISOFORM LEVEL
##     Each node is one isoform.  Colour encodes its seed direction:
##       up       → orange   (#E69F00)
##       down     → sky blue (#56B4E9)
##       both     → purple   (#CC79A7)
##       neighbour → grey    (#CCCCCC)  [or territory colour if colour_neighbours]
##
##     colour_neighbours = TRUE (territory map):
##       Each neighbour is coloured by its nearest seed group via unweighted
##       shortest-path distance (igraph::distances()).
##         dist_up < dist_down  → up territory   (orange)
##         dist_down < dist_up  → down territory (sky blue)
##         equidistant          → contested      (purple)
##         unreachable          → grey fallback
##
##   GENE LEVEL
##     Each node is a gene collapsed from multiple isoforms.  A binary
##     direction label is meaningless (every gene will have both up and down
##     isoforms).  Instead a continuous dominance score is computed per gene:
##
##       score = (n_up_isoforms - n_down_isoforms) /
##               (n_up_isoforms + n_down_isoforms)   ∈ [-1, +1]
##
##       +1  all contributing isoforms are up-regulated  → orange
##        0  equal up and down (or a pure neighbour)     → white / neutral
##       -1  all contributing isoforms are down-regulated → sky blue
##
##     Rendered on a diverging orange–white–blue colour scale.
##     Neighbour genes (no seed isoforms) receive score = 0 (neutral white).
##     colour_neighbours is silently ignored at gene level (dominance already
##     encodes the full picture).
##
## Exported: plot_component

#' @include isonetanalyzer_class.R
NULL
#' @include diff_pulldown_analysis_class.R
NULL
#' @include component_analysis.R
NULL


# ==============================================================================
#  Exported function
# ==============================================================================

#' Plot a single connected component from a ComponentAnalysis result
#'
#' Extracts one connected component by its \code{component_id} (as assigned in
#' \code{@@component_table}) and renders it using the same cross-gene plot
#' style as \code{\link{plot_network}} on a
#' \linkS4class{DifferentialPulldownAnalysis} object.
#'
#' Component IDs are assigned in descending order of node count (edge count as
#' tiebreaker), so \code{component_id = 1} is always the largest component.
#' These IDs match the \code{component_id} column in all four
#' \linkS4class{ComponentAnalysis} tables.
#'
#' @section Colour modes (\code{colour_by}):
#' \describe{
#'   \item{\code{"gene"}}{Each node is coloured by the pulldown gene it belongs
#'     to.  Shared (multi-gene) nodes are dark grey.}
#'   \item{\code{"hub_degree"}}{Shared nodes are coloured by degree on a
#'     sequential (Mako) scale; all other nodes are light grey.}
#'   \item{\code{"direction"}}{Encodes the regulation direction of each node.
#'     Behaviour differs by analysis level — see below.}
#' }
#'
#' @section Direction mode — isoform level:
#' Seed nodes are coloured by their caller-supplied direction: up (orange),
#' down (sky blue), or present in both directions (purple).
#' Neighbour behaviour is controlled by \code{colour_neighbours}:
#' \describe{
#'   \item{\code{FALSE} (default)}{Neighbours are light grey.}
#'   \item{\code{TRUE} — territory map}{Each neighbour is coloured by the
#'     direction of its nearest seed group, computed as unweighted shortest-path
#'     distance within the component.  Neighbours equidistant from both seed
#'     groups are shown in purple (equidistant territory).}
#' }
#'
#' @section Direction mode — gene level:
#' At gene level each node represents a collapsed gene rather than a single
#' isoform, so a binary up/down label is not meaningful (most genes contribute
#' both up- and down-labelled isoforms). Instead a continuous
#' \strong{dominance score} is computed per gene:
#' \deqn{score = \frac{n_{up} - n_{down}}{n_{up} + n_{down}} \in [-1, +1]}
#' where \eqn{n_{up}} and \eqn{n_{down}} are the numbers of seed isoforms from
#' the input table that are up- and down-labelled for that gene. The score is
#' mapped to a diverging orange–white–sky-blue colour scale:
#' \itemize{
#'   \item \strong{+1} (pure up) → orange
#'   \item \strong{0} (balanced or neighbour) → white / neutral
#'   \item \strong{-1} (pure down) → sky blue
#' }
#' Neighbour genes (no seed isoforms in this component) receive a score of 0
#' and appear white/neutral.  \code{colour_neighbours} is ignored at gene level.
#'
#' @param obj                An \linkS4class{IsoNetAnalyzeR} session object.
#' @param analysis           Character scalar.  Label of the parent
#'   \linkS4class{DifferentialPulldownAnalysis}; use \code{list_analyses(obj)}
#'   to see available labels.
#' @param component_analysis Character scalar.  Label of the
#'   \linkS4class{ComponentAnalysis}; use \code{list_analyses(obj)} then
#'   inspect the result to see available component analysis labels.
#' @param component_id       Integer scalar.  The component to plot, matching
#'   the \code{component_id} column in \code{get_component_overview()}.
#' @param colour_by          Character scalar: \code{"gene"} (default),
#'   \code{"hub_degree"}, or \code{"direction"}.
#' @param colour_neighbours  Logical.  Only used when
#'   \code{colour_by = "direction"} and \code{level = "isoform"}.  When
#'   \code{TRUE} neighbour nodes are coloured by nearest-seed territory instead
#'   of neutral grey.  Default \code{FALSE}.
#' @param layout             ggraph layout string.  Default \code{"stress"}.
#' @param node_size          Numeric.  Base node size.  Default \code{3}.
#' @param label_size         Numeric.  Label text size.  Default \code{2.5}.
#' @param edge_alpha         Numeric in \[0, 1\].  Edge opacity.  Default
#'   \code{0.4}.
#' @param title              Character scalar or \code{NULL} (auto-generated).
#' @param node_labels        Character scalar.  \code{"symbol"} (default)
#'   labels nodes with gene symbols, \code{"id"} with transcript (or, at gene
#'   level, gene) identifiers, and \code{"none"} draws no labels.  See
#'   \code{\link{plot_network}}.
#' @param ...                Additional arguments forwarded to
#'   \code{ggraph::ggraph()}.
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
#' obj <- run_component_analysis(obj, "t2d", label = "comp")
#'
#' # Colour seeds by direction
#' p <- plot_component(obj, "t2d", "comp", component_id = 1,
#'                     colour_by = "direction")
#'
#' # Territory map: neighbours coloured by the seed directions reaching them
#' p <- plot_component(obj, "t2d", "comp", component_id = 1,
#'                     colour_by = "direction", colour_neighbours = TRUE)
#'
#' # Without node labels
#' p <- plot_component(obj, "t2d", "comp", component_id = 1,
#'                     node_labels = "none")
#'
#' @export
plot_component <- function(obj,
                           analysis,
                           component_analysis,
                           component_id,
                           colour_by         = c("gene", "hub_degree", "direction"),
                           colour_neighbours = FALSE,
                           layout            = "stress",
                           node_size         = 3,
                           label_size        = 2.5,
                           edge_alpha        = 0.4,
                           title             = NULL,
                           node_labels       = c("symbol", "id", "none"),
                           ...) {

  # -- Input validation --------------------------------------------------------
  if (!is(obj, "IsoNetAnalyzeR"))
    stop("`obj` must be an IsoNetAnalyzeR object.", call. = FALSE)

  if (!is.character(analysis) || length(analysis) != 1L ||
      nchar(trimws(analysis)) == 0L)
    stop("`analysis` must be a single non-empty character string.", call. = FALSE)

  if (!analysis %in% names(obj@analyses))
    stop(sprintf(
      "No analysis labelled '%s'. Available: %s",
      analysis, paste(names(obj@analyses), collapse = ", ")
    ), call. = FALSE)

  diff_result <- obj@analyses[[analysis]]
  if (!is(diff_result, "DifferentialPulldownAnalysis"))
    stop(sprintf(
      "Analysis '%s' is not a DifferentialPulldownAnalysis object.", analysis
    ), call. = FALSE)

  if (!is.character(component_analysis) || length(component_analysis) != 1L ||
      nchar(trimws(component_analysis)) == 0L)
    stop("`component_analysis` must be a single non-empty character string.",
         call. = FALSE)

  if (!component_analysis %in% names(diff_result@component_analyses))
    stop(sprintf(
      "No component analysis labelled '%s' in analysis '%s'. Available: %s",
      component_analysis, analysis,
      paste(names(diff_result@component_analyses), collapse = ", ")
    ), call. = FALSE)

  ca <- diff_result@component_analyses[[component_analysis]]

  component_id <- .validate_positive_integer(component_id, "component_id")

  n_comps <- nrow(ca@component_table)
  if (component_id > n_comps)
    stop(sprintf(
      "`component_id` = %d is out of range; this analysis has %d component(s).",
      component_id, n_comps
    ), call. = FALSE)

  colour_by   <- match.arg(colour_by)
  node_labels <- match.arg(node_labels)

  if (!is.logical(colour_neighbours) || length(colour_neighbours) != 1L ||
      is.na(colour_neighbours))
    stop("`colour_neighbours` must be a single non-NA logical value.",
         call. = FALSE)

  if (colour_neighbours && colour_by != "direction")
    warning("`colour_neighbours` is ignored when `colour_by` is not 'direction'.",
            call. = FALSE)

  if (colour_neighbours && ca@level == "gene")
    message("`colour_neighbours` is ignored at gene level; dominance score is used instead.")

  .check_plot_pkgs()

  # -- Resolve parameters from ComponentAnalysis -------------------------------
  direction  <- ca@direction
  level      <- ca@level
  genes_used <- if (length(ca@genes) == 0L) {
    names(diff_result@gene_results)
  } else {
    ca@genes
  }

  gene_results_subset <- diff_result@gene_results[genes_used]

  # -- Rebuild and annotate merged graph ---------------------------------------
  # Mirrors run_component_analysis exactly so component ordering is identical.
  built <- .build_cross_gene_graph(gene_results_subset, direction, level)
  g     <- built$graph

  if (igraph::vcount(g) == 0L)
    stop("The merged cross-gene graph is empty; nothing to plot.", call. = FALSE)

  g <- .annotate_seed_direction(g, gene_results_subset, direction, ca@level)

  # -- Extract exactly the component stored by run_component_analysis() --------
  # Component analyses may have explicit component filters. Re-decomposing and
  # sorting the unfiltered graph here can therefore assign a different graph to
  # the stored component_id. The node table is the authoritative membership
  # record and also preserves valid one-node/zero-edge components.
  component_nodes <- ca@node_table$node[
    ca@node_table$component_id == component_id
  ]
  component_nodes <- intersect(component_nodes, igraph::V(g)$name)
  if (length(component_nodes) == 0L)
    stop("The stored component has no vertices in the rebuilt graph.",
         call. = FALSE)
  cg <- igraph::induced_subgraph(g, vids = component_nodes)
  seed_names <- igraph::V(cg)$name[!is.na(igraph::V(cg)$seed_direction)]

  # -- Auto-generate title -----------------------------------------------------
  if (is.null(title)) {
    ct        <- .df(ca@component_table)
    ct_row    <- ct[ct$component_id == component_id, ]
    mixed_tag <- if (isTRUE(ct_row$mixed)) "  \u00b7  mixed" else ""
    title <- sprintf(
      "Component %d: %s / %s  (%s level%s)",
      component_id, analysis, ca@label, level, mixed_tag
    )
  }

  # -- Dispatch to renderer ----------------------------------------------------
  if (colour_by == "direction") {
    if (level == "gene") {
      # Gene level: dominance score on a diverging scale
      .render_component_dominance_plot(
        g           = cg,
        input_table = .df(diff_result@input_table),
        layout      = layout,
        node_size   = node_size,
        label_size  = label_size,
        edge_alpha  = edge_alpha,
        title       = title,
        node_labels = node_labels,
        ...
      )
    } else {
      # Isoform level: discrete direction colours ± territory map
      .render_component_direction_plot(
        g                 = cg,
        colour_neighbours = colour_neighbours,
        layout            = layout,
        node_size         = node_size,
        label_size        = label_size,
        edge_alpha        = edge_alpha,
        title             = title,
        node_labels       = node_labels,
        ...
      )
    }
  } else {
    .render_cross_gene_plot(
      g          = cg,
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
}


# ==============================================================================
#  Internal: territory assignment (isoform level)
# ==============================================================================

#' Assign each neighbour node to the territory of its nearest seed group
#'
#' Uses \code{igraph::distances()} to compute the shortest-path distance from
#' every node to the up-seed group and to the down-seed group.  Returns a named
#' character vector with one entry per node:
#'   "up"        closer to up seeds
#'   "down"      closer to down seeds
#'   "both"      equidistant (contested zone)
#'   "neighbour" unreachable from any seed (grey fallback)
#'
#' Seed nodes are pre-assigned their own direction and are never overridden.
#'
#' @param g igraph for one component, with \code{seed_direction} attribute.
#' @return Named character vector, one entry per vertex.
#' @noRd
.assign_territory <- function(g) {

  vnames   <- igraph::V(g)$name
  seed_dir <- igraph::V(g)$seed_direction

  up_seeds   <- vnames[!is.na(seed_dir) & seed_dir %in% c("up",   "both")]
  down_seeds <- vnames[!is.na(seed_dir) & seed_dir %in% c("down", "both")]

  # Initialise: seeds keep their own direction, neighbours start NA
  territory <- seed_dir
  names(territory) <- vnames

  neighbour_idx <- which(is.na(seed_dir))
  if (length(neighbour_idx) == 0L)
    return(ifelse(is.na(territory), "neighbour", territory))

  # Minimum distance from each seed group to every node
  dist_up <- if (length(up_seeds) > 0L) {
    apply(
      igraph::distances(g, v = up_seeds, to = igraph::V(g),
                        mode = "all", algorithm = "unweighted"),
      MARGIN = 2, FUN = min
    )
  } else {
    stats::setNames(rep(Inf, igraph::vcount(g)), vnames)
  }

  dist_down <- if (length(down_seeds) > 0L) {
    apply(
      igraph::distances(g, v = down_seeds, to = igraph::V(g),
                        mode = "all", algorithm = "unweighted"),
      MARGIN = 2, FUN = min
    )
  } else {
    stats::setNames(rep(Inf, igraph::vcount(g)), vnames)
  }

  names(dist_up)   <- vnames
  names(dist_down) <- vnames

  # Classify neighbours only (seeds keep their pre-assigned direction)
  for (idx in neighbour_idx) {
    v   <- vnames[idx]
    d_u <- dist_up[v]
    d_d <- dist_down[v]

    territory[v] <- if (is.infinite(d_u) && is.infinite(d_d)) {
      "neighbour"      # unreachable fallback
    } else if (is.infinite(d_d)) {
      "up"             # only reachable from up seeds
    } else if (is.infinite(d_u)) {
      "down"           # only reachable from down seeds
    } else if (d_u < d_d) {
      "up"
    } else if (d_d < d_u) {
      "down"
    } else {
      "both"           # equidistant → contested zone
    }
  }

  territory
}


# ==============================================================================
#  Internal: dominance score computation (gene level)
# ==============================================================================

#' Compute per-gene dominance scores for the direction mode at gene level
#'
#' For each gene node in the component, counts the number of up- and
#' down-labelled seed isoforms from the original input table that map to that
#' gene, then computes:
#'
#'   score = (n_up - n_down) / (n_up + n_down)   in \[-1, +1\]
#'
#' Gene nodes with no seed isoforms (pure neighbours) receive score = 0.
#'
#' At gene level the vertex \code{name} attribute holds the gene ID (after
#' contraction by \code{.contract_to_gene}), so we can join directly against
#' the \code{gene_id} column of the input table.
#'
#' @param g           igraph for one component at gene level.
#' @param input_table The \code{@@input_table} data frame from the parent
#'   \linkS4class{DifferentialPulldownAnalysis}.  Must have columns
#'   \code{gene_id} and \code{direction} ("up"/"down").
#'
#' @return Named numeric vector of dominance scores, one per vertex.
#' @noRd
.compute_dominance_scores <- function(g, input_table) {

  gene_ids <- igraph::V(g)$name

  # The analysis input table contains resolved protein representatives.
  # Dominance must not double-count protein-identical submitted transcripts.

  scores <- vapply(gene_ids, function(gid) {
    rows  <- input_table[input_table$gene_id == gid, , drop = FALSE]
    n_up   <- sum(rows$direction == "up")
    n_down <- sum(rows$direction == "down")
    total  <- n_up + n_down
    if (total == 0L) 0  else (n_up - n_down) / total
  }, numeric(1))

  names(scores) <- gene_ids
  scores
}


# ==============================================================================
#  Internal: isoform-level direction renderer
# ==============================================================================

#' Render a component plot coloured by seed direction (isoform level)
#'
#' colour_neighbours = FALSE: seeds coloured, neighbours neutral grey.
#' colour_neighbours = TRUE:  territory map via nearest-seed distance.
#'
#' Colour palette (CVD-safe Okabe-Ito):
#'   up / up territory   → orange   (#E69F00)
#'   down / down terr.   → sky blue (#56B4E9)
#'   both / contested    → purple   (#CC79A7)
#'   unreachable nbr.    → grey     (#CCCCCC)
#'
#' @noRd
.render_component_direction_plot <- function(g, colour_neighbours, layout,
                                             node_size, label_size,
                                             edge_alpha, title,
                                             node_labels = "symbol", ...) {

  seed_dir <- igraph::V(g)$seed_direction
  is_seed  <- !is.na(seed_dir)

  colour_role <- if (colour_neighbours) {
    .assign_territory(g)
  } else {
    ifelse(is_seed, seed_dir, "neighbour")
  }

  # Shape encodes actual seed status regardless of territory colouring
  shape_role <- ifelse(is_seed, "seed", "neighbour")

  all_levels <- c("up", "down", "both", "neighbour")
  igraph::V(g)$colour_role <- factor(colour_role, levels = all_levels)
  igraph::V(g)$shape_role  <- factor(shape_role,
                                     levels = c("seed", "neighbour"))
  g <- .add_node_labels(g, node_labels)

  n_multi  <- sum(igraph::V(g)$source_gene == "multi-gene", na.rm = TRUE)
  subtitle <- sprintf(
    "%d nodes  \u00b7  %d edges  \u00b7  %d multi-gene node(s)",
    igraph::vcount(g), igraph::ecount(g), n_multi
  )

  dir_colours <- c(up = "#E69F00", down = "#56B4E9",
                   both = "#CC79A7", neighbour = "#CCCCCC")

  dir_labels <- if (colour_neighbours) {
    c(up = "Up territory", down = "Down territory",
      both = "Equidistant", neighbour = "Unreachable")
  } else {
    c(up = "Up seed", down = "Down seed",
      both = "Up & down seed", neighbour = "Neighbour")
  }
  colour_title <- if (colour_neighbours) "Nearest seed" else "Direction"

  shape_vals  <- c(seed = 18L, neighbour = 16L)
  shape_names <- c(seed = "Seed", neighbour = "Neighbour")
  size_vals   <- c(seed = node_size * 1.6, neighbour = node_size)

  p <- ggraph::ggraph(g, layout = .plot_layout(g, layout), ...) +
    ggraph::geom_edge_link(alpha = edge_alpha, colour = "grey60") +
    ggraph::geom_node_point(
      ggplot2::aes(colour = .data[["colour_role"]],
                   shape  = .data[["shape_role"]],
                   size   = .data[["shape_role"]])
    ) +
    ggplot2::scale_colour_manual(name = colour_title, values = dir_colours,
                                 labels = dir_labels, drop = TRUE) +
    ggplot2::scale_shape_manual(name  = "Node", values = shape_vals,
                                labels = shape_names, drop = TRUE) +
    ggplot2::scale_size_manual(name   = "Node", values = size_vals,
                               labels = shape_names, drop = TRUE) +
    .label_geom(label_size, node_labels) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle   = ggplot2::element_text(size = 9,  colour = "grey40"),
      legend.position = "right"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1,
                                     override.aes = list(shape = 16L,
                                                         size = 3)),
      shape  = ggplot2::guide_legend(order = 2),
      size   = "none"
    )

  p
}


# ==============================================================================
#  Internal: gene-level dominance renderer
# ==============================================================================

#' Render a component plot coloured by dominance score (gene level)
#'
#' Each gene node is coloured by its dominance score:
#'   score = (n_up_isoforms - n_down_isoforms) / (n_up_isoforms + n_down_isoforms)
#'
#' Mapped to a diverging orange–white–sky-blue scale:
#'   +1 (pure up)      → orange   (#E69F00)
#'    0 (balanced/nbr) → white
#'   -1 (pure down)    → sky blue (#56B4E9)
#'
#' Seed genes (score != 0) are drawn as diamonds; pure neighbour genes
#' (score == 0 and no seed isoforms) as circles.
#'
#' @param g           igraph for one component at gene level.
#' @param input_table data.frame from DifferentialPulldownAnalysis@@input_table.
#' @param layout      ggraph layout string.
#' @param node_size   Numeric.
#' @param label_size  Numeric.
#' @param edge_alpha  Numeric.
#' @param title       Character.
#' @param ...         Forwarded to ggraph::ggraph().
#'
#' @return A ggplot object, invisibly.
#' @noRd
.render_component_dominance_plot <- function(g, input_table, layout,
                                             node_size, label_size,
                                             edge_alpha, title,
                                             node_labels = "symbol", ...) {

  scores  <- .compute_dominance_scores(g, input_table)
  is_seed <- !is.na(igraph::V(g)$seed_direction)

  igraph::V(g)$dominance <- scores
  g <- .add_node_labels(g, node_labels)

  # Seeds (any seed isoform) → diamond; pure neighbours → circle
  igraph::V(g)$node_shape <- ifelse(is_seed, 18L, 16L)
  igraph::V(g)$node_size  <- ifelse(is_seed, node_size * 1.6, node_size)

  n_multi  <- sum(igraph::V(g)$source_gene == "multi-gene", na.rm = TRUE)
  subtitle <- sprintf(
    "%d nodes  \u00b7  %d edges  \u00b7  %d multi-gene node(s)",
    igraph::vcount(g), igraph::ecount(g), n_multi
  )

  # Diverging scale: orange (+1) → white (0) → sky blue (-1)
  # Colours match the up/down palette used throughout the package (Okabe-Ito).
  scale_col <- ggplot2::scale_colour_gradient2(
    low      = "#56B4E9",   # sky blue (down-dominated)
    mid      = "grey85",     # balanced / neighbour
    high     = "#E69F00",   # orange   (up-dominated)
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Dominance\n(up to down)",
    guide    = ggplot2::guide_colorbar(order = 1)
  )

  p <- ggraph::ggraph(g, layout = .plot_layout(g, layout), ...) +
    ggraph::geom_edge_link(alpha = edge_alpha, colour = "grey60") +
    ggraph::geom_node_point(
      ggplot2::aes(colour = .data[["dominance"]]),
      shape = igraph::V(g)$node_shape,
      size  = igraph::V(g)$node_size
    ) +
    scale_col +
    .label_geom(label_size, node_labels) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      caption  = "Diamond = seed gene  \u00b7  Circle = neighbour gene"
    ) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9,  colour = "grey40"),
      plot.caption  = ggplot2::element_text(size = 8,  colour = "grey50"),
      legend.position = "right"
    )

  p
}
