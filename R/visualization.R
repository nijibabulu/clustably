#' @importFrom utils globalVariables
#' @importFrom ggplot2 ggproto GeomViolin
#' @importFrom rlang "%||%"
#'
NULL

#' Plot the cross-classification table of cell identities
#'
#' @param x a seurat object with cells in common with y
#' @param y a seurat object with cells in common with x
#' @param names names to assign to the two cell identities
#' @param scale if true, scale the data
#' @param palette a RColorBrewer palette to use for the coloring
#'
#' @importFrom Seurat Idents
#' @importFrom dplyr mutate
#' @importFrom forcats fct_relevel
#' @importFrom tibble as_tibble
#' @importFrom purrr set_names
#' @importFrom ggplot2 geom_tile aes_string scale_fill_distiller
#'
#' @return an ggplot heatmap of
#'
#' @export
PlotCrossClassification <- function(x,y, names=c("x","y"), scale=F, palette="OrRd") {
  identsTable(Idents(x),Idents(y),names=c("x","y")) %>%
    scale(scale=scale, center=F) %>%
    as_tibble %>%
    mutate(x=fct_relevel(x, labelSort), y=fct_relevel(y, labelSort)) %>%
    rename(!!!set_names( c("x","y"), names)) %>%
    ggplot(aes_string(names[1], names[2], fill="n")) +
    geom_tile() +
    scale_fill_distiller(palette = "OrRd", direction=1)
}

#' Plot the frequency of cross labeling
#'
#' @param obj a seurat object for which has been jackknifed or bootstrapped
#' @param consensusName the slot under which the consensus labeling has been stored
#' @param frequencyName the slot in `misc` under which the labeling frequency has been stored
#' @param ... parameters to pass to Seurat::FeaturePlot
#'
#' @importFrom Seurat FeaturePlot
#' @importFrom stringr str_glue
#' @importFrom purrr flatten_chr flatten_dbl map_chr map_dfc map2_dfc modify_at reduce
#'
#' @return plots the frequency of crosslabeled cells
#'
#' @export
DimCrossClassificationPlot <- function(obj,
                                       consensusName="jackknifeConsensus",
                                       frequencyName="jackknifeFrequency",
                                       ...) {
  .toTitle <- function(label) str_glue("Crosslabeled {label}")
  .addCrossClassification <- function(.obj, label) {
    freqsLabel <- flatten_dbl(crossClassFreqs[which(labels == label),])
    .obj[[.toTitle(label)]] <- freqsLabel
    .obj
  }
  freqs <- map_dfc(obj@misc[[frequencyName]], ~.x/sum(.x))
  labels <- rownames(obj@misc[[frequencyName]])
  cons <- flatten_chr(obj[[consensusName]])
  crossClassFreqs <- map2_dfc(freqs, cons, ~modify_at(.x, which(labels == .y), ~0))
  obj <- reduce(labels, .addCrossClassification, .init=obj)
  Seurat::FeaturePlot(obj, map_chr(labels, .toTitle))
}



#' @importFrom seriation seriate get_order
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid grid.grabExpr
#' @importFrom stringr str_glue
#' @importFrom ggplotify as.ggplot
DoClassificationHeatmap <- function(mat, scale, center, normalize, ordering_method, column_title, row_title, legend_title, show_row_names, ...) {
  if(scale && normalize) {
    stop("scale and normalize can not both be TRUE.")
  }
  if(scale) {
    mat <- scale(mat, center=center)
  } else if(normalize) {
    mat <- mat/rowSums(mat)
  }
  if(ordering_method == "BEA_TSP") {
    if(length(mat) > 1e2) {
      message(str_glue("Attempting to seriate a matrix with {length(mat)} cells. This could take a while. Consider using 'GW' or 'default' for the ordering method"))
    }
    ord = seriate(max(mat)-mat, method="BEA_TSP", verbose=F)
    row_order = get_order(ord, 1)
    column_order = get_order(ord, 2)
  } else if(ordering_method == "GW") {
    row_order = seriate(dist(mat), method="GW") %>% get_order()
    column_order = seriate(dist(t(mat)), method="GW") %>% get_order()
  } else {
    row_order = column_order = NULL
  }
  breaks <- seq(round(min(mat), digits=1), round(max(mat), digits=1), length.out=3)
  hm = Heatmap(mat,  col=brewer.pal(9, "OrRd"),
          show_column_dend=F, show_row_dend = F,
          row_order = row_order, column_order = column_order,
          column_names_side = "bottom", row_names_side="left",
          column_title = column_title, row_title = row_title,
          show_row_names=show_row_names,
          heatmap_legend_param = list(title=legend_title, at=breaks,
                                      direction="horizontal", title_position="lefttop"),
          ...)
  as.ggplot(grid.grabExpr(draw(hm, heatmap_legend_side="bottom")))
}

#' Plot cross-classification frequency in a heatmap
#'
#' @param obj a seurat object which has been jackknifed or bootstrapped
#' @param consensusName the slot under which the consensus labeling has been stored
#' @param frequencyName the slot in `misc` under which the labeling frequency has been stored
#' @param center center the data
#' @param scale scale the frequencies
#' @param plotTitle title of the plot
#' @param rowTitle title along the y axis of the plot
#' @param orderingMethod the method to use for ordering the rows and columns of the heatmap.
#'                       Note that the default, "BEA_TSP" can take a long time. "GW" is a hierarchical clustering method.
#'                       "default" will use the default method for `Heatmap`
#' @param ... parameters to pass to Heatmap
#'
#' @importFrom tidyr gather spread
#' @importFrom dplyr rename group_by summarize
#' @importFrom tibble as_tibble rownames_to_column column_to_rownames
#' @importFrom stringr str_replace str_to_title
#'
#' @export
PlotCrossClassificationFrequency <- function(obj,  consensusName="jackknifeConsensus", frequencyName="jackknifeFrequency",
                                             center=FALSE, scale=FALSE, normalize=TRUE, orderingMethod=c("BEA_TSP", "GW", "default"),
                                             plotTitle="Cross Classification", rowTitle="Consensus", ...) {
  orderingMethod <- match.arg(orderingMethod)
  consVar <- list(cons=consensusName)
  cons <- rownames_to_column(obj[[consensusName]], "cell") %>% rename(!!!consVar)
  freqs <- obj@misc[[frequencyName]] %>%  rownames_to_column("cluster") %>% as_tibble()
  consFreqs <- freqs %>% gather("cell", "count", -cluster) %>% full_join(cons, by="cell")
  consFreqsMat <- consFreqs %>% group_by(cons, cluster) %>% summarize(n=sum(count)) %>% spread(cluster,n) %>%
    column_to_rownames("cons") %>% as.matrix()

  freqTitle <- str_replace(frequencyName, "([[:upper:]])", " \\1") %>% str_to_title()
  DoClassificationHeatmap(consFreqsMat, scale, center, normalize, ordering_method=orderingMethod,
                          column_title=plotTitle, row_title=rowTitle,
                          show_row_names=TRUE, legend_title = freqTitle, ...)

}



#' Plot background classification frequency of individual cells in a heatmap
#'
#' @param obj a seurat object which has been jackknifed or bootstrapped
#' @param consensusName the slot under which the consensus labeling has been stored
#' @param frequencyName the slot in `misc` under which the labeling frequency has been stored
#' @param center center the data
#' @param scale scale the frequencies
#' @param plotTitle title of the plot
#' @param rowTitle title along the y axis of the plot
#' @param orderingMethod the method to use for ordering the rows and columns of the heatmap.
#'                       Note that the default, "BEA_TSP" can take a long time. "GW" is a hierarchical clustering method.
#'                       "default" will use the default method for `Heatmap`
#' @param ... parameters to pass to `Heatmap``
#'
#' @importFrom tidyr gather spread
#' @importFrom dplyr rename group_by summarize
#' @importFrom tibble as_tibble rownames_to_column column_to_rownames
#'
#' @export
PlotCellClassifications <- function(obj,  consensusName="jackknifeConsensus", frequencyName="jackknifeFrequency",
                                             center=FALSE, scale=FALSE, normalize=TRUE, orderingMethod=c("BEA_TSP", "GW", "default"),
                                             plotTitle="Cell Classifications", rowTitle="Cells", ...) {

  orderingMethod <- match.arg(orderingMethod)
  freqs <- obj@misc[[frequencyName]] %>% rownames_to_column("cluster") %>% as_tibble()
  freqsT <- freqs %>%
    gather( key = 'cell', value = 'freq', 2:ncol(.) ) %>%
    spread( key = 1, value = freq ) %>% column_to_rownames("cell") %>% as.matrix()
  freqTitle <- str_replace(frequencyName, "([[:upper:]])", " \\1") %>% str_to_title()
  DoClassificationHeatmap(freqsT, scale, center, normalize, ordering_method = orderingMethod,
                          column_title=plotTitle, row_title=rowTitle,
                          show_row_names=FALSE, legend_title = freqTitle, ...)
}

#' Plot cluster-wise confidence distributions
#'
#' @param obj a seurat object which has been jackknifed or bootstrapped
#' @param consensusName the slot under which the consensus labeling has been stored
#' @param confidenceName store the confidence statistics under this name
#' @param boxplot use a boxplot. If combined with a violin plot, the box will be reduced to width 0.1
#' @param violin use a violin plot
#' @param jitter include jitter points in the plot (use ony with violin plot)
#' @param plotTitle a title for the plot
#'
#' @importFrom ggplot2 ggplot geom_violin geom_boxplot geom_jitter guides stat_summary
#' @importFrom cowplot theme_cowplot
#' @importFrom purrr flatten_chr flatten_dbl
#' @importFrom tibble tibble
#' @importFrom dplyr case_when
#'
#' @return a ggplot2 boxplot of cluster-wise confidence
#'
#' @export
PlotConfidenceDistributions <- function(obj,
                                        consensusName="jackknifeConsensus",
                                        confidenceName="jackknifeConfidence",
                                        boxplot=TRUE,
                                        violin=FALSE,
                                        jitter=FALSE,
                                        plotTitle=NULL) {
  tbl <- tibble(Consensus=obj[[consensusName]][,1],
                Confidence=obj[[confidenceName]][,1])
  p <- ggplot(tbl, aes(x=Consensus, y=Confidence, fill=Consensus))
  if(violin) {
    p <- p + geom_violin(scale="width")
  }
  if(boxplot) {
    p <- p + {
      if(violin) geom_boxplot(width=0.1, outlier.size=0, fill="white")
      else       geom_boxplot(outlier.size=0.5)
    }
  }
  if(jitter) {
    if(boxplot) { stop("Boxplot and jitter do not make sense together.") }
    p <- p + geom_jitter(size=0.5)
  }
  # if(annotateSize) {
  #   p <- p + stat_summary(fun.data=.stat_n, geom="text", size=3)
  # }
  p <- p  + guides(fill="none") + theme_cowplot()
  if(!is.null(plotTitle)) {
    p <- p + labs(title=plotTitle)
  }

  p

}


# Modified Seurat code to include confidence information

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dimensional reduction plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Get the names of objects within a Seurat object that are of a certain class
#
# @param object A Seurat object
# @param classes.keep A vector of names of classes to get
#
# @return A vector with the names of objects within the Seurat object that are of class \code{classes.keep}
#
.FilterObjects <- function(object, classes.keep = c('Assay', 'DimReduc')) {
  slots <- na.omit(object = Filter(
    f = function(x) {
      return(class(x = slot(object = object, name = x)) == 'list')
    },
    x = slotNames(x = object)
  ))
  slots <- grep(pattern = 'tools', x = slots, value = TRUE, invert = TRUE)
  slots <- grep(pattern = 'misc', x = slots, value = TRUE, invert = TRUE)
  slots.objects <- unlist(
    x = lapply(
      X = slots,
      FUN = function(x) {
        return(names(x = slot(object = object, name = x)))
      }
    ),
    use.names = FALSE
  )
  object.classes <- sapply(
    X = slots.objects,
    FUN = function(i) {
      return(class(x = object[[i]]))
    }
  )
  object.classes <- object.classes[object.classes %in% classes.keep]
  return(names(x = object.classes))
}

# Find the default DimReduc
#
# Searches for DimReducs matching 'umap', 'tsne', or 'pca', case-insensitive, and
# in that order. Priority given to DimReducs matching the DefaultAssay or assay specified
# (eg. 'pca' for the default assay weights higher than 'umap' for a non-default assay)
#
# @param object A Seurat object
# @param assay Name of assay to use; defaults to the default assay of the object
#
# @return The default DimReduc, if possible
#
.DefaultDimReduc <- function(object, assay = NULL) {
  assay <- assay %||% DefaultAssay(object = object)
  drs.use <- c('umap', 'tsne', 'pca')
  dim.reducs <- .FilterObjects(object = object, classes.keep = 'DimReduc')
  drs.assay <- Filter(
    f = function(x) {
      return(DefaultAssay(object = object[[x]]) == assay)
    },
    x = dim.reducs
  )
  if (length(x = drs.assay) > 0) {
    index <- lapply(
      X = drs.use,
      FUN = grep,
      x = drs.assay,
      ignore.case = TRUE
    )
    index <- Filter(f = length, x = index)
    if (length(x = index) > 0) {
      return(drs.assay[min(index[[1]])])
    }
  }
  index <- lapply(
    X = drs.use,
    FUN = grep,
    x = dim.reducs,
    ignore.case = TRUE
  )
  index <- Filter(f = length, x = index)
  if (length(x = index) < 1) {
    stop(
      "Unable to find a DimReduc matching one of '",
      paste(drs.use[1:(length(x = drs.use) - 1)], collapse = "', '"),
      "', or '",
      drs.use[length(x = drs.use)],
      "', please specify a dimensional reduction to use",
      call. = FALSE
    )
  }
  return(dim.reducs[min(index[[1]])])
}

# Make a theme for facet plots
#
# @inheritParams SeuratTheme
# @export
#
# @rdname SeuratTheme
# @aliases FacetTheme
#
.FacetTheme <- function(...) {
  return(theme(
    strip.background = element_blank(),
    strip.text = element_text(face = 'bold'),
    # Validate the theme
    validate = TRUE,
    ...
  ))
}



#' Dimensional reduction plot
#'
#' Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a
#' cell and it's positioned based on the cell embeddings determined by the reduction technique. By
#' default, cells are colored by their identity class (can be changed with the group.by parameter).
#'
#' @param object Seurat object
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param cells Vector of cells to plot (default is all cells)
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
#' or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
#' By default, ggplot2 assigns colors
#' @param pt.size Adjust point size for plotting
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param confidence One of the confidence values to plot (e.g. jackknifeConfidence). Assumes this is attached to the
#' Seurat object.
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by;
#' see \code{\link{FetchData}} for more details
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param order Specify the order of plotting for the idents. This can be
#' useful for crowded plots if points of interest are being buried. Provide
#' either a full list of valid idents or a subset to be plotted last (on top)
#' @param label Whether to label the clusters
#' @param label.size Sets size of labels
#' @param repel Repel labels
#' @param cells.highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply
#' pass a vector instead of a list. If set, colors selected cells to the color(s)
#' in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#' will also resize to the size(s) passed to \code{sizes.highlight}
#' @param cols.highlight A vector of colors to highlight the cells as; will
#' repeat to the length groups in cells.highlight
#' @param sizes.highlight Size of highlighted cells; will repeat to the length
#' groups in cells.highlight
#' @param na.value Color value for NA points when using custom scale
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param ncol Number of columns for display when combining plots
#' @param ... Extra parameters passed on to \code{\link{CombinePlots}}
#'
#' @return A ggplot object
#'
#' @importFrom rlang !! "%||%"
#' @importFrom ggplot2 facet_wrap vars sym
#' @importFrom Seurat Embeddings LabelClusters CombinePlots "Idents<-"
#'
#' @export
#'
#' @note For the old \code{do.hover} and \code{do.identify} functionality, please see
#' \code{HoverLocator} and \code{CellSelector}, respectively.
#'
#' @aliases TSNEPlot PCAPlot ICAPlot
#' @seealso \code{\link{FeaturePlot}} \code{\link{HoverLocator}}
#' \code{\link{CellSelector}} \code{link{FetchData}}
#'
#' @examples
#' DimPlot(object = pbmc_small)
#' DimPlot(object = pbmc_small, split.by = 'ident')
#'
DimConfidencePlot <- function(
  object,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = NULL,
  confidence = NULL,
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  order = NULL,
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  combine = TRUE,
  ncol = NULL,
  ...
) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  reduction <- reduction %||% .DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[['ident']] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  if (!is.null(x = confidence)) {
    data[, confidence] <- object[[confidence]]
  }
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      plot <- .SingleDimPlot(
        data = data[, c(dims, x, split.by, shape.by, confidence)],
        dims = dims,
        col.by = x,
        cols = cols,
        confidence = confidence,
        pt.size = pt.size,
        shape.by = shape.by,
        order = order,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value
      )
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = x,
          repel = repel,
          size = label.size,
          split.by = split.by
        )
      }
      if (!is.null(x = split.by)) {
        plot <- plot + .FacetTheme() +
          facet_wrap(
            facets = vars(!!sym(x = split.by)),
            ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split.by]))
            } else {
              ncol
            }
          )
      }
      return(plot)
    }
  )
  if (combine) {
    plots <- CombinePlots(
      plots = plots,
      ncol = if (!is.null(x = split.by) && length(x = group.by) > 1) {
        1
      } else {
        ncol
      },
      ...
    )
  }
  return(plots)
}


# Automagically calculate a point size for ggplot2-based scatter plots
#
# It happens to look good
#
# @param data A data frame being passed to ggplot2
#
# @return The "optimal" point size for visualizing these data
#
# @examples
# df <- data.frame(x = rnorm(n = 10000), y = runif(n = 10000))
# AutoPointSize(data = df)
#
.AutoPointSize <- function(data) {
     return(min(1583 / nrow(x = data), 1))
}

# Plot a single dimension
#
# @param data Data to plot
# @param dims A two-length numeric vector with dimensions to use
# @param pt.size Adjust point size for plotting
# @param col.by ...
# @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
# or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
# By default, ggplot2 assigns colors
# @param shape.by If NULL, all points are circles (default). You can specify any cell attribute
# (that can be pulled with FetchData) allowing for both different colors and different shapes on
# cells.
# @param order Specify the order of plotting for the idents. This can be useful for crowded plots if
# points of interest are being buried. Provide either a full list of valid idents or a subset to be
# @param confidence One of the confidence values to plot (e.g. jackknifeConfidence). Assumes this is attached to the
# Seurat object.
# plotted last (on top).
# @param label Whether to label the clusters
# @param repel Repel labels
# @param label.size Sets size of labels
# @param cells.highlight A list of character or numeric vectors of cells to
# highlight. If only one group of cells desired, can simply
# pass a vector instead of a list. If set, colors selected cells to the color(s)
# in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#  will also resize to the size(s) passed to \code{sizes.highlight}
# @param cols.highlight A vector of colors to highlight the cells as; will
# repeat to the length groups in cells.highlight
# @param sizes.highlight Size of highlighted cells; will repeat to the length
# groups in cells.highlight
# @param na.value Color value for NA points when using custom scale.
#
#' @importFrom cowplot theme_cowplot
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 ggplot aes_string labs geom_text guides geom_point
#' scale_color_brewer scale_color_manual element_rect guide_legend
#'
.SingleDimPlot <- function(
  data,
  dims,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  shape.by = NULL,
  order = NULL,
  confidence = NULL,
  label = FALSE,
  repel = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50'
) {
  pt.size <- pt.size %||% .AutoPointSize(data = data)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = sizes.highlight %||% pt.size,
      cols.highlight = cols.highlight,
      col.base = cols[1] %||% '#C3C3C3',
      pt.size = pt.size
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(
        order,
        setdiff(x = unique(x = data[, col.by]), y = order)
      ))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    # col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = '^\\d', x = col.by)) {
      # Do something for numbers
      col.by <- paste0('x', col.by)
    } else if (grepl(pattern = '-', x = col.by)) {
      # Do something for dashes
      col.by <- gsub(pattern = '-', replacement = '.', x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  aes_params = aes_string(
    x = dims[1],
    y = dims[2],
    color = paste0("`", col.by, "`"),
    shape = shape.by
  )
  if(!is.null(confidence)) {
    aes_params <- c(aes_params, aes_string(alpha=confidence))
    class(aes_params) <- "uneval"
  }
  plot <- ggplot(data = data) +
    geom_point(
      mapping = aes_params,
      size = pt.size
    ) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(color = NULL)
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot,
      id = col.by,
      repel = repel,
      size = label.size
    )
  }
  if (!is.null(x = cols)) {
    plot <- plot + if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale_color_brewer(palette = cols, na.value = na.value)
    } else {
      scale_color_manual(values = cols, na.value = na.value)
    }
  }
  plot <- plot + theme_cowplot()
  return(plot)
}
