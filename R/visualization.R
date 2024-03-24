#' @importFrom ggplot2 ggproto
NULL

####
# Spatial plots ####
####

####
## Spatial Identity Plot ####
####

#' vrSpatialPlot
#'
#' Plotting identification of spatially resolved cells, spots, and ROI on associated images from multiple assays in a VoltRon object.
#'
#' @param object VoltRon object
#' @param group.by a grouping label for the spatial entities
#' @param plot.segments plot segments instead of points
#' @param group.ids a subset of categories defined with in the grouping label \code{group.by}
#' @param n.tile should points be aggregated into tiles before visualization (see \code{geom_tile}). Applicable only for cells and molecules
#' @param assay the assay name
#' @param graph.name if not NULL, the spatial graph is with name \code{graph.name} is visualized as well
#' @param reduction Used by \code{vrSpatialPlotInteractive}, to visualize an embedding alongside with the spatial plot.
#' @param ncol column wise number of plots, for \code{ggarrange}
#' @param nrow row wise number of plots, for \code{ggarrange}
#' @param font.size font sizes
#' @param pt.size point size
#' @param cell.shape the shape of the points representing cells, see \code{help(geom_point)}
#' @param alpha alpha level for cells/spots/ROIs
#' @param label if TRUE, the labels of the ROI assays will be visualized
#' @param background the background of the plot, either "image" for overlaying the image of the assay, or "black" or "white" background (suitable for IF based assays)
#' @param reg if TRUE, the registered coordinates will be used
#' @param crop whether to crop an image of a spot assay
#' @param legend.pt.size the size of points at the legend
#' @param scale.image should the image be scaled down to a low resolution (width: 1000px)
#' @param legend.loc the location of the legend, default is "right"
#' @param common.legend whether to use a common legend for all plots
#' @param collapse whether to combine all ggplots
#' @param interactive if TRUE, run interactive plot
#'
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#'
#' @export
#'
vrSpatialPlot <- function(object, group.by = "Sample", plot.segments = FALSE, group.ids = NULL, n.tile = 0, assay = NULL, graph.name = NULL,
                          reduction = "umap", ncol = 2, nrow = NULL,
                          font.size = 2, pt.size = 2, cell.shape = 21, alpha = 1, label = FALSE, background = NULL, reg = FALSE,
                          crop = FALSE, legend.pt.size = 2, scale.image = TRUE, legend.loc = "right", common.legend = TRUE, collapse = TRUE, interactive = FALSE) {

  # check object for zarr
  if(is.character(object)){
    if(grepl(".zarr$", object)){

      return(vrSpatialPlotVitessce(zarr.file = object, group.by = group.by, plot.segments = plot.segments, group.ids = group.ids, assay = assay,
                                      reduction = reduction, background = background, reg = reg,  crop = crop))
    }
  }

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # interactive plotting
  if(interactive){
    if(length(assay_names) > 1){
      stop("Only one assay can be visualized with the interactive plot")
    } else{
      gg <- vrSpatialPlot(object, group.by = group.by, plot.segments = plot.segments, group.ids = group.ids, n.tile = n.tile, assay = assay,
                          graph.name = graph.name, reduction = reduction, ncol = ncol, nrow = nrow, font.size = font.size, pt.size = pt.size,
                          cell.shape = cell.shape, alpha = alpha, label = label, background = background, reg = reg,
                          crop = crop, legend.pt.size = legend.pt.size, scale.image = FALSE, legend.loc = legend.loc, common.legend = common.legend, collapse = collapse,
                          interactive = FALSE)
      return(vrSpatialPlotInteractive(plot_g = gg))
    }
  }

  # get entity type and metadata
  metadata <- Metadata(object, assay = assay)

  # configure titles
  plot_title <- as.list(apply(sample.metadata[assay_names,], 1, function(x) paste(x["Sample"], x["Layer"], sep = ", ")))
  names(plot_title) <- assay_names

  # for each assay
  i <- 1
  gg <- list()
  for(assy in assay_names){

    # get assay
    cur_assay <- object[[assy]]
    if(inherits(metadata, 'data.table')){
      cur_metadata <- metadata[assay_id == assy,]
    } else {
      assy_id <- paste0(assy,"$")
      cur_metadata <- metadata[grepl(assy_id, rownames(metadata)),]
    }

    # get graph
    if(!is.null(graph.name)){
      if(graph.name %in% vrGraphNames(object)){
        graph <- vrGraph(object, assay = assy, graph.type = graph.name)
      } else {
        stop("the graph with named '", graph.name, "' was not found in the graph list!")
      }
    } else {
      graph <- NULL
    }

    # check group.by
    if(!group.by %in% colnames(metadata))
      stop("The column '", group.by, "' was not found in the metadata!")
    levels_group.by <- as.character(unique(metadata[[group.by]][!is.na(metadata[[group.by]])]))
    if(all(!is.na(as.numeric(levels_group.by)))){
      levels_group.by <- sort(as.numeric(levels_group.by))
    }
    cur_metadata[[group.by]] <- factor(cur_metadata[[group.by]], levels = levels_group.by)

    # adjust group.ids
    if(is.null(group.ids)){
      group.ids <- unique(metadata[[group.by]])
    }

    # check group.id
    levels_group.ids <- as.character(group.ids)
    if(all(!is.na(as.numeric(levels_group.ids)))){
      levels_group.ids <- sort(as.numeric(levels_group.ids))
    }
    # group.ids <- levels(factor(group.ids, levels = levels_group.ids))
    group.ids <- factor(group.ids, levels = levels_group.ids)

    # visualize
    p_title <- plot_title[[assy]]
    gg[[i]] <- vrSpatialPlotSingle(assay = cur_assay, metadata = cur_metadata,
                                   group.by = group.by, plot.segments = plot.segments, group.ids = group.ids, n.tile = n.tile, graph = graph, font.size = font.size, pt.size = pt.size,
                                   alpha = alpha, cell.shape = cell.shape, plot_title = p_title, background = background, reg = reg,
                                   crop = crop, legend.pt.size = legend.pt.size, scale.image = scale.image)
    i <- i + 1
  }

  # return a list of plots or a single one
  if(collapse){
    if(length(assay_names) > 1){
      if(length(gg) < ncol) ncol <- length(gg)
      return(ggpubr::ggarrange(plotlist = gg, ncol = ncol, nrow = ceiling(length(gg)/ncol), common.legend = common.legend, legend = legend.loc))
    } else {
      return(gg[[1]])
    }
  } else {
    gg
  }
}

#' vrSpatialPlotSingle
#'
#' Plotting a single assay from a VoltRon object. We plot the identification of spatially resolved cells, spots, and ROI on associated images.
#'
#' @param assay vrAssay object
#' @param metadata the metadata associated with the assay
#' @param group.by a grouping label for the spatial entities
#' @param plot.segments plot segments instead of points
#' @param group.ids a subset of categories defined with in the grouping label \code{group.by}
#' @param n.tile should points be aggregated into tiles before visualization (see \code{geom_tile}). Applicable only for cells and molecules
#' @param graph if not NULL, the graph is added to the plot
#' @param font.size font sizes
#' @param pt.size point size
#' @param cell.shape the shape of the points representing cells, see \code{help(geom_point)}
#' @param alpha alpha level for cells/spots/ROIs
#' @param plot_title the title of the single plot
#' @param background the background of the plot, either "image" for overlaying the image of the assay, or "black" or "white" background (suitable for IF based assays)
#' @param reg if TRUE, the registered coordinates will be used
#' @param crop whether to crop an image of a spot assay
#' @param legend.pt.size the size of points at the legend
#' @param scale.image should the image be scaled down to a low resolution (width: 1000px)
#'
#' @import ggplot2
#' @importFrom igraph get.data.frame
#'
#' @noRd
vrSpatialPlotSingle <- function(assay, metadata, group.by = "Sample", plot.segments = FALSE, group.ids = NULL, n.tile = 0, graph = NULL,
                                font.size = 2, pt.size = 2, cell.shape = 21, alpha = 1, plot_title = NULL, background = NULL,
                                reg = FALSE, crop = FALSE, legend.pt.size = 2, scale.image = TRUE){

  # plot
  g <- ggplot()
  scale_factors <- 1

  # add image
  if(is.null(background))
    background <- vrMainImage(assay)
  if(length(background) == 2) {
    channel <- background[2]
  } else {
    channel <- NULL
  }
  background <- background[1]
  if(background %in% vrImageNames(assay)){
    image <- vrImages(assay, name = background, channel = channel)
    if(!is.null(image)){
      info <- image_info(image)
      if(info$width > 1000 && scale.image){
        image <- magick::image_resize(image, geometry = "1000x")
        scale_factors <- info$width/1000
        info <- magick::image_info(image)
      }
      g <- g +
        ggplot2::annotation_raster(image, 0, info$width, info$height, 0, interpolate = FALSE)
    } else {
      info <- NULL
    }
    image_name <- background
  } else {
    info <- NULL
    image_name <- vrMainImage(assay)
  }

  # data
  coords <- as.data.frame(vrCoordinates(assay, image_name = image_name, reg = reg))
  coords <- coords/scale_factors
  segments <- vrSegments(assay, image_name = image_name)

  # plotting features
  if(!group.by %in% colnames(metadata))
    stop("The column '", group.by, "' was not found in the metadata!")
  if(inherits(metadata, "data.table")){
    coords[[group.by]] <- metadata[,get(names(metadata)[which(colnames(metadata) == group.by)])]
  } else {
    coords[[group.by]] <- metadata[,group.by]
  }
  if(!is.null(group.ids)){
    if(length(setdiff(group.ids,  coords[[group.by]])) > 0){
      # warning("Some groups defined in group.ids does not exist in group.by!")
      coords <- droplevels(coords[coords[[group.by]] %in% group.ids,])
    } else if(length(setdiff(group.ids,  coords[[group.by]])) == length(group.ids)){
      stop("None of the groups defined in group.ids exist in group.by!")
    } else {
      segments <- segments[coords[[group.by]] %in% group.ids]
      coords <- droplevels(coords[coords[[group.by]] %in% group.ids,])
    }
  }

  # change levels of groups
  coords[[group.by]] <- factor(coords[[group.by]], levels = group.ids)
  
  # set up the limits
  if(vrAssayTypes(assay) == "spot"){
    if(crop){
      g <- g +
        coord_fixed(xlim = range(coords$x), ylim = range(coords$y))
    } else {
      if(!is.null(info)){
        g <- g +
          coord_fixed(xlim = c(0,info$width), ylim = c(0,info$height))
      }
    }
  } else {
    if(crop){
      g <- g +
        coord_fixed(xlim = range(coords$x), ylim = range(coords$y))
    } else {
      if(!is.null(info)){
        g <- g +
          xlim(0,info$width) + ylim(0, info$height)
      }
    }
  }

  # visualize based on points type
  if(vrAssayTypes(assay) == "ROI"){
    polygon_data <- NULL
    circle_data <- NULL
    for(i in 1:length(segments)){
      cur_data <- as.data.frame(cbind(segments[[i]][,c("x","y")], names(segments)[i], coords[[group.by]][i]))
      if(nrow(segments[[i]]) > 1){
        colnames(cur_data) <- c("x", "y", "segment", "group.by")
        cur_data[,c("x", "y")] <- cur_data[,c("x", "y")]/scale_factors
        polygon_data <- as.data.frame(rbind(polygon_data, cur_data))
      } else {
        colnames(cur_data) <- c("x", "y", "rx", "ry", "segment", "group.by")
        cur_data[,c("x", "y","rx", "ry")] <- cur_data[,c("x", "y","rx", "ry")]/scale_factors
        circle_data <- as.data.frame(rbind(circle_data,  cur_data))
      }
    }
    if(!is.null(polygon_data)){
      g <- g +
        geom_polygon(aes(x = x, y = y, fill = group.by, group = segment), data = polygon_data, alpha = alpha)
    }
    if(!is.null(circle_data)){
      g <- g +
        ggforce::geom_ellipse(aes(x0 = as.numeric(x), y0 = as.numeric(y), a = as.numeric(rx), b = as.numeric(ry), angle = 0,
                                  fill = group.by, group = segment), data = circle_data, lwd = 0, alpha = alpha)
    }
    g <- g +
      scale_fill_manual(values = scales::hue_pal()(length(levels(coords[[group.by]]))), labels = levels(coords[[group.by]]), drop = FALSE) +
      guides(fill = guide_legend(title = group.by))

  # spot visualization
  } else if(vrAssayTypes(assay) == "spot"){
    g <- g +
      geom_spot(mapping = aes_string(x = "x", y = "y", fill = group.by), coords, shape = 21, alpha = alpha, spot.radius = vrAssayParams(assay, param = "vis.spot.radius")) +
      scale_fill_manual(values = scales::hue_pal()(length(levels(coords[[group.by]]))), labels = levels(coords[[group.by]]), drop = FALSE) +
      guides(fill = guide_legend(override.aes=list(shape = 21, size = 4, lwd = 0.1)))

  # cell visualization
  } else if(vrAssayTypes(assay) %in% c("cell", "tile")) {

      if(plot.segments){

        if(length(segments) == 0) {
          stop("No Segments are available in this assay!")
        } else {
          polygon_data <- do.call(rbind,segments)
          polygon_data[,c("x", "y")] <- polygon_data[,c("x", "y")]/scale_factors
          len_segments <- sapply(segments, nrow, simplify = TRUE)
          polygon_data <- data.frame(polygon_data, segment = rep(names(segments), len_segments), group.by = rep(coords[[group.by]], len_segments))
          g <- g +
            geom_polygon(aes(x = x, y = y, fill = group.by, group = segment), data = polygon_data, alpha = alpha) +
            scale_fill_manual(values = scales::hue_pal()(length(levels(coords[[group.by]]))), labels = levels(coords[[group.by]]), drop = FALSE) +
            guides(fill = guide_legend(title = group.by))
        }
      } else {

        # add points
        if(n.tile == 0){
          g <- g +
            geom_point(mapping = aes_string(x = "x", y = "y", fill = group.by, color = group.by), coords, shape = cell.shape, size = rel(pt.size), alpha = alpha) +
            scale_fill_manual(values = scales::hue_pal()(length(levels(coords[[group.by]]))), labels = levels(coords[[group.by]]), drop = FALSE) +
            scale_color_manual(values = scales::hue_pal()(length(levels(coords[[group.by]]))), labels = levels(coords[[group.by]]), drop = FALSE) +
            guides(color = guide_legend(override.aes=list(size = legend.pt.size)))
        } else {
          g <- vrSpatialPlotSingleTiling(g = g, data = coords, n.tile = n.tile, alpha = alpha)
        }

        g <- g +
          guides(color = guide_legend(override.aes=list(size = legend.pt.size)))

        # add if a graph exists
        if(!is.null(graph)){
          graph.df <- igraph::get.data.frame(graph)
          graph.df$from.x <- coords$x[match(graph.df$from, rownames(coords))]
          graph.df$from.y <- coords$y[match(graph.df$from, rownames(coords))]
          graph.df$to.x <- coords$x[match(graph.df$to, rownames(coords))]
          graph.df$to.y <- coords$y[match(graph.df$to, rownames(coords))]
          g <- g +
            geom_segment(data = graph.df, mapping = aes(x=from.x,xend = to.x, y=from.y,yend = to.y), alpha = 0.5, color = ifelse(background == "black", "grey", "black"))
        }
      }
  } else if(vrAssayTypes(assay) == "molecule") {

    if(n.tile == 0){
      g <- g +
        geom_point(mapping = aes_string(x = "x", y = "y", fill = group.by, color = group.by), coords, shape = cell.shape, size = rel(pt.size), alpha = alpha) +
        scale_fill_manual(values = scales::hue_pal()(length(levels(coords[[group.by]]))), labels = levels(coords[[group.by]]), drop = FALSE) +
        scale_color_manual(values = scales::hue_pal()(length(levels(coords[[group.by]]))), labels = levels(coords[[group.by]]), drop = FALSE) +
        guides(color = guide_legend(override.aes=list(size = legend.pt.size)))
    } else {
      # coords_orig <- as.data.frame(vrCoordinates(assay, image_name = image_name, reg = reg))
      # coords_orig <- coords_orig/scale_factors
      # coords_orig[[group.by]] <- NA
      # coords_orig[rownames(coords), group.by] <- coords[[group.by]]
      g <- vrSpatialPlotSingleTiling(g = g, data = coords, n.tile = n.tile, alpha = alpha)
    }

  } else {
    stop("Only ROIs, spots, cells, molecules and tiles can be visualized with vrSpatialPlot!")
  }

  # more visualization parameters
  g <- g +
    ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,0,0)),
                                panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                                axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                                legend.margin = margin(0,0,0,0))

  # # set up the limits
  # if(vrAssayTypes(assay) == "spot"){
  #   if(crop){
  #     g <- g +
  #       coord_fixed(xlim = range(coords$x), ylim = range(coords$y))
  #   } else {
  #     if(!is.null(info)){
  #       g <- g +
  #         coord_fixed(xlim = c(0,info$width), ylim = c(0,info$height))
  #     }
  #   }
  # } else {
  #   if(crop){
  #     g <- g +
  #       coord_fixed(xlim = range(coords$x), ylim = range(coords$y))
  #   } else {
  #     if(!is.null(info)){
  #       g <- g +
  #         xlim(0,info$width) + ylim(0, info$height)
  #     }
  #   }
  # }

  # background
  if(any(background %in% c("white","black"))){
    g <- g +
      theme(panel.background = element_rect(fill = background, colour = background, size = 0.5, linetype = "solid"))
  } else{
    if(is.null(info)){
      g <- g +
        theme(panel.background = element_rect(fill = "lightgrey", colour = "lightgrey", size = 0.5, linetype = "solid"))
    } else{
      g <- g +
        theme(panel.background = element_blank())
    }
  }

  # return data
  return(g)
}

#' vrSpatialPlotSingleTiling
#'
#' Plotting a tiled version of the vrSpatialPlot
#'
#' @param g the ggplot figure
#' @param data the data frame with coordinates and group identities
#' @param n.tile should points be aggregated into tiles before visualization (see \code{geom_tile}). Applicable only for cells and molecules
#' @param alpha alpha level for cells/spots/ROIs
#'
#' @import ggplot2
#'
#' @noRd
vrSpatialPlotSingleTiling <- function(g, data, n.tile, alpha = 1){

  # gplot <- g + geom_hex(data = data, mapping = aes(x = x, y = y), bins = n.tile, alpha = alpha)
  gplot <- g + stat_bin_2d(mapping = aes(x = x, y = y), data = data, bins = n.tile, drop = FALSE, alpha = alpha)
  hex_count_data <- ggplot_build(gplot)$data
  hex_count_data <- hex_count_data[[length(hex_count_data)]]
  midpoint <- max(hex_count_data$count)/2
  gplot <- gplot +
    scale_fill_gradientn(name = "Count",
                         colors=c("dodgerblue2", "white", "yellow3"),
                         values=scales::rescale(c(0, midpoint, max(hex_count_data$count))), limits = c(0, max(hex_count_data$count)))

  # return
  gplot
}

####
## Spatial Feature Plot ####
####

#' vrSpatialFeaturePlot
#'
#' Plotting single/multiple features of spatially resolved cells, spots, and ROI on associated images from multiple assays in a VoltRon object.
#'
#' @param object VoltRon object
#' @param features a set of features, either from the rows of rawdata, normdata or columns of the metadata
#' @param group.by a grouping label for the spatial entities
#' @param plot.segments plot segments instead of points
#' @param norm if TRUE, the normalized data is used
#' @param log if TRUE, data features (excluding metadata features) will be log transformed
#' @param assay the assay name
#' @param graph.name if not NULL, the spatial graph is with name \code{graph.name} is visualized as well
#' @param ncol column wise number of plots, for \code{ggarrange}
#' @param nrow row wise number of plots, for \code{ggarrange}
#' @param font.size font sizes
#' @param pt.size point size
#' @param title.size title size of legend and plot
#' @param alpha alpha level for cells/spots/ROIs
#' @param keep.scale whether unify all scales for all features or not
#' @param label if TRUE, labels of ROIs will be visualized too
#' @param background the background of the plot, either "image" for overlaying the image of the assay, or "black" or "white" background (suitable for IF based assays)
#' @param reg if TRUE, the registered coordinates will be used
#' @param crop whether to crop an image of a spot assay
#' @param common.legend whether to use a common legend for all plots
#' @param collapse whether to combine all ggplots
#'
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#'
#' @export
#'
vrSpatialFeaturePlot <- function(object, features, group.by = "label", plot.segments = FALSE, norm = TRUE, log = FALSE, assay = NULL, graph.name = NULL, ncol = 2, nrow = NULL,
                         font.size = 2, pt.size = 2, title.size = 10, alpha = 0.6, keep.scale = "feature", label = FALSE, background = NULL, reg = FALSE,
                         crop = FALSE, common.legend = FALSE, collapse = TRUE) {

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get entity type and metadata
  metadata <- Metadata(object, assay = assay)

  # calculate limits for plotting, all for making one scale, feature for making multiple
  limits <- Map(function(feat){
    range_feat <- Map(function(assy){
      spatialpoints <- vrSpatialPoints(object[[assy]])
      if(feat %in% vrFeatures(object, assay = assy)){
        data <- vrData(object[[assy]], features = feat, norm = norm)
        if(log)
          data <- log(data)
        return(range(data, na.rm = TRUE, finite = TRUE))
      } else {
        metadata <- Metadata(object, assay = assy)
        if(feat %in% colnames(metadata)){
          return(range(metadata[,feat], na.rm = TRUE, finite = TRUE))
        } else {
          stop("Feature '", feat, "' cannot be found in data or metadata!")
        }
      }
    }, assay_names)
    if(keep.scale == "all"){
      range_feat_all <- c(do.call(min, range_feat), do.call(max, range_feat))
      range_feat <- Map(function(assy) return(range_feat_all), assay_names)
    }
    if(!keep.scale %in% c("all", "feature")){
      stop("keep.scale should be either 'all' or 'feature', check help(vrSpatialFeaturePlot)")
    }
    return(range_feat)
  }, features)

  # configure titles
  assay_title <- as.list(apply(sample.metadata[assay_names,], 1, function(x) paste(x["Sample"], x["Layer"], sep = ", ")))
  names(assay_title) <- assay_names
  feature_title <- as.list(features)
  names(feature_title) <- features
  plot_title <- assay_title
  legend_title <- feature_title

  # for each feature
  i <- 1
  gg <- list()
  for(assy in assay_names){

    # for each assay
    for(feat in features){

      # get assay
      cur_assay <- object[[assy]]
      cur_metadata <- metadata[grepl(paste0(assy, "$"), rownames(metadata)),]

      # get graph
      if(!is.null(graph.name)){
        if(graph.name %in% vrGraphNames(object)){
          graph <- vrGraph(object, assay = assy, graph.type = graph.name)
        } else {
          stop("the graph with named '", graph.name, "' was not found in the graph list!")
        }
      } else {
        graph <- NULL
      }

      # visualize
      p_title <- plot_title[[assy]]
      l_title <- legend_title[[feat]]
      gg[[i]] <- vrSpatialFeaturePlotSingle(assay = cur_assay, metadata = cur_metadata, feature = feat, plot.segments = plot.segments, graph = graph, limits = limits[[feat]][[assy]],
                              group.by = group.by, norm = norm, log = log, font.size = font.size, pt.size = pt.size, title.size = title.size, alpha = alpha,
                              label = label, plot_title = p_title, legend_title = l_title, background = background, reg = reg, crop = crop)
      i <- i + 1
    }
  }

  if(collapse){
    # return a list of plots or a single one
    if(length(features) > 1 && length(assay_names) > 1){
      return(ggpubr::ggarrange(plotlist = gg, ncol = length(features), nrow = length(assay_names)))
    } else if(length(features) > 1 && length(assay_names) == 1){
      if(length(gg) < ncol) ncol <- length(gg)
      return(ggpubr::ggarrange(plotlist = gg, ncol = ncol, nrow = ceiling(length(gg)/ncol)))
    } else if(length(features) == 1 && length(assay_names) > 1){
      if(length(gg) < ncol) ncol <- length(gg)
      return(ggpubr::ggarrange(plotlist = gg, ncol = ncol, nrow = ceiling(length(gg)/ncol), common.legend = common.legend, legend = "right"))
    } else {
      return(gg[[1]])
    }
  } else {
    return(gg)
  }
}

#' vrSpatialFeaturePlotSingle
#'
#' A single Spatial Feature plot of VoltRon object
#'
#' @param assay vrAssay object
#' @param metadata the metadata associated with the assay
#' @param feature a feature, either from the rows of rawdata, normdata or columns of the metadata
#' @param plot.segments plot segments instead of points
#' @param graph if not NULL, the graph is added to the plot
#' @param limits limits of the legend of the plot
#' @param group.by a grouping label for the spatial entities
#' @param norm if TRUE, the normalized data is used
#' @param log if TRUE, data features (excluding metadata features) will be log transformed
#' @param font.size font sizes
#' @param pt.size point size
#' @param title.size title size of legend and plot
#' @param alpha alpha level for cells/spots/ROIs
#' @param label if TRUE, labels of ROIs will be visualized too
#' @param plot_title the main title of the single plot
#' @param legend_title the legend title of the single plot
#' @param background the background of the plot, either "image" for overlaying the image of the assay, or "black" or "white" background (suitable for IF based assays)
#' @param reg if TRUE, the registered coordinates will be used
#' @param crop whether to crop an image of a spot assay
#'
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggforce geom_ellipse
#' @importFrom igraph get.data.frame
#' @importFrom dplyr arrange
#'
#' @noRd
vrSpatialFeaturePlotSingle <- function(assay, metadata, feature, plot.segments = FALSE, graph = NULL, limits, group.by = "label", norm = TRUE, log = FALSE,
                               font.size = 2, pt.size = 2, title.size = 10, alpha = 0.6, label = FALSE, plot_title = NULL,
                               legend_title = NULL, background = NULL, reg = FALSE, crop = FALSE){

  # plot
  g <- ggplot()
  scale_factors <- 1

  # add image
  if(is.null(background))
    background <- vrMainImage(assay)
  if(length(background) == 2) {
    channel <- background[2]
  } else {
    channel <- NULL
  }
  background <- background[1]
  if(background %in% vrImageNames(assay)){
    image <- vrImages(assay, name = background, channel = channel)
    if(!is.null(image)){
      info <- image_info(image)
      if(info$width > 1000){
        image <- magick::image_resize(image, geometry = "1000x")
        scale_factors <- info$width/1000
        info <- magick::image_info(image)
      }
      g <- g +
        ggplot2::annotation_raster(image, 0, info$width, info$height, 0, interpolate = FALSE)
    } else {
      info <- NULL
    }
    image_name <- background
  } else {
    info <- NULL
    image_name <- vrMainImage(assay)
  }

  # data
  coords <- as.data.frame(vrCoordinates(assay, image_name = image_name, reg = reg))
  coords <- coords/scale_factors
  segments <- vrSegments(assay, image_name = image_name)
  data_features <- feature[feature %in% vrFeatures(assay)]
  if(length(data_features) > 0){
    normdata <- vrData(assay, features = feature, norm = norm)
    if(log)
      normdata <- log(normdata)
  }

  # get data
  if(feature %in% data_features){
    coords$score <- normdata[feature,]
  } else {
    coords$score <- metadata[,feature]
  }

  # get image information and plotting features
  midpoint <- sum(limits)/2

  # add points or segments
  if(vrAssayTypes(assay) == "ROI" && !is.null(segments)){
    polygon_data <- NULL
    circle_data <- NULL
    for(i in 1:length(segments)){
      cur_data <- as.data.frame(cbind(segments[[i]], names(segments)[i], coords$score[i]))
      if(nrow(segments[[i]]) > 1){
        colnames(cur_data) <- c("x", "y", "segment", "score")
        cur_data[,c("x", "y")] <- cur_data[,c("x", "y")]/scale_factors
        polygon_data <- as.data.frame(rbind(polygon_data, cur_data))
      } else {
        colnames(cur_data) <- c("x", "y", "rx", "ry", "segment", "score")
        cur_data[,c("x", "y","rx", "ry")] <- cur_data[,c("x", "y","rx", "ry")]/scale_factors
        circle_data <- as.data.frame(rbind(circle_data,  cur_data))
      }
    }
    if(!is.null(geom_polygon)){
      g <- g +
        geom_polygon(aes(x = x, y = y, fill = score, group = segment), data = polygon_data, alpha = alpha)
    }
    if(!is.null(circle_data)){
      g <- g +
        ggforce::geom_ellipse(aes(x0 = as.numeric(x), y0 = as.numeric(y), a = as.numeric(rx), b = as.numeric(ry), angle = 0,
                         fill = score, group = segment), data = circle_data, lwd = 0, alpha = alpha)
    }
    g <- g +
      scale_fill_gradientn(name = legend_title,
                           colors=c("dodgerblue2", "white", "yellow3"),
                           values=scales::rescale(c(limits[1], midpoint, limits[2])), limits = limits)
  } else if(vrAssayTypes(assay) == "spot"){
    g <- g +
      geom_spot(mapping = aes(x = x, y = y, fill = score), coords, shape = 21, alpha = alpha, spot.radius = vrAssayParams(assay, param = "vis.spot.radius")) +
      scale_fill_gradientn(name = legend_title,
                             colors=c("dodgerblue3", "yellow", "red"),
                             values=scales::rescale(c(limits[1], midpoint, limits[2])), limits = limits)
  } else if(vrAssayTypes(assay) %in% c("cell", "tile")) {

    if(plot.segments){

      if(length(segments) == 0) {
        stop("No Segments are available in this assay!")
      } else {
        polygon_data <- do.call(rbind,segments)
        polygon_data[,c("x", "y")] <- polygon_data[,c("x", "y")]/scale_factors
        len_segments <- sapply(segments, nrow, simplify = TRUE)
        polygon_data <- data.frame(polygon_data, segment = rep(names(segments), len_segments), score = rep(coords$score, len_segments))
        g <- g +
          geom_polygon(aes(x = x, y = y, fill = score, group = segment), data = polygon_data, alpha = alpha)
      }
      g <- g +
        scale_fill_gradientn(name = legend_title,
                             colors=c("dodgerblue2", "white", "yellow3"),
                             values=scales::rescale(c(limits[1], midpoint, limits[2])), limits = limits)
    } else {
      g <- g +
        geom_point(mapping = aes(x = x, y = y, colour = score), dplyr::arrange(coords,score), shape = 16, size = rel(pt.size), alpha = alpha) +
        scale_colour_gradientn(name = legend_title,
                               colors=c("dodgerblue2", "white", "yellow3"),
                               values=scales::rescale(c(limits[1], midpoint, limits[2])), limits = limits)

      # add if a graph exists
      if(!is.null(graph)){
        graph.df <- igraph::get.data.frame(graph)
        graph.df$from.x <- coords$x[match(graph.df$from, rownames(coords))]
        graph.df$from.y <- coords$y[match(graph.df$from, rownames(coords))]
        graph.df$to.x <- coords$x[match(graph.df$to, rownames(coords))]
        graph.df$to.y <- coords$y[match(graph.df$to, rownames(coords))]
        g <- g +
          geom_segment(data = graph.df, mapping = aes(x=from.x,xend = to.x, y=from.y,yend = to.y), alpha = 0.5, color = ifelse(background == "black", "grey", "black"))
      }
    }
  }

  # more visualization parameters
  g <- g +
    ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,0,0), size = title.size),
                                panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                                axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                                legend.key.size = unit(title.size, "points"), legend.title = element_text(size=title.size),
                                legend.margin = margin(0,0,0,0))

  # set up the limits
  if(vrAssayTypes(assay) == "spot"){
    if(crop){
      g <- g +
        coord_fixed(xlim = range(coords$x), ylim = range(coords$y))
    } else {
      if(!is.null(info)){
        g <- g +
          coord_fixed(xlim = c(0,info$width), ylim = c(0,info$height))
      }
    }
  } else {
    if(crop){
      g <- g +
        coord_fixed(xlim = range(coords$x), ylim = range(coords$y))
    } else {
      if(!is.null(info)){
        g <- g +
          xlim(0,info$width) + ylim(0, info$height)
      }
    }
  }

  # background
  # if(any(background %in% c("white","black"))){
  #   g <- g +
  #     theme(panel.background = element_rect(fill = background, colour = background, size = 0.5, linetype = "solid"))
  # } else if(background %in% vrImageNames(assay)){
  #   g <- g +
  #     theme(panel.background = element_blank())
  # } else {
  #   g <- g +
  #     theme(panel.background = element_rect(fill = "lightgrey", colour = "lightgrey", size = 0.5, linetype = "solid"))
  #   warning("background image ", background, " is not found in ", vrAssayNames(assay), "\n")
  # }
  if(any(background %in% c("white","black"))){
    g <- g +
      theme(panel.background = element_rect(fill = background, colour = background, size = 0.5, linetype = "solid"))
  } else{
    if(is.null(info)){
      g <- g +
        theme(panel.background = element_rect(fill = "lightgrey", colour = "lightgrey", size = 0.5, linetype = "solid"))
    } else{
      g <- g +
        theme(panel.background = element_blank())
    }
  }

  # visualize labels
  if(label){
    if(group.by %in% colnames(metadata)){
      coords[[group.by]] <- metadata[,group.by]
    } else {
      stop("The column ", group.by, " was not found in the metadata!")
    }
    g <- g + ggrepel::geom_label_repel(mapping = aes_string(x = "x", y = "y", label = group.by), coords,
                                box.padding = 0.5, size = font.size, direction = "both", seed = 1)
  }

  # return data
  return(g)
}

####
## Spatial Auxiliary ####
####

#' @import ggplot2
#' @importFrom grid pointsGrob unit gpar
#' @importFrom rlang list2
#'
#' @noRd
NULL

geom_spot <- function (mapping = NULL, data = NULL, stat = "identity", position = "identity",
                           ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSpot,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = rlang::list2(na.rm = na.rm, ...))
}

GeomSpot <- ggplot2::ggproto("GeomSpot",
                      Geom,
                      required_aes = c("x", "y"),
                      non_missing_aes = c("size", "shape", "colour"),
                      default_aes = aes(
                        shape = 21,
                        colour = "black",
                        size = 1.5,
                        fill = NA,
                        spot.radius = 1,
                        alpha = NA, stroke = 0.5
                      ),
                      draw_panel = function(self, data, panel_params, coord, na.rm = FALSE) {
                        if (is.character(data$shape)) {
                          data$shape <- translate_shape_string(data$shape)
                        }
                        coords <- coord$transform(data, panel_params)
                        stroke_size <- coords$stroke
                        stroke_size[is.na(stroke_size)] <- 0
                        xrange <- panel_params$x.range[2] - panel_params$x.range[1]
                        yrange <- panel_params$y.range[2] - panel_params$y.range[1]
                        mainrange <- min(xrange, yrange)
                        spot.radius <- data$spot.radius/mainrange
                        ggname("geom_spot",
                               grid::pointsGrob(
                                 coords$x, coords$y,
                                 pch = coords$shape,
                                 size = grid::unit(spot.radius, "npc"),
                                 gp = grid::gpar(
                                   col = ggplot2::alpha(coords$colour, coords$alpha),
                                   fill = ggplot2::alpha(coords$fill, coords$alpha),
                                   fontsize = coords$size * .pt + stroke_size * .stroke / 2,
                                   lwd = coords$stroke * .stroke / 2
                                 )
                               )
                        )
                      },

                      draw_key = draw_key_point
)

####
# Embeddings plots ####
####

####
## Embedding Identity Plot ####
####

#' vrEmbeddingPlot
#'
#' Plotting embeddings of cells and spots on associated images from multiple assays in a VoltRon object.
#'
#' @param object VoltRon object
#' @param embedding the embedding type, i.e. pca, umap etc.
#' @param group.by a grouping label for the spatial entities
#' @param assay the assay name
#' @param ncol column wise number of plots, for \code{ggarrange}
#' @param nrow row wise number of plots, for \code{ggarrange}
#' @param font.size font size of labels, if label is TRUE
#' @param pt.size point size
#' @param label if TRUE, the labels of the ROI assays will be visualized
#' @param common.legend whether to use a common legend for all plots
#' @param collapse whether to combine all ggplots
#'
#' @import ggplot2
#' @importFrom stats aggregate
#'
#' @export
#'
vrEmbeddingPlot <- function(object, embedding = "pca", group.by = "Sample", assay = NULL, ncol = 2, nrow = NULL,
                            font.size = 5, pt.size = 1, label = FALSE, common.legend = TRUE, collapse = TRUE) {

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get entity type and metadata
  metadata <- Metadata(object, assay = assay)

  # grep assays from metadata
  assy_id <- paste(paste0(assay_names,"$"), collapse = "|")
  if(inherits(metadata, "data.table")){
    metadata <- subset(metadata, subset = assay_id %in% assay_names)
  } else {
    assy_id <- paste(paste0(assay_names,"$"), collapse = "|")
    metadata <- metadata[grepl(assy_id, rownames(metadata)),]
  }

  # plotting features
  datax <- data.frame(vrEmbeddings(object, assay = assay_names, type = embedding))
  datax <- datax[,1:2]
  colnames(datax) <- c("x", "y")
  if(group.by %in% colnames(metadata)){
    if(inherits(metadata, "data.table")){
      datax[[group.by]] <- metadata[,get(names(metadata)[which(colnames(metadata) == group.by)])]
    } else {
      datax[[group.by]] <- as.factor(metadata[,group.by])
    }
  } else {
    stop("Column ", group.by, " cannot be found in metadata!")
  }

  # plot
  g <- ggplot()

  # add points or segments
  g <- g +
    geom_point(mapping = aes_string(x = "x", y = "y", color = group.by), datax, shape = 16, size = pt.size) +
    guides(color = guide_legend(override.aes=list(size = 2)))

  # more visualization parameters
  g <- g +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,0,0)),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          legend.margin = margin(0,0,0,0), panel.background = element_blank()) +
    xlab(paste0(toupper(embedding), "_1")) + ylab(paste0(toupper(embedding), "_2"))

  # labels
  if(label){
    datax_group <- stats::aggregate(datax[,c("x","y")], list(datax[[group.by]]),mean)
    colnames(datax_group) <- c(group.by, "x", "y")
    g <- g +
      geom_text(mapping = aes_string(x = "x", y = "y", label = group.by), datax_group, size = font.size)
  }

  # return data
  return(g)
}

####
## Embedding Feature Plot ####
####

#' vrEmbeddingFeaturePlot
#'
#' Plotting features of spatially resolved cells and spots on embeddings from multiple assays in a VoltRon object.
#'
#' @param object VoltRon object
#' @param embedding the embedding type, i.e. pca, umap etc.
#' @param features a set of features, either from the rows of rawdata, normdata or columns of the metadata
#' @param norm if TRUE, the normalized data is used
#' @param log if TRUE, data features (excluding metadata features) will be log transformed
#' @param assay the assay name
#' @param ncol column wise number of plots, for \code{ggarrange}
#' @param nrow row wise number of plots, for \code{ggarrange}
#' @param font.size font sizes
#' @param pt.size point size
#' @param keep.scale whether unify all scales for all features or not
#' @param common.legend whether to use a common legend for all plots
#' @param collapse whether to combine all ggplots
#'
#' @import ggplot2
#'
#' @export
#'
vrEmbeddingFeaturePlot <- function(object, embedding = "pca", features = NULL, norm = TRUE, log = FALSE, assay = NULL, ncol = 2, nrow = NULL,
                                   font.size = 2, pt.size = 1, keep.scale = "feature", common.legend = TRUE, collapse = TRUE) {

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # features
  if(is.null(features))
    stop("You have to define at least one feature")

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get entity type and metadata
  metadata <- Metadata(object, assay = assay)

  # get data
  data_features <- features[features %in% vrFeatures(object)]
  if(length(data_features) > 0){
    normdata <- vrData(object, features = data_features, assay = assay_names, norm = norm)
    if(log)
      normdata <- log(normdata)
  }

  # get embedding
  datax <- data.frame(vrEmbeddings(object, assay = assay_names, type = embedding))
  datax <- datax[,1:2]
  colnames(datax) <- c("x", "y")

  # calculate limits for plotting, all for making one scale, feature for making multiple
  limits <- Map(function(feat){
    if(feat %in% vrFeatures(object)){
      return(range(normdata[feat, ]))
    } else {
      if(feat %in% colnames(metadata)){
        return(range(metadata[,feat]))
      } else {
        stop("Feature '", feat, "' cannot be found in data or metadata!")
      }
    }
  }, features)

  # configure titles
  feature_title <- as.list(features)
  names(feature_title) <- features
  legend_title <- feature_title

  # for each feature
  i <- 1
  gg <- list()
  for(feat in features){

    # get data
    if(feat %in% vrFeatures(object)){
      datax$score <- normdata[feat,]
    } else {
      datax$score <- metadata[,feat]
    }

    # get image information and plotting features
    midpoint <- sum(limits[[feat]])/2

    # plot
    g <- ggplot()

    # add points or segments
    # datax <- datax[sample(1:nrow(datax)),]
    g <- g +
      geom_point(mapping = aes(x = x, y = y, color = score), dplyr::arrange(datax,score), shape = 16, size = pt.size) +
      scale_color_gradientn(name = legend_title[[feat]],
                            colors=c("lightgrey", "blue"),
                            values=scales::rescale(c(limits[[feat]][1], limits[[feat]][2])), limits = limits[[feat]])

    # more visualization parameters
    g <- g +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,0,0)),
            panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
            legend.margin = margin(0,0,0,0), panel.background = element_blank()) +
      xlab(paste0(toupper(embedding), "_1")) + ylab(paste0(toupper(embedding), "_2"))
    gg[[i]] <- g
    i <- i + 1
  }

  if(collapse){
    # return a list of plots or a single one
    if(length(features) > 1){
      if(length(gg) < ncol) ncol <- length(gg)
      return(ggpubr::ggarrange(plotlist = gg, ncol = ncol, nrow = ceiling(length(gg)/ncol)))
    } else {
      return(gg[[1]])
    }
  } else {
    return(gg)
  }
}

####
# Scatter Plot ####
####

#' vrScatterPlot
#'
#' get a scatter plot between two features
#'
#' @param object a VoltRon object
#' @param feature.1 first feature
#' @param feature.2 second feature
#' @param norm if TRUE, the normalize data will be used
#' @param assay assay name
#' @param pt.size point size
#' @param font.size font size
#' @param group.by a column from metadata to label points
#' @param label whether labels are visualized or not
#' @param trend inserting a trend line two the scatter plot
#'
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#'
#' @export
#'
vrScatterPlot <- function(object, feature.1, feature.2, norm = TRUE, assay = NULL,
                               pt.size = 2, font.size = 2, group.by = "label", label = FALSE, trend = FALSE){

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # check the number of features
  if(is.null(feature.1) | is.null(feature.2))
    stop("Please provide both 'feature.1' and 'feature.2'")

  # check the number of features
  if((length(feature.1) != 1 | length(feature.2) != 1))
    stop("Both 'feature.1' and 'feature.2' should be of length 1.")

  # data
  normdata <- vrData(object, assay = assay, norm = norm)

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get entity type and metadata
  # get entity type and metadata
  # if(is.null(assay.type)){
  #   assay_types <- vrAssayTypes(object, assay = assay)
  #   if(length(unique(assay_types)) == 1){
  #     metadata <- Metadata(object, type = unique(assay_types))
  #   } else {
  #     stop("Please select assay.type as 'cell', 'spot' or 'ROI'")
  #   }
  # } else {
  #   metadata <- Metadata(object, type = assay.type)
  # }
  metadata <- Metadata(object, assay = assay)

  # get data
  data_feature <- sapply(c(feature.1, feature.2), function(feat){
    if(feat %in% rownames(normdata)){
      return(normdata[feat,])
    } else {
      return(metadata[,feat])
    }
  })
  data_feature <- as.data.frame(data_feature)

  # plot
  g <- ggplot()

  # plot scatter
  g <- g +
    geom_point(mapping = aes_string(x = feature.1, y = feature.2), data = data_feature, size = pt.size)

  # visualize labels
  if(label){
    data_feature[[group.by]] <- metadata[,group.by]
    g <- g + ggrepel::geom_label_repel(mapping = aes_string(x = feature.1, y = feature.2, label = group.by), data_feature,
                                box.padding = 0.5, size = font.size, direction = "both", seed = 1)
  }

  # visualize trend
  if(trend){
    g <- g + geom_smooth()
  }

  g
}

####
# Heatmap Plot ####
####

#' HeatmapPlot
#'
#' @param object VoltRon object
#' @param assay assay name
#' @param features a set of features to be visualized
#' @param group.by a column from metadata to seperate columns of the heatmap
#' @param norm if TRUE, the normalized data is used
#' @param scaled if TRUE, the data will be scaled before visualization
#' @param show_row_names if TRUE, row names of the heatmap will be shown
#' @param cluster_rows if TRUE, the rows of the heatmap will be clustered
#' @param show_heatmap_legend if TRUE, the heatmap legend is shown
#' @param outlier.quantile quantile for detecting outliers whose values are set to the quantile, change to lower values to adjust large number of outliers, default: 0.99
#' @param highlight.some if TRUE, some rows will be showed at random, reproducible by \code{seed} arguement
#' @param n_highlight the number of shown row labels, if \code{show_row_names} is TRUE
#' @param font.size font size
#' @param seed the seed for \code{set.seed}
#' @param ... additional parameters passed to \code{getVariableFeatures}
#'
#' @importFrom scales viridis_pal
#' @importFrom stats quantile
#'
#' @export
#'
vrHeatmapPlot <- function(object, assay = NULL, features = NULL, group.by = "clusters",
                          norm = TRUE, scaled = TRUE, show_row_names = NULL, cluster_rows = TRUE, show_heatmap_legend = FALSE,
                          outlier.quantile = 0.99, highlight.some = FALSE, n_highlight = 30, font.size = 13.2, seed = 1, ...){

  if (!requireNamespace('ComplexHeatmap'))
    stop("Please install ComplexHeatmap package to use the Heatmap function")

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # seed
  set.seed(seed)

  # data
  heatmapdata <- vrData(object, assay = assay, norm = norm)

  # features
  if(is.null(features)){
    if(nrow(vrFeatureData(object)) > 0){
      features <- getVariableFeatures(object, assay = assay, ...)
    } else {
      features <- vrFeatures(object, assay = assay)
    }
  } else {
    nonmatching_features <- setdiff(features, vrFeatures(object))
    features <- features[features %in% vrFeatures(object)]
    if(length(features) == 0){
      stop("None of the provided features are found in the assay")
    }
    if(length(nonmatching_features))
      message("the following features are not found in the assay: ", paste(nonmatching_features, collapse = ", "))
  }
  heatmapdata <- heatmapdata[features, ]

  # get entity type and metadata
  # get entity type and metadata
  # if(is.null(assay.type)){
  #   assay_types <- vrAssayTypes(object, assay = assay)
  #   if(length(unique(assay_types)) == 1){
  #     metadata <- Metadata(object, type = unique(assay_types))
  #   } else {
  #     stop("Please select assay.type as 'cell', 'spot' or 'ROI'")
  #   }
  # } else {
  #   metadata <- Metadata(object, type = assay.type)
  # }
  metadata <- Metadata(object, assay = assay)
  metadata <- metadata[colnames(heatmapdata),]

  # scaling, optional
  if(scaled){
    heatmapdata <- apply(heatmapdata, 1, scale)
    heatmapdata <- t(heatmapdata)
    legend_title <- "Scaled \n Exp."
  } else {
    legend_title <- "Norm. \n Exp."
  }

  # manage data for plotting
  if(group.by %in% colnames(metadata)){
    heatmapdata <- heatmapdata[,order(metadata[[group.by]], decreasing = FALSE)]
    labels_ordered <- metadata[[group.by]][order(metadata[[group.by]], decreasing = FALSE)]
    labels_ordered_table <- table(labels_ordered)
    col_split = factor(labels_ordered, levels = names(labels_ordered_table))
  } else {
    stop("Column '", group.by, "' cannot be found in metadata!")
  }

  # update limits
  limits <- stats::quantile(as.vector(heatmapdata), probs = c(1-outlier.quantile, outlier.quantile))
  heatmapdata[heatmapdata > limits[2]] <- limits[2]
  heatmapdata[heatmapdata < limits[1]] <- limits[1]

  # highlight some rows
  if(highlight.some){
    ind <- sample(1:nrow(heatmapdata), n_highlight, replace = FALSE)
    ha <- ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = ind, labels = rownames(heatmapdata)[ind], padding = 1,
                                                                        labels_gp = gpar(fontsize = font.size)))
  } else{
    ha <- NULL
  }

  # visualize
  if(is.null(show_row_names))
    show_row_names <- (nrow(heatmapdata) < 30)
  legend_at <- seq(min(heatmapdata), max(heatmapdata), (max(heatmapdata)-min(heatmapdata))/5)
  legend_label <- round(legend_at, 2)
  ComplexHeatmap::Heatmap(heatmapdata,
                          show_row_names = show_row_names, show_row_dend = FALSE, row_names_gp = gpar(fontsize = font.size),
                          show_column_names = FALSE, column_title_rot = 45, column_title_gp = gpar(fontsize = font.size),
                          column_split = col_split, cluster_columns = FALSE, cluster_rows = cluster_rows,
                          show_heatmap_legend = show_heatmap_legend,
                          heatmap_legend_param = list(title = legend_title, at = legend_at, labels = legend_label),
                          right_annotation = ha,
                          col = scales::viridis_pal()(100))
}

####
# Violin Plot ####
####

#' vrViolinPlot
#'
#' @param object A VoltRon object
#' @param features a set of features to be visualized
#' @param assay assay name
#' @param group.by a column from metadata to seperate columns of the heatmap
#' @param norm if TRUE, the normalized data is used
#' @param points if TRUE, measures are visualized as points as well.
#' @param ncol column wise number of plots, for \code{ggarrange}
#' @param nrow row wise number of plots, for \code{ggarrange}
#' @param ... additional parameters passed to \code{getVariableFeatures}
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @export
#'
vrViolinPlot <- function(object, features = NULL, assay = NULL, group.by = "Sample", norm = TRUE, points = TRUE, ncol = 2, nrow = NULL, ...){

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # features
  if(is.null(features))
    stop("You have to define at least one feature")

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get data
  if(any(features %in% vrFeatures(object))){
    selected_features <- features[features %in% vrFeatures(object)]
    violindata <- vrData(object, features = selected_features, assay = assay, norm = norm)
  }

  # get entity type and metadata
  # get entity type and metadata
  # if(is.null(assay.type)){
  #   assay_types <- vrAssayTypes(object, assay = assay)
  #   if(length(unique(assay_types)) == 1){
  #     metadata <- Metadata(object, type = unique(assay_types))
  #   } else {
  #     stop("Please select assay.type as 'cell', 'spot' or 'ROI'")
  #   }
  # } else {
  #   metadata <- Metadata(object, type = assay.type)
  # }
  metadata <- Metadata(object, assay = assay)

  # get feature data
  datax <- lapply(features, function(x){
    if(x %in% vrFeatures(object)){
      return(violindata[x,])
    } else if(x %in% colnames(metadata)){
      return(metadata[,x])
    } else {
      stop("Please provide feature names which reside in either the data or metadata slots!")
    }
  })
  datax <- do.call(cbind, datax)
  colnames(datax) <- features

  # assays
  assays <- stringr::str_extract(rownames(metadata), "Assay[0-9]+$")
  assay_title <- apply(sample.metadata[assays,], 1, function(x) paste(x["Sample"], x["Layer"], x["Assay"], sep = "|"))

  # violin plot
  ggplotdatax <- data.frame(datax,
                      group.by =  as.factor(metadata[[group.by]]),
                      assay_title = assay_title,
                      spatialpoints = rownames(metadata))
  ggplotdatax <- reshape2::melt(ggplotdatax, id.var = c("group.by", "assay_title", "spatialpoints"))
  gg <- ggplot(ggplotdatax, aes(x = group.by, y = value, color = group.by)) +
    geom_violin()

  # visualize points on violin
  if(points){
    gg <- gg +
      geom_point(size = 0.5, position = position_jitter())
  }

  # theme
  gg <- gg +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ylab("") + xlab(group.by) +
    guides(fill = guide_legend(show = FALSE), color = guide_legend(title = group.by, override.aes=list(size = 2)))

  if(length(features) > 1){
    if(length(gg) < ncol) ncol <- length(gg)
    gg <- gg + facet_wrap(.~variable, ncol = ncol, nrow = ceiling(length(gg)/ncol), scales = "free_y")
  } else {
    gg <- gg + labs(title = features)
  }

  return(gg)
}

####
# ROI Plots ####
####

#' vrBarPlot
#'
#' @param object A VoltRon object
#' @param features a set of features to be visualized
#' @param assay assay name
#' @param x.label labels of the x axis
#' @param group.by a column from metadata to seperate columns of the heatmap
#' @param split.by the column to split the barplots by
#' @param norm if TRUE, the normalized data is used
#' @param log if TRUE, data features (excluding metadata features) will be log transformed
#' @param ncol column wise number of plots, for \code{ggarrange}
#' @param nrow row wise number of plots, for \code{ggarrange}
#' @param ... additional parameters passed to \code{getVariableFeatures}
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @export
#'
vrBarPlot <- function(object, features = NULL, assay = NULL, x.label = NULL, group.by = "Sample", split.by = NULL, norm = TRUE, log = FALSE, ncol = 2, nrow = NULL, ...){

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # features
  if(is.null(features))
    stop("You have to define at least one feature")

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get data
  barplotdata <- vrData(object, assay = assay, norm = norm)
  if(log)
    barplotdata <- log(barplotdata)

  # get entity type and metadata
  assay_types <- vrAssayTypes(object, assay = assay)
  if(unique(assay_types) %in% c("spot","cell")){
    stop("vrBarPlot can only be used for ROI assays")
  } else {
    # metadata <- Metadata(object, assay = assay, type = "ROI")
    metadata <- Metadata(object, type = "ROI")
    assy_id <- paste(paste0(assay_names,"$"), collapse = "|")
    metadata <- metadata[grepl(assy_id, rownames(metadata)),]
  }

  # get feature data
  datax <- lapply(features, function(x){
    if(x %in% rownames(barplotdata)){
      return(barplotdata[x,])
    } else if(x %in% colnames(metadata)){
      return(metadata[,x])
    } else{
      stop("Feature '", x, "' cannot be found in data or metadata!")
    }
  })
  datax <- as.data.frame(do.call(cbind, datax))
  colnames(datax) <- features

  # violin plot
  assays <- stringr::str_extract(rownames(metadata), "Assay[0-9]+$")
  assay_title <- apply(sample.metadata[assays,], 1, function(x) paste(x["Sample"], x["Layer"], x["Assay"], sep = "|"))

  # labels and groups
  if(is.null(x.label)) {
    x.labels <- factor(rownames(metadata))
  } else {
    if(x.label %in% colnames(metadata)){
      x.labels <- factor(metadata[[x.label]])
    } else {
      stop("Column '", x.label, "' cannot be found in metadata!")
    }
  }
  if(group.by %in% colnames(metadata)){
    group.by.col <- factor(metadata[[group.by]])
  } else {
    stop("Column '", group.by, "' cannot be found in metadata!")
  }

  # plotting data
  if(is.null(split.by)){
    ggplotdatax <- data.frame(datax,
                              x.label = x.labels,
                              group.by = group.by.col,
                              assay_title = assay_title,
                              spatialpoints = rownames(metadata))
    ggplotdatax <- reshape2::melt(ggplotdatax, id.var = c("x.label", "assay_title", "group.by", "spatialpoints"))
    gg <- ggplot(ggplotdatax, aes(x = x.label, y = value,
                                  fill = factor(group.by, levels = unique(group.by)))) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      ylab("") + xlab(x.label) +
      guides(fill = guide_legend(title = group.by))

    if(length(features) > 1){
      if(length(gg) < ncol) ncol <- length(gg)
      gg <- gg + facet_wrap(.~variable, ncol = ncol, nrow = ceiling(length(features)/ncol), scales = "free")
      return(gg)
    } else {
      gg <- gg + labs(title = features)
      return(gg)
    }
  } else {

    # check split.by column name
    if(split.by %in% colnames(metadata)){
      split.by.col <- factor(metadata[[split.by]])
    } else {
      stop("Column '", split.by, "' cannot be found in metadata!")
    }

    # make ggplot
    ggplotdatax <- data.frame(datax,
                              x.label =  x.labels,
                              group.by = group.by.col,
                              split.by = split.by.col,
                              assay_title = assay_title,
                              spatialpoints = rownames(metadata))
    ggplotdatax <- reshape2::melt(ggplotdatax, id.var = c("x.label", "assay_title", "group.by", "split.by", "spatialpoints"))
    gg <- ggplot(ggplotdatax, aes(x = x.label, y = value,
                                  fill = factor(group.by, levels = unique(group.by)))) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      ylab("") + xlab(x.label) +
      guides(fill = guide_legend(title = group.by))

    if(length(features) > 1){
      # gg <- gg + facet_grid(variable~split.by, scales = "free_x", space = "free")
      gg <- gg + facet_grid(variable~split.by, scales = "free", space = "free_x")
      return(gg)
    } else {
      if(length(gg) < ncol) ncol <- length(gg)
      gg <- gg + facet_wrap(.~split.by, ncol = ncol, nrow = ceiling(length(unique(split.by.col))/ncol), scales = "free_x")
      return(gg)
    }
  }
}

#' vrPercentagePlot
#'
#' @param object A VoltRon object
#' @param assay assay name
#' @param x.label labels of the x axis
#' @param split.by the column to split the barplots by
#' @param split.method either facet_grid or facet_wrap, not affected if \code{split.by} is \code{NULL}
#' @param ncol column wise number of plots, for \code{ggarrange}
#' @param nrow row wise number of plots, for \code{ggarrange}
#' @param ... additional parameters passed to \code{getVariableFeatures}
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @export
#'
vrProportionPlot <- function(object, assay = NULL, x.label = NULL, split.by = NULL, split.method = "facet_wrap", ncol = 2, nrow = NULL, ...){

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get data
  barplotdata <- vrData(object, assay = assay, norm = FALSE)

  # get entity type and metadata
  assay_types <- vrAssayTypes(object, assay = assay)
  if(unique(assay_types) %in% c("spot","cell")){
    stop("vrProportionPlot can only be used for ROI assays")
  } else {
    # metadata <- Metadata(object, assay = assay, type = "ROI")
    metadata <- Metadata(object, type = "ROI")
    assy_id <- paste(paste0(assay_names,"$"), collapse = "|")
    metadata <- metadata[grepl(assy_id, rownames(metadata)),]
  }

  # violin plot
  assays <- stringr::str_extract(rownames(metadata), "Assay[0-9]+$")
  assay_title <- apply(sample.metadata[assays,], 1, function(x) paste(x["Sample"], x["Layer"], x["Assay"], sep = "|"))

  if(is.null(x.label)) {
    x.label <- factor(rownames(metadata))
  } else {
    x.label <- factor(metadata[[x.label]])
  }

  ggplotdatax <- data.frame(t(barplotdata),
                            x.label =  x.label,
                            assay_title = assay_title,
                            spatialpoints = rownames(metadata))
  ggplotdatax <- reshape2::melt(ggplotdatax, id.var = c("x.label", "assay_title", "spatialpoints"))
  ggplotdatax <- ggplotdatax[ggplotdatax$value > 0,]
  gg <- ggplot(ggplotdatax, aes(x = x.label, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ylab("") + xlab("") +
    guides(fill = guide_legend(title = ""))

  if(is.null(split.by)){
    ggplotdatax <- data.frame(t(barplotdata),
                              x.label =  x.label,
                              assay_title = assay_title,
                              spatialpoints = rownames(metadata))
    ggplotdatax <- reshape2::melt(ggplotdatax, id.var = c("x.label", "assay_title", "spatialpoints"))
    ggplotdatax <- ggplotdatax[ggplotdatax$value > 0,]
    gg <- ggplot(ggplotdatax, aes(x = x.label, y = value, fill = variable)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      ylab("") + xlab("") +
      guides(fill = guide_legend(title = ""))
  } else {

    # check split.by column name
    if(split.by %in% colnames(metadata)){
      split.by.col <- factor(metadata[[split.by]])
    } else {
      stop("Column '", split.by, "' cannot be found in metadata!")
    }
    ggplotdatax <- data.frame(t(barplotdata),
                              x.label =  x.label,
                              assay_title = assay_title,
                              split.by = split.by.col,
                              spatialpoints = rownames(metadata))
    ggplotdatax <- reshape2::melt(ggplotdatax, id.var = c("x.label", "assay_title", "split.by", "spatialpoints"))
    ggplotdatax <- ggplotdatax[ggplotdatax$value > 0,]
    gg <- ggplot(ggplotdatax, aes(x = x.label, y = value, fill = variable)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      ylab("") + xlab("") +
      guides(fill = guide_legend(title = ""))

    # split
    if(length(gg) < ncol) ncol <- length(gg)
    if(split.method == "facet_wrap"){
      gg <- gg + facet_wrap(.~split.by, ncol = ncol, nrow = ceiling(length(unique(split.by.col))/ncol), scales = "free_x")
    } else {
      gg <- gg + facet_grid(.~split.by, scales = "free_x", space = "free_x", )
    }
  }
  return(gg)
}
