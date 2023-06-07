####
# Spatial plots ####
####

####
## Spatial Identity Plot ####
####

#' SpatPlot
#'
#' Plotting identification of spatially resolved cells, spots, and ROI on associated images from multiple assays in a SpaceRover object.
#'
#' @param object SpaceRover object
#' @param group.by a grouping label for the spatial entities
#' @param assay the assay name
#' @param assay.type the assay type name: 'cell', 'spot' or 'ROI'
#' @param ncol column wise number of plots, for \code{ggarrange}
#' @param nrow row wise number of plots, for \code{ggarrange}
#' @param font.size font sizes
#' @param pt.size point size
#' @param alpha alpha level for cells/spots/ROIs
#' @param label if TRUE, the labels of the ROI assays will be visualized
#' @param background the background of the plot, either "image" for overlaying the image of the assay, or "black" or "white" background (suitable for IF based assays)
#' @param crop whether to crop an image of a spot assay
#' @param common.legend whether to use a common legend for all plots
#' @param collapse whether to combine all ggplots
#'
#' @importFrom ggpubr ggarrange
#' @export
#'
SpatPlot <- function(object, group.by = "label", assay = NULL, assay.type = NULL, ncol = 2, nrow = NULL,
                     font.size = 2, pt.size = 2, alpha = 0.6, label = FALSE, background = "image",
                     crop = FALSE, common.legend = TRUE, collapse = TRUE) {

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # list of plots
  gg <- list()

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get entity type and metadata
  if(is.null(assay.type)){
    # assay_types <- unlist(lapply(assay_names, function(x) object[[x]]@type))
    assay_types <- AssayTypes(object, assay = assay)
    if(length(unique(assay_types)) == 1){
      metadata <- Metadata(object, type = unique(assay_types))
    } else {
      stop("Please select assay.type as 'cell', 'spot' or 'ROI'")
    }
  } else {
    metadata <- Metadata(object, type = assay.type)
  }

  # configure titles
  plot_title <- as.list(apply(sample.metadata[assay_names,], 1, function(x) paste(x["Sample"], x["Layer"], sep = ", ")))
  names(plot_title) <- assay_names

  # for each assay
  i <- 1
  for(assy in assay_names){

    # get assay
    cur_assay <- object[[assy]]
    assy_id <- paste0(assy,"$")
    cur_metadata <- metadata[grepl(assy_id, rownames(metadata)),]

    # visualize
    p_title <- plot_title[[assy]]
    gg[[i]] <- SpatPlotSingle(assay = cur_assay, metadata = cur_metadata,
                              group.by = group.by, font.size = font.size, pt.size = pt.size, alpha = alpha,
                              plot_title = p_title, background = background, crop = crop)
    i <- i + 1
  }

  # return a list of plots or a single one
  if(collapse){
    if(length(assay_names) > 1){
      if(length(gg) < ncol) ncol <- length(gg)
      return(ggarrange(plotlist = gg, ncol = ncol, nrow = ceiling(length(gg)/ncol), common.legend = common.legend, legend = "right"))
    } else {
      return(gg[[1]])
    }
  } else {
    gg
  }
}

#' SpatPlotSingle
#'
#' Plotting a single assay from a SpaceRover object. We plot the identification of spatially resolved cells, spots, and ROI on associated images.
#'
#' @param assay srAssay object
#' @param metadata the metadata associated with the assay
#' @param group.by a grouping label for the spatial entities
#' @param font.size font sizes
#' @param pt.size point size
#' @param alpha alpha level for cells/spots/ROIs
#' @param plot_title the title of the single plot
#' @param background the background of the plot, either "image" for overlaying the image of the assay, or "black" or "white" background (suitable for IF based assays)
#' @param crop whether to crop an image of a spot assay
#'
#' @import ggplot2
#'
SpatPlotSingle <- function(assay, metadata, group.by = "label", font.size = 2, pt.size = 2, alpha = 0.6, plot_title = NULL, background = "image", crop = FALSE){

  # data
  info <- image_info(assay@image)
  image <- assay@image
  coords <- as.data.frame(assay@coords)
  normdata <- assay@normdata

  # plotting features
  coords[[group.by]] <- as.factor(metadata[,group.by])

  # get image information and plotting features
  info <- image_info(image)

  # plot
  g <- ggplot()

  # add image
  if(background == "image") {
    g <- g +
      ggplot2::annotation_raster(image, 0, info$width, info$height, 0, interpolate = FALSE)
  }

  # add points or segments
  if(assay@type == "spot"){
    g <- g +
      coord_fixed(xlim = c(0,info$width), ylim = c(0,info$height)) +
      geom_spot(mapping = aes_string(x = "x", y = "y", fill = group.by), coords, shape = 21, alpha = alpha, spot.radius = assay@params[["spot.radius"]]) +
      guides(fill = guide_legend(override.aes=list(shape = 21, size = 4, lwd = 0.1)))
  } else if(assay@type == "cell") {
    g <- g +
      geom_point(mapping = aes_string(x = "x", y = "y", fill = group.by, color = group.by), coords, shape = 21, size = rel(pt.size), alpha = alpha)
  } else {
    stop("Only spots and cells can be visualized with SpatPlot!")
  }

  # more visualization parameters
  g <- g +
    ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,0,0)),
                                panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                                axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                                legend.margin = margin(0,0,0,0))

  # set up the limits
  if(assay@type == "spot"){
    if(crop){
      g <- g +
        coord_fixed(xlim = range(coords$x), ylim = range(coords$y))
    } else {
      g <- g +
        coord_fixed(xlim = c(0,info$width), ylim = c(0,info$height))
    }
  } else {
    g <- g +
      xlim(0,info$width) + ylim(0, info$height)
  }

  # background
  if(any(background %in% c("white","black"))){
    g <- g +
      theme(panel.background = element_rect(fill = background, colour = background, size = 0.5, linetype = "solid"))
  } else{
    g <- g +
      theme(panel.background = element_blank())
  }

  # return data
  return(g)
}

####
## Spatial Feature Plot ####
####

#' SpatFeatPlot
#'
#' Plotting single/multiple features of spatially resolved cells, spots, and ROI on associated images from multiple assays in a SpaceRover object.
#'
#' @param object SpaceRover object
#' @param features a set of features, either from the rows of rawdata, normdata or columns of the metadata
#' @param group.by a grouping label for the spatial entities
#' @param assay the assay name
#' @param assay.type the assay type name: 'cell', 'spot' or 'ROI'
#' @param ncol column wise number of plots, for \code{ggarrange}
#' @param nrow row wise number of plots, for \code{ggarrange}
#' @param font.size font sizes
#' @param pt.size point size
#' @param alpha alpha level for cells/spots/ROIs
#' @param keep.scale whether unify all scales for all features or not
#' @param label if TRUE, labels of ROIs will be visualized too
#' @param background the background of the plot, either "image" for overlaying the image of the assay, or "black" or "white" background (suitable for IF based assays)
#' @param crop whether to crop an image of a spot assay
#' @param common.legend whether to use a common legend for all plots
#' @param collapse whether to combine all ggplots
#'
#' @importFrom ggpubr ggarrange
#' @export
#'
SpatFeatPlot <- function(object, features, group.by = "label", assay = NULL, assay.type = NULL, ncol = 2, nrow = NULL,
                         font.size = 2, pt.size = 2, alpha = 0.6, keep.scale = "feature", label = FALSE, background = "image",
                         crop = FALSE, common.legend = TRUE, collapse = TRUE) {

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # list of plots
  gg <- list()

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get entity type and metadata
  if(is.null(assay.type)){
    assay_types <- AssayTypes(object, assay = assay)
    if(length(unique(assay_types)) == 1){
      assay.type <- unique(assay_types)
      metadata <- Metadata(object, type = assay.type)
    } else {
      stop("Please select assay.type as 'cell', 'spot' or 'ROI'")
    }
  } else {
    metadata <- Metadata(object, type = assay.type)
  }

  # calculate limits for plotting, all for making one scale, feature for making multiple
  limits <- Map(function(feat){
    range_feat <- Map(function(assy){
      normdata <- Data(object[[assy]], norm = TRUE)
      metadata <- Metadata(object, type = assay.type)
      metadata <- metadata[grepl(assy, rownames(metadata)),]
      if(feat %in% rownames(normdata)){
        range(normdata[feat, ])
      } else {
        if(feat %in% colnames(metadata)){
          range(metadata[,feat])
        } else {
          stop("Feature ", feat, " cannot be found in data or metadata!")
        }
      }
    }, assay_names)
    if(keep.scale == "all"){
      range_feat_all <- c(do.call(min, range_feat), do.call(max, range_feat))
      range_feat <- Map(function(assy) return(range_feat_all), assay_names)
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
  for(feat in features){

    # for each assay
    for(assy in assay_names){

      # get assay
      cur_assay <- object[[assy]]
      cur_metadata <- metadata[grepl(assy, rownames(metadata)),]

      # visualize
      p_title <- plot_title[[assy]]
      l_title <- legend_title[[feat]]
      gg[[i]] <- SpatFeatPlotSingle(assay = cur_assay, metadata = cur_metadata, feature = feat, limits = limits[[feat]][[assy]],
                              group.by = group.by, font.size = font.size, pt.size = pt.size, alpha = alpha,
                              label = label, plot_title = p_title, legend_title = l_title, background = background, crop = crop)
      i <- i + 1
    }
  }

  if(collapse){
    # return a list of plots or a single one
    if(length(features) > 1 && length(assay_names) > 1){
      return(ggarrange(plotlist = gg, ncol = length(features), nrow = length(assay_names)))
    } else if(length(features) > 1 && length(assay_names) == 1){
      if(length(gg) < ncol) ncol <- length(gg)
      return(ggarrange(plotlist = gg, ncol = ncol, nrow = ceiling(length(gg)/ncol)))
    } else if(length(features) == 1 && length(assay_names) > 1){
      if(length(gg) < ncol) ncol <- length(gg)
      return(ggarrange(plotlist = gg, ncol = ncol, nrow = ceiling(length(gg)/ncol), common.legend = common.legend, legend = "right"))
    } else {
      return(gg[[1]])
    }
  } else {
    return(gg)
  }
}

#' SpatFeatPlotSingle
#'
#' A single Spatial Feature plot of spacerover objects
#'
#' @param assay srAssay object
#' @param metadata the metadata associated with the assay
#' @param feature a feature, either from the rows of rawdata, normdata or columns of the metadata
#' @param limits limits of the legend of the plot
#' @param group.by a grouping label for the spatial entities
#' @param font.size font sizes
#' @param pt.size point size
#' @param alpha alpha level for cells/spots/ROIs
#' @param label if TRUE, labels of ROIs will be visualized too
#' @param plot_title the main title of the single plot
#' @param legend_title the legend title of the single plot
#' @param background the background of the plot, either "image" for overlaying the image of the assay, or "black" or "white" background (suitable for IF based assays)
#' @param crop whether to crop an image of a spot assay
#'
#' @import ggplot2
#'
SpatFeatPlotSingle <- function(assay, metadata, feature, limits, group.by = "label",
                               font.size = 2, pt.size = 2, alpha = 0.6, label = FALSE, plot_title = NULL, legend_title = NULL, background = "image", crop = FALSE){

  # data
  info <- image_info(assay@image)
  image <- assay@image
  coords <- as.data.frame(assay@coords)
  normdata <- assay@normdata

  # get data
  if(feature %in% rownames(normdata)){
    coords$score <- normdata[feature,]
  } else {
    coords$score <- metadata[,feature]
  }

  # get image information and plotting features
  info <- image_info(image)
  midpoint <- sum(limits)/2

  # plot
  g <- ggplot()

  # add image
  if(background == "image") {
    g <- g +
      ggplot2::annotation_raster(image, 0, info$width, info$height, 0, interpolate = FALSE)
  }

  # add points or segments
  if(assay@type == "ROI" && !is.null(assay@segments)){
    polygon_data <- NULL
    circle_data <- NULL
    for(i in 1:length(assay@segments)){
      cur_data <- as.data.frame(cbind(assay@segments[[i]], names(assay@segments)[i], coords$score[i]))
      if(nrow(assay@segments[[i]]) > 1){
        colnames(cur_data) <- c("x", "y", "segment", "score")
        polygon_data <- as.data.frame(rbind(polygon_data, cur_data))
      } else {
        colnames(cur_data) <- c("x", "y", "rx", "ry", "segment", "score")
        circle_data <- as.data.frame(rbind(circle_data,  cur_data))
      }
    }
    if(!is.null(geom_polygon)){
      g <- g +
        geom_polygon(aes(x = x, y = y, fill = score, group = segment), data = polygon_data, alpha = alpha)
    }
    if(!is.null(circle_data)){
      g <- g +
        geom_ellipse(aes(x0 = as.numeric(x), y0 = as.numeric(y), a = as.numeric(rx), b = as.numeric(ry), angle = 0,
                         fill = score, group = segment), data = circle_data, lwd = 0, alpha = alpha)
    }
    g <- g +
      scale_fill_gradientn(name = legend_title,
                           colors=c("dodgerblue2", "white", "yellow3"),
                           values=scales::rescale(c(limits[1], midpoint, limits[2])), limits = limits)
  } else if(assay@type == "spot"){
    g <- g +
      geom_spot(mapping = aes(x = x, y = y, fill = score), coords, shape = 21, alpha = alpha, spot.radius = assay@params[["spot.radius"]]) +
      scale_fill_gradientn(name = legend_title,
                             colors=c("dodgerblue3", "yellow", "red"),
                             values=scales::rescale(c(limits[1], midpoint, limits[2])), limits = limits)
  } else {
    g <- g +
      geom_point(mapping = aes(x = x, y = y, fill = score), coords, shape = 21, size = rel(pt.size), alpha = alpha) +
      scale_fill_gradientn(name = legend_title,
                             colors=c("dodgerblue2", "white", "yellow3"),
                             values=scales::rescale(c(limits[1], midpoint, limits[2])), limits = limits)
  }

  # more visualization parameters
  g <- g +
    ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,0,0)),
                                panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                                axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                                legend.margin = margin(0,0,0,0))

  # set up the limits
  if(assay@type == "spot"){
    if(crop){
      g <- g +
        coord_fixed(xlim = range(coords$x), ylim = range(coords$y))
    } else {
      g <- g +
        coord_fixed(xlim = c(0,info$width), ylim = c(0,info$height))
    }
  } else {
    g <- g +
      xlim(0,info$width) + ylim(0, info$height)
  }

  # background
  if(any(background %in% c("white","black"))){
    g <- g +
      theme(panel.background = element_rect(fill = background, colour = background, size = 0.5, linetype = "solid"))
  } else{
    g <- g +
      theme(panel.background = element_blank())
  }

  # visualize labels
  if(label){
    coords[[group.by]] <- metadata[,group.by]
    g <- g + geom_label_repel(mapping = aes_string(x = "x", y = "y", label = group.by), coords,
                                box.padding = 0.5, size = font.size, direction = "both", seed = 1)
  }

  # return data
  return(g)
}

####
## Spatial Auxiiary ####
####

geom_spot <- function (mapping = NULL, data = NULL, stat = "identity", position = "identity",
                           ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSpot,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = rlang::list2(na.rm = na.rm, ...))
}

GeomSpot <- ggproto("GeomSpot",
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
                                 size = unit(spot.radius, "npc"),
                                 gp = grid::gpar(
                                   col = alpha(coords$colour, coords$alpha),
                                   fill = alpha(coords$fill, coords$alpha),
                                   fontsize = coords$size * .pt + stroke_size * .stroke / 2,
                                   lwd = coords$stroke * .stroke / 2
                                 )
                               )
                        )
                      },

                      draw_key = draw_key_point
)

####
# Scatter Plot ####
####

ScatterFeaturePlot <- function(object, feature.1, feature.2, norm = TRUE, assay = NULL, assay.type = NULL,
                               pt.size = 2, font.size = 2, group.by = "label", label = FALSE, trend = FALSE){

  # check the number of features
  if(is.null(feature.1) | is.null(feature.2))
    stop("Please provide both 'feature.1' and 'feature.2'")

  # check the number of features
  if((length(feature.1) != 1 | length(feature.2) != 1))
    stop("Both 'feature.1' and 'feature.2' should be of length 1.")

  # data
  normdata <- Data(object, assay = assay, norm = norm)

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get entity type and metadata
  if(is.null(assay.type)){
    assay_types <- AssayTypes(object, assay = assay)
    if(length(unique(assay_types)) == 1){
      assay.type <- unique(assay_types)
      metadata <- Metadata(object, type = assay.type)
    } else {
      stop("Please select assay.type as 'cell', 'spot' or 'ROI'")
    }
  } else {
    metadata <- Metadata(object, type = assay.type)
  }

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
    g <- g + geom_label_repel(mapping = aes_string(x = feature.1, y = feature.2, label = group.by), data_feature,
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
#' @param object spacerover object
#' @param assay
#' @param assay.type
#' @param group.by
#' @param norm
#'
#' @import ComplexHeatmap
#'
HeatmapPlot <- function(object, assay = NULL, assay.type = NULL, group.by = "clusters", norm = TRUE){

  # data
  heatmapdata <- Data(object, assay = assay, norm = norm)

  # get entity type and metadata
  if(is.null(assay.type)){
    assay_types <- AssayTypes(object, assay = assay)
    if(length(unique(assay_types)) == 1){
      assay.type <- unique(assay_types)
      metadata <- Metadata(object, type = assay.type)
    } else {
      stop("Please select assay.type as 'cell', 'spot' or 'ROI'")
    }
  } else {
    metadata <- Metadata(object, type = assay.type)
  }
  metadata <- metadata[colnames(heatmapdata),]

  # manage data for plotting
  heatmapdata_scale <- apply(heatmapdata, 1, scale)
  heatmapdata_scale <- t(heatmapdata_scale)
  heatmapdata_scale_order <- heatmapdata_scale[,order(metadata[[group.by]], decreasing = FALSE)]
  labels_ordered <- metadata[[group.by]][order(metadata[[group.by]], decreasing = FALSE)]
  labels_ordered_table <- table(labels_ordered)
  col_split = factor(labels_ordered, levels = names(labels_ordered_table))
  ComplexHeatmap::Heatmap(heatmapdata_scale_order, show_row_names = TRUE, show_column_names = FALSE,
          column_split = col_split, cluster_columns = FALSE)
}
