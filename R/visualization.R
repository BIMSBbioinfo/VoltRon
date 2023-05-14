####
# Feature plots ####
####

#' SpatFeatPlot
#'
#' Functions for plotting spatial objects
#'
#' @param object SpaceRover object
#' @param features a vector of features to visualize
#' @param sample sample names
#' @param group.by a grouping label for the spatial entities
#' @param assay the assay type name
#' @param ncol ncol
#' @param nrow nrow
#' @param font.size font sizes
#' @param pt.size point size
#' @param alpha alpha
#' @param keep.scale whether unify all scales for all features or not
#' @param label
#'
#' @export
#'
SpatFeatPlot <- function(object, features, group.by = "label", assay = "GeoMx", ncol = 2, nrow = NULL, font.size = 2, pt.size = 10, alpha = 1, keep.scale = "feature", label = FALSE) {

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # list of plots
  gg <- list()

  # get assay names
  if(assay %in% sample.metadata$Assay){
    assay_names <- rownames(sample.metadata)[sample.metadata$Assay %in% assay]
  } else {
    if(assay %in% rownames(sample.metadata)) {
      assay_names <- assay
    } else {
      stop("Assay name or type is not found in the object")
    }
  }


  # calculate limits for plotting, all for making one scale, feature for making multiple
  limits <- Map(function(feat){
    range_feat <- Map(function(assy){
      normdata <- object[[assy]]@normdata
      metadata <- Metadata(object, type = "ROI")
      metadata <- metadata[grepl(assy, rownames(metadata)),]
      if(feat %in% rownames(normdata)){
        range(normdata[feat, ])
      } else {
        range(metadata[,feat])
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
  # if(length(features) > 1 && length(assay_names) > 1){
  #   plot_title <- assay_title
  #   legend_title <- feature_title
  # } else if(length(features) > 1 && length(assay_names) == 1){
  #   plot_title <- feature_title
  #   legend_title <- as.list(rep("Log.Exp", length(features)))
  #   names(legend_title) <- features
  # } else if(length(features) == 1 && length(assay_names) > 1){
  #   plot_title <- assay_title
  #   legend_title <- feature_title
  # } else {
  #   plot_title <- as.list(assay_title)
  #   names(plot_title) <- assay_names
  #   legend_title <- as.list(features)
  #   names(legend_title) <- features
  # }
  plot_title <- assay_title
  legend_title <- feature_title

  # for each feature
  i <- 1
  for(feat in features){

    # for each assay
    for(assy in assay_names){

      # get assay
      cur_assay <- object[[assy]]
      metadata <- Metadata(object, type = "ROI")
      metadata <- metadata[grepl(assy, rownames(metadata)),]

      # visualize
      p_title <- plot_title[[assy]]
      l_title <- legend_title[[feat]]
      gg[[i]] <- SpatFeatPlotSingle(assay = cur_assay, metadata = metadata, feature = feat, limits = limits[[feat]][[assy]],
                              group.by = group.by, font.size = font.size, pt.size = pt.size, alpha = alpha,
                              label = label, plot_title = p_title, legend_title = l_title)
      i <- i + 1
    }
  }

  # return a list of plots or a single one
  if(length(features) > 1 && length(assay_names) > 1){
    return(ggarrange(plotlist = gg, ncol = length(features), nrow = length(assay_names)))
  } else if(length(features) > 1 && length(assay_names) == 1){
    return(ggarrange(plotlist = gg, ncol = ncol, nrow = ceiling(length(gg)/ncol)))
  } else if(length(features) == 1 && length(assay_names) > 1){
    return(ggarrange(plotlist = gg, ncol = ncol, nrow = ceiling(length(gg)/ncol)))
  } else {
    return(gg[[1]])
  }
}

#' SpatFeatPlotSingle
#'
#' A single Spatial Feature plot of spacerover objects
#'
#' @param assay
#' @param metadata
#' @param feature
#' @param limits
#' @param group.by
#' @param font.size
#' @param pt.size
#' @param alpha
#' @param label
#' @param plot_title
#' @param legend_title
#'
SpatFeatPlotSingle <- function(assay, metadata, feature, limits, group.by = "label",
                               font.size = 2, pt.size = 10, alpha = 1, label = FALSE, plot_title = NULL, legend_title = NULL){

  # data
  info <- image_info(assay@image)
  image <- assay@image
  coords <- as.data.frame(assay@coords)
  normdata <- assay@normdata

  # plotting features
  coords[[group.by]] <- metadata[,group.by]

  # get data
  if(feature %in% rownames(normdata)){
    coords$score <- normdata[feature,]
  } else {
    coords$score <- metadata[,feature]
  }

  # get image information and plotting features
  info <- image_info(image)
  midpoint <- sum(limits)/2

  # add image
  g <- ggplot() +
    ggplot2::annotation_raster(image, 0, info$width, info$height, 0, interpolate = FALSE)

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
  } else {
    g <- g +
      geom_point(mapping = aes(x = x, y = y, fill = score), coords, shape = 21, size = pt.size, alpha = alpha)
  }

  # adjust gradient
  g <- g +
    scale_fill_gradientn(name = legend_title,
                         colors=c("dodgerblue2", "white", "yellow3"),
                         values=scales::rescale(c(limits[1], midpoint, limits[2])), limits = limits)

  # more visualization parameters
  g <- g +
    ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,-10,0)), panel.background = element_blank(), panel.grid.minor = element_blank(),
                                axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                                legend.margin = margin(0,0,0,-10)) +
    xlim(0,info$width) + ylim(0, info$height)

  # visualize labels
  if(label){
    g <- g + geom_label_repel(mapping = aes_string(x = "x", y = "y", label = group.by), coords,
                                box.padding = 0.5, size = font.size, direction = "both", seed = 1)
  }

  # return data
  return(g)
}

####
# FeatPlot ####
####

FeatPlot <- function(object, features, group.by = "label", assay = "GeoMx") {

}
