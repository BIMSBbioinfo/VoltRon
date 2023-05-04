####
# ROI plotting ####
####

#' SFPlotMulti
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
#' @param keep.scale whether unify all scales for all features or not
#'
#' @export
#'
SFPlot <- function(object, features, group.by = "label", assay = "GeoMx", ncol = 2, nrow = NULL, font.size = 2, keep.scale = "feature", label = FALSE) {

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
      range(object[[assy]]@normdata[feat, ])
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
  if(length(features) > 1 && length(assay_names) > 1){
    plot_title <- assay_title
    legend_title <- feature_title
  } else if(length(features) > 1 && length(assay_names) == 1){
    plot_title <- feature_title
    legend_title <- as.list(rep("Log.Exp", length(features)))
    names(legend_title) <- features
  } else if(length(features) == 1 && length(assay_names) > 1){
    plot_title <- assay_title
    legend_title <- feature_title
  } else {
    plot_title <- as.list(assay_title)
    names(plot_title) <- assay_names
    legend_title <- as.list(features)
    names(legend_title) <- features
  }

  # for each feature
  i <- 1
  for(feat in features){

    # for each assay
    for(assy in assay_names){

      # get assay
      cur_assay <- object[[assy]]

      # normalize
      info <- image_info(cur_assay@image)
      image <- cur_assay@image
      coords <- as.data.frame(cur_assay@coords)
      normdata <- cur_assay@normdata

      # get data
      coords$score <- normdata[feat,]

      # plotting features
      group.by_labels <- Metadata(object, type = "ROI")
      coords[[group.by]] <- group.by_labels[grepl(assy, rownames(group.by_labels)),group.by]
      p_title <- plot_title[[assy]]
      l_title <- legend_title[[feat]]

      # visualize
      gg[[i]] <- SFPlotSingle(coords = coords, image = image, feature = feat, limits = limits[[feat]][[assy]],
                              group.by = group.by, font.size = font.size, keep.scale = keep.scale,
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

SFPlotSingle <- function(coords, image, feature, limits, group.by = "label", font.size = 2, keep.scale = "feature", label = FALSE, plot_title = NULL, legend_title = NULL){

  # get image information and plotting features
  info <- image_info(image)
  midpoint <- sum(limits)/2

  # visualize with ggplot
  # g <- ggplot() +
  #   ggplot2::coord_fixed(expand = FALSE, xlim = c(0, info$width), ylim = c(0, info$height)) +
  #   ggplot2::annotation_raster(image, 0, info$width, info$height, 0, interpolate = TRUE) +
  #   geom_point(mapping = aes(x = x, y = y, fill = score), coords, shape = 21, size = 10) +
  #   scale_fill_gradientn(name = legend_title,
  #                        colors=c("dodgerblue2", "white", "yellow3"),
  #                        values=scales::rescale(c(limits[1], midpoint, limits[2])), limits = limits) +
  #   NoAxes() +
  #   ggtitle(plot_title) +
  #   theme(plot.title = element_text(hjust = 0.5))

  # visualize with ggplot
  g <- ggplot() +
    ggplot2::annotation_raster(image, 0, info$width, info$height, 0, interpolate = FALSE) +
    geom_point(mapping = aes(x = x, y = y, fill = score), coords, shape = 21, size = 10) +
    scale_fill_gradientn(name = legend_title,
                         colors=c("dodgerblue2", "white", "yellow3"),
                         values=scales::rescale(c(limits[1], midpoint, limits[2])), limits = limits) +
    ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5, margin=margin(0,0,-15,0)), panel.background = element_blank(), panel.grid.minor = element_blank(),
                                axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),
                                axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                                legend.margin = margin(0,0,0,-20)) +
    xlim(0,info$width) + ylim(0, info$height)

  # visualize labels
  if(label){
    g <- g + geom_label_repel(mapping = aes_string(x = "x", y = "y", label = group.by), coords,
                                box.padding = 0.5, size = font.size, direction = "y", seed = 1)
  }

  # return data
  return(g)
}

####
# old functions ####
####

#' SFPlotold
#'
#' Functions for plotting spatial objects
#'
#' @param object SpaceRover object
#' @param features a vector of features to visualize
#' @param ncol ncol
#' @param nrow nrow
#' @param sample sample names matching spatial points (i.e. entities)
#' @param font.size font sizes
#' @param keep.scale whether unify all scales for all features or not
#' @param mode plotting mode for different spatial entities, ROI, Spot or cell
#' @param split.by split by a metadata element
#'
#' @export
#'
SFPlotold <- function(object, features, ncol = 3, nrow = NULL, sample = NULL, font.size = 2, keep.scale = "feature", mode = "ROI") {

  # get assay
  assay <- GeoMxR1_subset[[object@sample.metadata$Sample, object@sample.metadata$Layer]][[MainAssay(object)]]

  # normalize
  info <- image_info(assay@image)
  image <- assay@image
  coords <- as.data.frame(assay@coords)
  normdata <- assay@normdata

  # Check keep.scale param for valid entries
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature", "all"))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }

  # calculate limits for plotting
  if(keep.scale == "all") {
    normdata <- normdata[features,, drop = FALSE]
    limits <- lapply(1:length(features), function(x) range(normdata))
  } else if(keep.scale == "feature"){
    normdata <- normdata[features,, drop = FALSE]
    limits <- lapply(features, function(feat) range(normdata[feat,]))
  }
  midpoint <- lapply(limits, function(x) sum(x)/2)

  # visualize for either one feature or multiple features
  if(length(features) > 1){
    gg <- list()
    for(i in 1:length(features)){
      coords$feature <- normdata[features[i],]
      g1 <- ggplot() +
        ggplot2::coord_fixed(expand = FALSE, xlim = c(0, info$width), ylim = c(0, info$height)) +
        ggplot2::annotation_raster(assay@image, 0, info$width, info$height, 0, interpolate = FALSE) +
        geom_point(mapping = aes(x = x, y = y, fill = feature), coords, shape = 21, size = 10) +
        scale_fill_gradientn(name = "Log.Exp.",
                             # low = "darkgreen", mid = "white", high = "yellow2",
                             colors=c("dodgerblue2", "white", "yellow3"),
                             values=scales::rescale(c(limits[[i]][1],midpoint[[i]],limits[[i]][2])), limits = limits[[i]]) +
        NoAxes() +
        ggtitle(paste(features[i])) +
        theme(plot.title = element_text(hjust = 0.5))
      gg[[i]]  <- g1 + geom_label_repel(mapping = aes(x = x, y = y, label = sample), coords,
                                        box.padding = 0.5, size = font.size, direction = "y", seed = 1)
    }
    if(length(gg) > ncol){
      nrow <- ceiling(length(gg)/ncol)
    } else{
      ncol <- length(gg)
    }
    return(ggarrange(plotlist = gg, ncol = ncol, nrow = nrow))
  } else {
    coords$score <- normdata[features,]
    g1 <- ggplot() +
      ggplot2::coord_fixed(expand = FALSE, xlim = c(0, info$width), ylim = c(0, info$height)) +
      ggplot2::annotation_raster(assay@image, 0, info$width, info$height, 0, interpolate = FALSE) +
      geom_point(mapping = aes(x = x, y = y, fill = score), coords, shape = 21, size = 10) +
      scale_fill_gradientn(name = "Log.Exp.",
                           # low = "darkgreen", mid = "white", high = "yellow2",
                           colors=c("dodgerblue2", "white", "yellow3"),
                           values=scales::rescale(c(limits[[1]][1],midpoint[[1]],limits[[1]][2])), limits = limits[[1]]) +
      NoAxes() +
      ggtitle(paste(features[1])) +
      theme(plot.title = element_text(hjust = 0.5))
    g1 <- g1 + geom_label_repel(mapping = aes(x = x, y = y, label = sample_names), coords,
                                box.padding = 0.5, size = font.size, direction = "y", seed = 1)
    return(g1)
  }
}

#' SFPlotMulti2
#'
#' Functions for plotting spatial objects
#'
#' @param object SpaceRover object
#' @param features a vector of features to visualize
#' @param sample sample names matching spatial points (i.e. entities)
#' @param split.by split by a metadata element
#' @param mode plotting mode for different spatial entities, ROI, Spot or cell
#' @param ncol ncol
#' @param nrow nrow
#' @param font.size font sizes
#' @param keep.scale whether unify all scales for all features or not
#'
#' @export
#'
SFPlotMulti2 <- function(object, features, sample = NULL, split.by = NULL, assay = "GeoMx", ncol = 3, nrow = NULL, font.size = 2, keep.scale = "feature") {

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # list of plots
  gg <- list()

  # normalize
  info <- image_info(assay@image)
  image <- assay@image
  coords <- as.data.frame(assay@coords)
  normdata <- assay@normdata

  # Check keep.scale param for valid entries
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature", "all"))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }

  # calculate limits for plotting, all for making one scale, feature for making multiple
  if(keep.scale == "all") {
    normdata <- normdata[features,, drop = FALSE]
    limits <- lapply(1:length(features), function(x) range(normdata))
  } else if(keep.scale == "feature"){
    normdata <- normdata[features,, drop = FALSE]
    limits <- lapply(features, function(feat) range(normdata[feat,]))
  }

  # visualize for either one feature or multiple features
  if(length(features) > 1){
    gg <- list()
    for(i in 1:length(features)){
      coords$score <- normdata[features[i],]
      gg[[i]] <- SFPlotSingle(coords, image, limits[[i]], sample, font.size, features[i])
    }
    if(length(gg) > ncol){
      nrow <- ceiling(length(gg)/ncol)
    } else{
      ncol <- length(gg)
    }
    return(ggarrange(plotlist = gg, ncol = ncol, nrow = nrow))
  } else {
    coords$score <- normdata[features,]
    g1 <- SFPlotSingle(coords, image, limits[[1]], sample, font.size, features[1])
    return(g1)
  }
}
