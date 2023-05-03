####
# ROI plotting ####
####

#' SFPlot
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
#'
#' @export
#'
SFPlot <- function(object, features, ncol = 3, nrow = NULL, sample = NULL, font.size = 2, keep.scale = "feature", mode = "ROI") {

  if(nrow(object@sample.metadata) > 1)
    stop("Plotting can only be performed on these SpaceRover objects with a single assay")

  # get assay
  assay <- GeoMxR1_subset[[object@sample.metadata$Sample, object@sample.metadata$Layer]][[MainAssay(object)]]

  # normalize
  info <- image_info(assay@image)
  coords <- as.data.frame(assay@coords)
  normdata <- assay@rawdata
  sizefactor <- colSums(normdata)
  sizefactor <- matrix(rep(sizefactor, nrow(normdata)), byrow = T, nrow = nrow(normdata))
  normdata <- (normdata/sizefactor)*1000000
  normdata <- log(normdata + 1)

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
    limits <- lapply(1:length(features), range)
  }
  midpoint <- lapply(limits, function(x) sum(x)/2)

  # visualize
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
                             values=scales::rescale(c(limits[[i]][1],midpoint[[i]],limits[[i]][2])), limits = limits[[i]]) + # scale_fill_continuous(name = "Log.Exp.", type = "viridis") +
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
    coords$feature <- normdata[features,]
    g1 <- ggplot() +
      ggplot2::coord_fixed(expand = FALSE, xlim = c(0, info$width), ylim = c(0, info$height)) +
      ggplot2::annotation_raster(assay@image, 0, info$width, info$height, 0, interpolate = FALSE) +
      geom_point(mapping = aes(x = x, y = y, fill = feature), coords, shape = 21, size = 10) +
      scale_fill_gradientn(name = "Log.Exp.",
                           # low = "darkgreen", mid = "white", high = "yellow2",
                           colors=c("dodgerblue2", "white", "yellow3"),
                           values=scales::rescale(c(limits[[1]][1],midpoint[[i]],limits[[1]][2])), limits = limits[[1]]) + # scale_fill_continuous(name = "Log.Exp.", type = "viridis") +
      NoAxes() +
      ggtitle(paste(features)) +
      theme(plot.title = element_text(hjust = 0.5))
    g1 <- g1 + geom_label_repel(mapping = aes(x = x, y = y, label = sample_names), coords,
                                box.padding = 0.5, size = font.size, direction = "y", seed = 1)
    return(g1)
  }
}
