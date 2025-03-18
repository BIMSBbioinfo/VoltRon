#' @include zzz.R
#'
NULL

####
# General ####
####

#' vrFeatures
#'
#' Get names of the features.
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#' 
#' @rdname vrFeatures
#' @export vrFeatures
#' @order 1
setGeneric("vrFeatures", function(object, ...) standardGeneric("vrFeatures"))

#' Feature Data
#'
#' Get and set feature data from the main.assay
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrFeatureData
#' @export vrFeatureData
#' @order 1
setGeneric("vrFeatureData", function(object, ...) standardGeneric("vrFeatureData"))

#' @rdname vrFeatureData
#' @export vrFeatureData<-
setGeneric("vrFeatureData<-", function(object, ..., value) standardGeneric("vrFeatureData<-"))

#' getFeatures
#'
#' Get variable features of assays
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname getFeatures
#' @export getFeatures
setGeneric("getFeatures", function(object, ...) standardGeneric("getFeatures"))

#' vrData
#'
#' Get data of assays
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrData
#' @export vrData
#' @order 1
setGeneric("vrData", function(object, ...) standardGeneric("vrData"))

#' generateTileData
#'
#' Generating data matrices for tile-based VoltRon objects from images
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname generateTileData
#' @export generateTileData
#' @order 1
setGeneric("generateTileData", function(object, ...) standardGeneric("generateTileData"))

#' changeSampleNames
#'
#' change sample names of VoltRon or other objects
#'
#' @param object a VoltRon or vrMetadata object.
#' @param ... arguments passed to other methods.
#'
#' @rdname changeSampleNames
#'
#' @noRd
setGeneric("changeSampleNames", function(object, ...) standardGeneric("changeSampleNames"))

#' changeAssayNames
#'
#' change assay names of VoltRon or other objects
#'
#' @param object a VoltRon, vrSample or vrLayer object.
#' @param ... arguments passed to other methods.
#'
#' @rdname changeAssayNames
#'
#' @noRd
setGeneric("changeAssayNames", function(object, ...) standardGeneric("changeAssayNames"))

####
# Assay ####
####

#' Main Assay
#'
#' Get and set the main assay of a VoltRon object
#' 
#' @param object a VoltRon object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrMainAssay
#' @export vrMainAssay
setGeneric("vrMainAssay", function(object, ...) standardGeneric("vrMainAssay"))

#' @param value new assay name
#' 
#' @rdname vrMainAssay
#' @export vrMainAssay<-
setGeneric("vrMainAssay<-", function(object, ..., value) standardGeneric("vrMainAssay<-"))

#' Update Assay
#'
#' update assays from vrAssay to vrAssayV2
#' 
#' @param object a VoltRon object.
#' @param ... arguments passed to other methods.
#'
#' @rdname updateAssay
#' @export updateAssay
setGeneric("updateAssay", function(object, ...) standardGeneric("updateAssay"))

#' Add Assay
#'
#' add a new assay (vrAssay object) to a VoltRon object
#'
#' @param object a VoltRon object.
#' @param ... arguments passed to other methods.
#'
#' @rdname addAssay
#' @export addAssay
setGeneric("addAssay", function(object, ...) standardGeneric("addAssay"))

#' Get Assay names
#'
#' Given a VoltRon, vrMetadata or vrAssay object, get/set names of assays.
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrAssayNames
#' @export vrAssayNames
#' @order 1
setGeneric("vrAssayNames", function(object, ...) standardGeneric("vrAssayNames"))

#' @rdname vrAssayNames
#' @noRd
setGeneric("vrAssayNames<-", function(object, ..., value) standardGeneric("vrAssayNames<-"))

#' Get assay types
#' 
#' Given a VoltRon or vrAssay object, get types of assays. 
#' Here, an assay type is of either tile, molecule, cell, spot or ROI.
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrAssayTypes
#' @export vrAssayTypes
#' @order 1
setGeneric("vrAssayTypes", function(object, ...) standardGeneric("vrAssayTypes"))

####
# Sample ####
####

#' Get Sample names
#'
#' Given a vrMetadata object, give names of samples
#'
#' @param object a vrMetadata object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrSampleNames
#' @export vrSampleNames
setGeneric("vrSampleNames", function(object, ...) standardGeneric("vrSampleNames"))

####
# Metadata ####
####

#' Metadata
#'
#' Get the metadata of a VoltRon object.
#'
#' @param object a VoltRon object.
#' @param ... arguments passed to other methods.
#'
#' @rdname Metadata
#' @export Metadata
setGeneric("Metadata", function(object, ...) standardGeneric("Metadata"))

#' @param value new metadata
#'
#' @rdname Metadata
setGeneric("Metadata<-", function(object, ..., value) standardGeneric("Metadata<-"))

####
# Processing ####
####

#' Normalize Data
#' 
#' Given a VoltRon or vrAssay object, normalize the raw count data.
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#' 
#' @rdname normalizeData
#' @export normalizeData
setGeneric("normalizeData", function(object, ...) standardGeneric("normalizeData"))

####
# Embedding ####
####

#' Get Embedding names
#'
#' Given a VoltRon or vrAssay object, give names of embeddings
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrEmbeddingNames
#' @export vrEmbeddingNames
#' @order 1
setGeneric("vrEmbeddingNames", function(object, ...) standardGeneric("vrEmbeddingNames"))

#' vrEmbeddings
#'
#' Given a VoltRon or vrAssay object, get embeddings of spatial points
#' 
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#' 
#' @rdname vrEmbeddings
#' @export vrEmbeddings
#' @order 1
setGeneric("vrEmbeddings", function(object, ...) standardGeneric("vrEmbeddings"))

#' @param value new embedding data
#'
#' @rdname vrEmbeddings
#' @export vrEmbeddings<-
setGeneric("vrEmbeddings<-", function(object, ..., value) standardGeneric("vrEmbeddings<-"))

####
# Feature ####
####

#' addFeature
#'
#' add a new feature set to a vrAssay (or VoltRon) object
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname addFeature
#' @export addFeature
#' @order 1
setGeneric("addFeature", function(object, ...) standardGeneric("addFeature"))

#' vrMainFeatureType
#'
#' Get the main feature type
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrMainFeatureType
#' @export vrMainFeatureType
#' @order 1
setGeneric("vrMainFeatureType", function(object, ...) standardGeneric("vrMainFeatureType"))

#' vrMainFeatureType
#'
#' Set the main image
#'
#' @param value the name of main feature set.
#'
#' @rdname vrMainFeatureType
#' @export vrMainFeatureType<-
setGeneric("vrMainFeatureType<-", function(object, ..., value) standardGeneric("vrMainFeatureType<-"))

#' vrFeatureTypeNames
#'
#' Get the names of the feature types
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrFeatureTypeNames
#' @export vrFeatureTypeNames
#' @order 1
setGeneric("vrFeatureTypeNames", function(object, ...) standardGeneric("vrFeatureTypeNames"))

####
# Spatial ####
####

#' vrSpatialPoints
#'
#' Get and set spatial entities.
#'
#' @param object a VoltRon, vrSample, vrLayer, vrAssay or vrSpatial object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrSpatialPoints
#' @export vrSpatialPoints
#' @order 1
#'
#' @examples
#' vrSpatialPoints(visium_data)
#' vrSpatialPoints(visium_data, assay = "Visium")
#' vrSpatialPoints(visium_data, assay = "Assay1")
setGeneric("vrSpatialPoints", function(object, ...) standardGeneric("vrSpatialPoints"))

#' @param value names for spatial points
#'
#' @rdname vrSpatialPoints
#' @export vrSpatialPoints<-
setGeneric("vrSpatialPoints<-", function(object, ..., value) standardGeneric("vrSpatialPoints<-"))

#' vrCoordinates
#'
#' Given a VoltRon, vrAssay or vrSpatial object, get and set the coordinates of assays and coordinate systems
#'
#' @param object a VoltRon, vrAssay or vrSpatial object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrCoordinates
#' @export vrCoordinates
#' @order 1
setGeneric("vrCoordinates", function(object, ...) standardGeneric("vrCoordinates"))

#' @param value coordinates
#' 
#' @rdname vrCoordinates
#' @export vrCoordinates<-
setGeneric("vrCoordinates<-", function(object, ..., value) standardGeneric("vrCoordinates<-"))

#' vrSegments
#'
#' Given a VoltRon, vrAssay or vrSpatial object, get and set the segment coordinates of assays and coordinate systems
#'
#' @param object a VoltRon, vrAssay or vrSpatial object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrSegments
#' @export vrSegments
#' @order 1
setGeneric("vrSegments", function(object, ...) standardGeneric("vrSegments"))

#' @param value segments
#' 
#' @rdname vrSegments
#' @export vrSegments<-
setGeneric("vrSegments<-", function(object, ..., value) standardGeneric("vrSegments<-"))

#' flipCoordinates
#'
#' Flip the coordinates of spatial points in the y axis direction. 
#'
#' @param object a VoltRon, vrAssay or vrSpatial object.
#' @param ... arguments passed to other methods.
#'
#' @rdname flipCoordinates
#' @export flipCoordinates
#' @order 1
setGeneric("flipCoordinates", function(object, ...) standardGeneric("flipCoordinates"))

####
# Image ####
####

#' vrImages
#'
#' Get images of VoltRon objects
#'
#' @param object a VoltRon, vrAssay or vrSpatial object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrImages
#' @export vrImages
#' @order 1
setGeneric("vrImages", function(object, ...) standardGeneric("vrImages"))

#' @param value a raster image or an image of \code{magick-image} object
#' 
#' @rdname vrImages
#' @export vrImages<-
setGeneric("vrImages<-", function(object, ..., value) standardGeneric("vrImages<-"))

#' vrImageNames
#'
#' Get names of all images
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrImageNames
#' @export vrImageNames
setGeneric("vrImageNames", function(object, ...) standardGeneric("vrImageNames"))

#' vrSpatialNames
#'
#' Get names of all spatial systems
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrSpatialNames
#' @export vrSpatialNames
setGeneric("vrSpatialNames", function(object, ...) standardGeneric("vrSpatialNames"))

#' vrImageChannelNames
#'
#' Get names of all image channels
#'
#' @param object a VoltRon, vrAssay or vrSpatial object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrImageChannelNames
#' @export vrImageChannelNames
#'
setGeneric("vrImageChannelNames", function(object, ...) standardGeneric("vrImageChannelNames"))

#' vrMainImage
#'
#' Get the main image
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrMainImage
#' @export vrMainImage
#' @order 1
setGeneric("vrMainImage", function(object, ...) standardGeneric("vrMainImage"))

#' vrMainSpatial
#'
#' Get the main spatial system name
#'
#' @param object a VoltRon or vrAssay object.
#' @param ... arguments passed to other methods.
#'
#' @rdname vrMainSpatial
#' @export vrMainSpatial
#' @order 1
setGeneric("vrMainSpatial", function(object, ...) standardGeneric("vrMainSpatial"))

#' vrMainImage
#'
#' Set the main image
#'
#' @param value the name of main spatial coordinate system.
#'
#' @rdname vrMainImage
#' @export vrMainImage<-
setGeneric("vrMainImage<-", function(object, ..., value) standardGeneric("vrMainImage<-"))

#' vrMainSpatial
#'
#' Set the main image
#'
#' @param value the name of main spatial coordinate system
#'
#' @rdname vrMainSpatial
#' @export vrMainSpatial<-
setGeneric("vrMainSpatial<-", function(object, ..., value) standardGeneric("vrMainSpatial<-"))

#' vrMainChannel
#'
#' Get and set the main channel name of the spatial system.
#' @param ... arguments passed to other methods.
#'
#' @param object a vrAssay or vrSpatial object
#'
#' @rdname vrMainChannel
#' @export vrMainChannel
#' @order 1
setGeneric("vrMainChannel", function(object, ...) standardGeneric("vrMainChannel"))

#' @param value the name of main channel of the spatial system
#'
#' @rdname vrMainChannel
#' @export vrMainChannel<-
setGeneric("vrMainChannel<-", function(object, ..., value) standardGeneric("vrMainChannel<-"))

#' resizeImage
#'
#' Resizing Magick images
#'
#' @param object a VoltRon, vrAssay or vrSpatial object.
#' @param ... arguments passed to other methods.
#'
#' @rdname resizeImage
#' @export resizeImage
#'
setGeneric("resizeImage", function(object, ...) standardGeneric("resizeImage"))

#' modulateImage
#'
#' Modulating Magick images
#'
#' @param object a VoltRon, vrAssay or vrSpatial object.
#' @param ... arguments passed to other methods.
#'
#' @rdname modulateImage
#' @export modulateImage
#'
setGeneric("modulateImage", function(object, ...) standardGeneric("modulateImage"))

#' combineChannels
#'
#' Combining channels into novel channels of the same image
#'
#' @param object a VoltRon, vrAssay or vrSpatial object.
#' @param ... arguments passed to other methods.
#'
#' @rdname combineChannels
#' @export combineChannels
#'
setGeneric("combineChannels", function(object, ...) standardGeneric("combineChannels"))

####
# Conversion ####
####

#' as.VoltRon
#'
#' Generic methods for conversion into a VoltRon object
#'
#' @param object a Seurat or SpatialExperiment object
#' @param ... Arguments passed to other methods
#'
#' @rdname as.VoltRon
#' @export as.VoltRon
as.VoltRon <- function(object, ...) {
  UseMethod(generic = 'as.VoltRon', object = object)
}

#' as.Zarr
#'
#' Generic methods to save VoltRon or magick-image objects as zarr files
#'
#' @param object a VoltRon or magick-image object
#' @param out_path output path to zarr file
#' @param image_id image name
#'
#' @rdname as.Zarr
#' @export as.Zarr
#'
as.Zarr <- function(object, out_path, image_id) {
  UseMethod(generic = 'as.Zarr', object = object)
}

####
# OnDisk ####
####

#' .shorten_assay_links
#'
#' @param object a vrAssay or vrAssayV2 object.
#' @param ... arguments passed to other methods.
#'
#' @noRd
setGeneric(".shorten_assay_links", function(object, ...) standardGeneric(".shorten_assay_links"))

#' .restore_absolute_assay_links
#'
#' @param object a vrAssay or vrAssayV2 object.
#' @param ... arguments passed to other methods.
#'
#' @noRd
setGeneric(".restore_absolute_assay_links", function(object, ...) standardGeneric(".restore_absolute_assay_links"))

#' .restore_absolute_links
#'
#' @param x an object of either DelayedArray or other variants
#' @param ... arguments passed to other methods.
#'
#' @noRd
setGeneric(".restore_absolute_links", function(x, ...) standardGeneric(".restore_absolute_links"))