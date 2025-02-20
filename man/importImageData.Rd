% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{importImageData}
\alias{importImageData}
\title{importImageData}
\usage{
importImageData(
  image,
  tile.size = 10,
  segments = NULL,
  image_name = "main",
  channel_names = NULL,
  ...
)
}
\arguments{
\item{image}{a single or a list of image paths or magick-image objects}

\item{tile.size}{the size of tiles}

\item{segments}{Either a list of segments or a GeoJSON file. This will result in a second assay in the VoltRon object to be created}

\item{image_name}{the image name of the Image assay, Default: main}

\item{channel_names}{the channel names of the images if multiple images are provided}

\item{...}{additional parameters passed to \link{formVoltRon}}
}
\description{
import an image as VoltRon object
}
\examples{
# single image
imgfile <- system.file("extdata", "DAPI.tif", package = "VoltRon")
vrdata <- importImageData(imgfile, image_name = "main")

# multiple images
imgfile <- c(system.file("extdata", "DAPI.tif", package = "VoltRon"), 
             system.file("extdata", "DAPI.tif", package = "VoltRon"))
vrdata <- importImageData(imgfile, image_name = "main", channel_name = c("DAPI", "DAPI2"))

}
