% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/image.R
\name{modulateImage}
\alias{modulateImage}
\alias{modulateImage.VoltRon}
\alias{modulateImage.vrAssay}
\title{modulateImage}
\usage{
modulateImage(object, ...)

\method{modulateImage}{VoltRon}(object, ...)

\method{modulateImage}{vrAssay}(object, brightness = 100, saturation = 100, hue = 100)
}
\arguments{
\item{object}{An object}

\item{...}{Arguments passed to other methods}

\item{brightness}{modulation of brightness as percentage of the current value (100 for no change)}

\item{saturation}{modulation of saturation as percentage of the current value (100 for no change)}

\item{hue}{modulation of hue is an absolute rotation of -180 degrees to +180 degrees from the current position corresponding to an argument range of 0 to 200 (100 for no change)}
}
\value{
\code{modulateImage}: The name of the default assay
}
\description{
Modulating Magick images
}