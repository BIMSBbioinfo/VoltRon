% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assay.R
\name{formAssay}
\alias{formAssay}
\title{formAssay}
\usage{
formAssay(
  data = NULL,
  coords,
  segments = list(),
  image = NULL,
  params = list(),
  type = "ROI",
  name = "Assay1",
  main_image = "image_1",
  ...
)
}
\arguments{
\item{data}{the count table}

\item{coords}{the coordinates of the spatial points}

\item{segments}{the segments of the spatial points (optional)}

\item{image}{the image of the data}

\item{params}{additional parameters of the object}

\item{type}{the type of the assay (cells, spots, ROIs)}

\item{name}{the name of the assay}

\item{main_image}{the name of the main_image}

\item{...}{additional arguements passed to \code{formImage}}
}
\description{
Create a vrAssay (VoltRon assay) object
}
