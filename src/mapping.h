#include "Rcpp.h"
#include <opencv2/opencv.hpp>
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/shape/shape_transformer.hpp"

#ifndef MAPPING_H
#define MAPPING_H

Rcpp::NumericMatrix applyMapping(Rcpp::NumericMatrix coords, Rcpp::List mapping);

#endif