#include "Rcpp.h"
#include <optional>
#include <opencv2/opencv.hpp>

#ifndef MATTE_MI_H
#define MATTE_MI_H

struct IntensityRange {
  double min;
  double max;
};

double cubicBSpline(double u);

double scaleToBinPosition(double value, IntensityRange range, 
                          double low, double high);

std::size_t roundToNearestEvenNonnegative(double x);

bool isValidRange(IntensityRange range);

double mattesMiFromValues(const double* fixedValues,
                          const double* movingValues,
                          std::size_t count,
                          std::size_t bins,
                          std::optional<IntensityRange> fixedRange = std::nullopt,
                          std::optional<IntensityRange> movingRange = std::nullopt);

cv::Mat1d MatteMIMap(const cv::Mat& fixed, const cv::Mat& moving,
                     const cv::Mat& mask, int bins);

double MatteMI(const cv::Mat& fixed, const cv::Mat& moving,
               const cv::Mat& mask, int bins);
  
#endif