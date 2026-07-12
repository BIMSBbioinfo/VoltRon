#include "Rcpp.h"
#include <opencv2/opencv.hpp>

#ifndef METRICS_H
#define METRICS_H

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
                          std::size_t bins = 64,
                          std::optional<IntensityRange> fixedRange = std::nullopt,
                          std::optional<IntensityRange> movingRange = std::nullopt);
  
#endif