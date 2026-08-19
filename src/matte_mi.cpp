#include "matte_mi.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

#include <opencv2/opencv.hpp>

namespace {

using Pixel = unsigned char;

struct IntensityRange {
  double min;
  double max;
};

struct ChunkSize {
  int height = 50;
  int width = 50;
};

struct GlobalCounts {
  std::array<std::size_t, 256> fixed{};
  std::array<std::size_t, 256> moving{};
  std::size_t validPairs = 0U;
};

// Source pixels remain CV_8U. All Mattes arithmetic remains double.
double cubicBSpline(double u) noexcept {
  u = std::abs(u);
  if (u < 1.0) {
    const double u2 = u * u;
    return (4.0 - 6.0 * u2 + 3.0 * u2 * u) / 6.0;
  }
  if (u < 2.0) {
    const double t = 2.0 - u;
    return t * t * t / 6.0;
  }
  return 0.0;
}

double scaleToBinPosition(
    double value,
    IntensityRange range,
    double low,
    double high) {
  if (!std::isfinite(range.min) ||
      !std::isfinite(range.max) ||
      !(range.max > range.min)) {
    throw std::invalid_argument("Invalid intensity range.");
  }
  value = std::clamp(value, range.min, range.max);
  return low + (value - range.min) * (high - low) /
    (range.max - range.min);
}

std::size_t roundToNearestEvenNonnegative(double x) noexcept {
  const double lowerDouble = std::floor(x);
  const double fraction = x - lowerDouble;
  const auto lower = static_cast<std::size_t>(lowerDouble);
  if (fraction < 0.5) return lower;
  if (fraction > 0.5) return lower + 1U;
  return (lower % 2U == 0U) ? lower : lower + 1U;
}

bool isValidRange(IntensityRange range) noexcept {
  return std::isfinite(range.min) &&
    std::isfinite(range.max) &&
    range.max > range.min;
}

int ceilDividePositive(int value, int divisor) noexcept {
  return value / divisor + ((value % divisor) != 0 ? 1 : 0);
}

void validateInputs(
    const cv::Mat& fixed,
    const cv::Mat& moving,
    const cv::Mat& mask,
    int bins) {
  if (fixed.empty() || moving.empty()) {
    throw std::invalid_argument(
      "fixed and moving images must not be empty.");
  }
  if (fixed.type() != CV_8UC1 || moving.type() != CV_8UC1) {
    throw std::invalid_argument(
      "The zero-conversion implementation expects fixed and moving "
      "to be CV_8UC1.");
  }
  if (fixed.size() != moving.size()) {
    throw std::invalid_argument(
      "fixed and moving must have the same dimensions.");
  }
  if (!mask.empty() &&
      (mask.type() != CV_8UC1 || mask.size() != fixed.size())) {
    throw std::invalid_argument(
      "mask must be empty or CV_8UC1 with the same dimensions.");
  }
  if (bins < 4) {
    throw std::invalid_argument(
      "bins must be >= 4 for cubic B-spline smoothing.");
  }
}

GlobalCounts collectGlobalCounts(
    const cv::Mat& fixed,
    const cv::Mat& moving,
    const cv::Mat& mask) {
  GlobalCounts out;
  for (int y = 0; y < fixed.rows; ++y) {
    const Pixel* fixedRow = fixed.ptr<Pixel>(y);
    const Pixel* movingRow = moving.ptr<Pixel>(y);
    const Pixel* maskRow = mask.empty() ? nullptr : mask.ptr<Pixel>(y);
    for (int x = 0; x < fixed.cols; ++x) {
      if (maskRow != nullptr && maskRow[x] == 0U) continue;
      ++out.fixed[fixedRow[x]];
      ++out.moving[movingRow[x]];
      ++out.validPairs;
    }
  }
  return out;
}

Pixel valueAtRank(
    const std::array<std::size_t, 256>& counts,
    std::size_t rank) {
  std::size_t cumulative = 0U;
  for (std::size_t value = 0U; value < counts.size(); ++value) {
    cumulative += counts[value];
    if (rank < cumulative) return static_cast<Pixel>(value);
  }
  throw std::out_of_range("Percentile rank is out of range.");
}

// Exact NumPy-style linear percentile for CV_8U values, without sorting pixels.
double percentileFromCounts(
    const std::array<std::size_t, 256>& counts,
    std::size_t count,
    double percentile) {
  if (count == 0U) {
    throw std::invalid_argument("Cannot calculate an empty percentile.");
  }
  if (!std::isfinite(percentile) || percentile < 0.0 || percentile > 100.0) {
    throw std::invalid_argument("Percentile must be in [0, 100].");
  }
  if (count == 1U) return static_cast<double>(valueAtRank(counts, 0U));

  const double index = static_cast<double>(count - 1U) * percentile / 100.0;
  const auto lowerIndex = static_cast<std::size_t>(std::floor(index));
  const auto upperIndex = static_cast<std::size_t>(std::ceil(index));
  const double fraction = index - static_cast<double>(lowerIndex);
  const double lower = static_cast<double>(valueAtRank(counts, lowerIndex));
  const double upper = static_cast<double>(valueAtRank(counts, upperIndex));
  return lower + fraction * (upper - lower);
}

IntensityRange minMaxRangeFromCounts(
    const std::array<std::size_t, 256>& counts,
    std::size_t count) {
  if (count == 0U) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return {nan, nan};
  }
  std::size_t minimum = 0U;
  while (minimum < counts.size() && counts[minimum] == 0U) ++minimum;
  std::size_t maximum = counts.size() - 1U;
  while (maximum > 0U && counts[maximum] == 0U) --maximum;
  return {
    static_cast<double>(minimum),
    static_cast<double>(maximum)
  };
}

void addMattesPair(
    double fixedValue,
    double movingValue,
    IntensityRange fixedRange,
    IntensityRange movingRange,
    std::size_t bins,
    std::vector<double>& jointHistogram) {
  const double fixedPosition = scaleToBinPosition(
    fixedValue, fixedRange, 0.0, static_cast<double>(bins - 1U));
  std::size_t fixedBin = roundToNearestEvenNonnegative(fixedPosition);
  fixedBin = std::min(fixedBin, bins - 1U);

  const double movingPosition = scaleToBinPosition(
    movingValue, movingRange, 1.0, static_cast<double>(bins - 2U));
  const auto baseBin = static_cast<std::ptrdiff_t>(std::floor(movingPosition));

  for (int offset = -1; offset <= 2; ++offset) {
    const std::ptrdiff_t movingBin = baseBin + offset;
    if (movingBin < 0 || movingBin >= static_cast<std::ptrdiff_t>(bins)) {
      continue;
    }
    const double weight = cubicBSpline(
      movingPosition - static_cast<double>(movingBin));
    if (weight <= 0.0) continue;
    jointHistogram[
      fixedBin * bins + static_cast<std::size_t>(movingBin)
    ] += weight;
  }
}

double mutualInformationFromHistogram(
    const std::vector<double>& jointHistogram,
    std::size_t bins) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double total = std::accumulate(
    jointHistogram.begin(), jointHistogram.end(), 0.0);
  if (!(total > 0.0) || !std::isfinite(total)) return nan;

  std::vector<double> px(bins, 0.0);
  std::vector<double> py(bins, 0.0);

  for (std::size_t fixedBin = 0; fixedBin < bins; ++fixedBin) {
    for (std::size_t movingBin = 0; movingBin < bins; ++movingBin) {
      const double pxy = jointHistogram[
        fixedBin * bins + movingBin
      ] / total;
      px[fixedBin] += pxy;
      py[movingBin] += pxy;
    }
  }

  double mi = 0.0;
  for (std::size_t fixedBin = 0; fixedBin < bins; ++fixedBin) {
    for (std::size_t movingBin = 0; movingBin < bins; ++movingBin) {
      const double pxy = jointHistogram[
        fixedBin * bins + movingBin
      ] / total;
      const double pxPy = px[fixedBin] * py[movingBin];
      if (pxy > 0.0 && pxPy > 0.0) {
        mi += pxy * std::log(pxy / pxPy);
      }
    }
  }
  return mi;
}

double mattesMiFromValues(
    const Pixel* fixedValues,
    const Pixel* movingValues,
    std::size_t count,
    std::size_t bins,
    std::optional<IntensityRange> fixedRange = std::nullopt,
    std::optional<IntensityRange> movingRange = std::nullopt) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  if (count != 0U && (fixedValues == nullptr || movingValues == nullptr)) {
    throw std::invalid_argument("Input value pointer is null.");
  }
  if (count < 2U) return nan;
  if (bins < 4U) {
    throw std::invalid_argument(
      "bins must be >= 4 for cubic B-spline smoothing.");
  }
  if (bins > std::numeric_limits<std::size_t>::max() / bins) {
    throw std::length_error("Histogram dimensions are too large.");
  }

  if (!fixedRange.has_value() || !movingRange.has_value()) {
    std::array<std::size_t, 256> fixedCounts{};
    std::array<std::size_t, 256> movingCounts{};
    for (std::size_t i = 0; i < count; ++i) {
      ++fixedCounts[fixedValues[i]];
      ++movingCounts[movingValues[i]];
    }
    if (!fixedRange.has_value()) {
      fixedRange = minMaxRangeFromCounts(fixedCounts, count);
    }
    if (!movingRange.has_value()) {
      movingRange = minMaxRangeFromCounts(movingCounts, count);
    }
  }

  if (!isValidRange(*fixedRange) || !isValidRange(*movingRange)) return nan;

  std::vector<double> jointHistogram(bins * bins, 0.0);
  for (std::size_t i = 0; i < count; ++i) {
    addMattesPair(
      static_cast<double>(fixedValues[i]),
      static_cast<double>(movingValues[i]),
      *fixedRange,
      *movingRange,
      bins,
      jointHistogram);
  }
  return mutualInformationFromHistogram(jointHistogram, bins);
}

} // namespace

cv::Mat1d MatteMIMap(
    const cv::Mat& fixed,
    const cv::Mat& moving,
    const cv::Mat& mask,
    int bins) {
  validateInputs(fixed, moving, mask, bins);

  constexpr ChunkSize chunkSize{};
  constexpr std::size_t minValidPixels = 100U;
  constexpr double minValidFraction = 0.10;
  constexpr double lowerPercentile = 0.5;
  constexpr double upperPercentile = 99.5;

  const GlobalCounts global = collectGlobalCounts(fixed, moving, mask);
  if (global.validPairs == 0U) {
    throw std::invalid_argument("The mask contains no valid pixels.");
  }

  const IntensityRange fixedRange{
    percentileFromCounts(global.fixed, global.validPairs, lowerPercentile),
    percentileFromCounts(global.fixed, global.validPairs, upperPercentile)
  };
  const IntensityRange movingRange{
    percentileFromCounts(global.moving, global.validPairs, lowerPercentile),
    percentileFromCounts(global.moving, global.validPairs, upperPercentile)
  };

  if (!isValidRange(fixedRange)) {
    throw std::invalid_argument("Invalid fixed intensity range.");
  }
  if (!isValidRange(movingRange)) {
    throw std::invalid_argument("Invalid moving intensity range.");
  }

  const int nRows = ceilDividePositive(fixed.rows, chunkSize.height);
  const int nCols = ceilDividePositive(fixed.cols, chunkSize.width);
  cv::Mat1d nmiMap(nRows, nCols);
  nmiMap.setTo(cv::Scalar(std::numeric_limits<double>::quiet_NaN()));

  const std::size_t maxChunkPixels =
    static_cast<std::size_t>(std::min(chunkSize.height, fixed.rows)) *
    static_cast<std::size_t>(std::min(chunkSize.width, fixed.cols));

  std::vector<Pixel> fixedChunkValues;
  std::vector<Pixel> movingChunkValues;
  fixedChunkValues.reserve(maxChunkPixels);
  movingChunkValues.reserve(maxChunkPixels);

  for (int row = 0; row < nRows; ++row) {
    for (int col = 0; col < nCols; ++col) {
      const int y0 = row * chunkSize.height;
      const int x0 = col * chunkSize.width;
      const int y1 = y0 + std::min(chunkSize.height, fixed.rows - y0);
      const int x1 = x0 + std::min(chunkSize.width, fixed.cols - x0);
      const std::size_t totalPixels =
        static_cast<std::size_t>(y1 - y0) *
        static_cast<std::size_t>(x1 - x0);

      fixedChunkValues.clear();
      movingChunkValues.clear();

      for (int y = y0; y < y1; ++y) {
        const Pixel* fixedRow = fixed.ptr<Pixel>(y);
        const Pixel* movingRow = moving.ptr<Pixel>(y);
        const Pixel* maskRow = mask.empty() ? nullptr : mask.ptr<Pixel>(y);
        for (int x = x0; x < x1; ++x) {
          if (maskRow != nullptr && maskRow[x] == 0U) continue;
          fixedChunkValues.push_back(fixedRow[x]);
          movingChunkValues.push_back(movingRow[x]);
        }
      }

      const std::size_t validPixels = fixedChunkValues.size();
      if (validPixels < minValidPixels) continue;

      const double validFraction =
        static_cast<double>(validPixels) /
        static_cast<double>(totalPixels);
      if (validFraction < minValidFraction) continue;

      nmiMap(row, col) = mattesMiFromValues(
        fixedChunkValues.data(),
        movingChunkValues.data(),
        validPixels,
        static_cast<std::size_t>(bins),
        std::optional<IntensityRange>{fixedRange},
        std::optional<IntensityRange>{movingRange});
    }
  }

  return nmiMap;
}

double MatteMI(
    const cv::Mat& fixed,
    const cv::Mat& moving,
    const cv::Mat& mask,
    int bins) {
  validateInputs(fixed, moving, mask, bins);

  const std::size_t binCount = static_cast<std::size_t>(bins);
  if (binCount > std::numeric_limits<std::size_t>::max() / binCount) {
    throw std::length_error("Histogram dimensions are too large.");
  }

  const GlobalCounts global = collectGlobalCounts(fixed, moving, mask);
  if (global.validPairs < 2U) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const IntensityRange fixedRange =
    minMaxRangeFromCounts(global.fixed, global.validPairs);
  const IntensityRange movingRange =
    minMaxRangeFromCounts(global.moving, global.validPairs);
  if (!isValidRange(fixedRange) || !isValidRange(movingRange)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::vector<double> jointHistogram(binCount * binCount, 0.0);
  for (int y = 0; y < fixed.rows; ++y) {
    const Pixel* fixedRow = fixed.ptr<Pixel>(y);
    const Pixel* movingRow = moving.ptr<Pixel>(y);
    const Pixel* maskRow = mask.empty() ? nullptr : mask.ptr<Pixel>(y);
    for (int x = 0; x < fixed.cols; ++x) {
      if (maskRow != nullptr && maskRow[x] == 0U) continue;
      addMattesPair(
        static_cast<double>(fixedRow[x]),
        static_cast<double>(movingRow[x]),
        fixedRange,
        movingRange,
        binCount,
        jointHistogram);
    }
  }

  return mutualInformationFromHistogram(jointHistogram, binCount);
}
