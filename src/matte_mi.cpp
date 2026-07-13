#include <Rcpp.h>

// OpenCV
#include <opencv2/opencv.hpp>

// Internal functions
#include "auxiliary.h"
#include "image.h"

// Namespaces
using namespace Rcpp;
using namespace std;
using namespace cv;

struct IntensityRange {
  double min;
  double max;
};

////
// Chunk wise Matte MI
////

double cubicBSpline(double u) {
  u = std::abs(u);
  
  if (u < 1.0) {
    const double u2 = u * u;
    const double u3 = u2 * u;
    
    return (4.0 - 6.0 * u2 + 3.0 * u3) / 6.0;
  }
  
  if (u < 2.0) {
    const double t = 2.0 - u;
    return (t * t * t) / 6.0;
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
  
  return low +
    (value - range.min) * (high - low) /
      (range.max - range.min);
}

std::size_t roundToNearestEvenNonnegative(double x) {
  const double lowerAsDouble = std::floor(x);
  const double fraction = x - lowerAsDouble;
  const auto lower = static_cast<std::size_t>(lowerAsDouble);
  
  if (fraction < 0.5) {
    return lower;
  }
  
  if (fraction > 0.5) {
    return lower + 1U;
  }
  
  return (lower % 2U == 0U) ? lower : lower + 1U;
}

bool isValidRange(IntensityRange range) {
  return std::isfinite(range.min) &&
    std::isfinite(range.max) &&
    range.max > range.min;
}

/**
 * Compute Mattes-style mutual information from paired sample values.
 *
 * Fixed samples use nearest-bin assignment.
 * Moving samples use cubic B-spline Parzen smoothing.
 *
 * @param fixedValues   Pointer to fixed-image sample values.
 * @param movingValues  Pointer to corresponding moving-image values.
 * @param countPair     Number of paired values.
 * @param bins          Number of histogram bins; must be >= 4.
 * @param fixedRange    Optional fixed intensity range.
 * @param movingRange   Optional moving intensity range.
 *
 * @return Mutual information in nats, or quiet NaN for degenerate data.
 */
double mattesMiFromValues(
    const double* fixedValues,
    const double* movingValues,
    std::size_t countPair,
    std::size_t bins = 64,
    std::optional<IntensityRange> fixedRange = std::nullopt,
    std::optional<IntensityRange> movingRange = std::nullopt) {
  
  const double nan =
    std::numeric_limits<double>::quiet_NaN();
  
  if (countPair != 0U &&
      (fixedValues == nullptr || movingValues == nullptr)) {
    throw std::invalid_argument("Input value pointer is null.");
  }
  
  /*
   * First pass:
   *   - count finite sample pairs
   *   - calculate automatic intensity ranges
   */
  std::size_t validCount = 0;
  
  double fixedMin =  std::numeric_limits<double>::infinity();
  double fixedMax = -std::numeric_limits<double>::infinity();
    
  double movingMin = std::numeric_limits<double>::infinity();
  double movingMax = -std::numeric_limits<double>::infinity();
  
  for (std::size_t i = 0; i < countPair; ++i) {
    
    if (!std::isfinite(fixedValues[i]) ||
        !std::isfinite(movingValues[i])) {
        continue;
    }
    
    ++validCount;
    
    fixedMin = std::min(fixedMin, fixedValues[i]);
    fixedMax = std::max(fixedMax, fixedValues[i]);
    
    movingMin = std::min(movingMin, movingValues[i]);
    movingMax = std::max(movingMax, movingValues[i]);
  }
  
  // Preserve the Python function's validation order.
  if (validCount < 2U) {
    return nan;
  }
  
  if (bins < 4U) {
    throw std::invalid_argument(
        "bins must be >= 4 for cubic B-spline smoothing.");
  }
  
  if (bins >
        std::numeric_limits<std::size_t>::max() / bins ||
        bins >
        static_cast<std::size_t>(
          std::numeric_limits<std::ptrdiff_t>::max())) {
    throw std::length_error(
        "Histogram dimensions are too large.");
  }
  
  /*
   * Row: fixed-image bin
   * Column: moving-image bin
   */
  std::vector<double> jointHistogram(
      bins * bins,
      0.0);
  
  /*
   * Second pass: construct the joint histogram.
   */
  for (std::size_t i = 0; i < countPair; ++i) {
    const double fixedValue = static_cast<double>(fixedValues[i]);
    const double movingValue = static_cast<double>(movingValues[i]);
    
    if (!std::isfinite(fixedValue) ||
        !std::isfinite(movingValue)) {
        continue;
    }
    
    /*
     * Fixed image:
     * map to [0, bins - 1] and use nearest-bin assignment.
     */
    const double fixedPosition =
      scaleToBinPosition(
        fixedValue,
        *fixedRange,
        0.0,
        static_cast<double>(bins - 1U));
    
    std::size_t fixedBin =
      roundToNearestEvenNonnegative(
        fixedPosition);
    
    fixedBin = std::min(
      fixedBin,
      bins - 1U);
    
    /*
     * Moving image:
     * map to [1, bins - 2] so the cubic kernel has
     * room at both ends of the histogram.
     */
    const double movingPosition =
      scaleToBinPosition(
        movingValue,
        *movingRange,
        1.0,
        static_cast<double>(bins - 2U));
    
    const auto baseBin =
      static_cast<std::ptrdiff_t>(
        std::floor(movingPosition));
    
    /*
     * A cubic B-spline contributes to at most four bins.
     */
    for (int offset = -1; offset <= 2; ++offset) {
      const std::ptrdiff_t movingBin =
        baseBin + offset;
      
      if (movingBin < 0 ||
          movingBin >=
          static_cast<std::ptrdiff_t>(bins)) {
        continue;
      }
      
      const double weight =
        cubicBSpline(movingPosition - static_cast<double>(movingBin));
      
      if (weight <= 0.0) {
        continue;
      }
      
      jointHistogram[
      fixedBin * bins +
        static_cast<std::size_t>(movingBin)
      ] += weight;
    }
  }
  
  const double total =
    std::accumulate(
      jointHistogram.begin(),
      jointHistogram.end(),
      0.0);
  
  if (!(total > 0.0)) {
    return nan;
  }
  
  /*
   * Marginal distributions.
   */
  std::vector<double> px(bins, 0.0);
  std::vector<double> py(bins, 0.0);
  
  for (std::size_t fixedBin = 0;
       fixedBin < bins;
       ++fixedBin) {
    
    for (std::size_t movingBin = 0;
         movingBin < bins;
         ++movingBin) {
      
      const double pxy =
      jointHistogram[
      fixedBin * bins + movingBin
      ] / total;
      
      px[fixedBin] += pxy;
      py[movingBin] += pxy;
    }
  }
  
  /*
   * MI = sum p(x,y) log(p(x,y) / (p(x)p(y)))
   */
  double mi = 0.0;
  
  for (std::size_t fixedBin = 0;
       fixedBin < bins;
       ++fixedBin) {
    
    for (std::size_t movingBin = 0;
         movingBin < bins;
         ++movingBin) {
      
      const double pxy =
      jointHistogram[
      fixedBin * bins + movingBin
      ] / total;
      
      const double pxPy =
        px[fixedBin] * py[movingBin];
      
      if (pxy > 0.0 && pxPy > 0.0) {
        mi +=
          pxy * std::log(pxy / pxPy);
      }
    }
  }
  
  // std::log is the natural logarithm, so MI is in nats.
  return mi;
}

////
// Matte MI
////

/**
 * Spatial chunk dimensions.
 *
 * The field order intentionally follows the Python API:
 *     ChunkSize{height, width}
 *
 * This avoids cv::Size's opposite (width, height) ordering.
 */
struct ChunkSize {
  int height = 50;
  int width = 50;
};

/**
 * Output of chunkedNmiMap().
 *
 * Despite the legacy nmiMap name, the values are Mattes-style mutual
 * information values, not normalized mutual information values.
 *
 * Matrix layouts:
 *   nmiMap  : CV_64FC1, shape [chunk rows, chunk columns]
 *   bounds  : CV_32SC4, each element is [y0, y1, x0, x1]
 *   centers : CV_64FC2, each element is [y center, x center]
 */
struct ChunkedNmiMapResult {
  cv::Mat1d nmiMap;
  cv::Mat_<cv::Vec4i> bounds;
  cv::Mat_<cv::Vec2d> centers;
};

int ceilDividePositive(int value, int divisor) noexcept {
  return value / divisor + ((value % divisor) != 0 ? 1 : 0);
}

double linearPercentileFromSorted(
    const std::vector<double>& sortedValues,
    double percentile) {
  if (sortedValues.empty()) {
    throw std::invalid_argument(
        "Cannot calculate a percentile of an empty array.");
  }
  
  if (!std::isfinite(percentile) ||
      percentile < 0.0 || percentile > 100.0) {
    throw std::invalid_argument(
        "Percentile must be finite and in [0, 100].");
  }
  
  if (sortedValues.size() == 1U) {
    return sortedValues.front();
  }
  
  // Matches NumPy's default linear percentile interpolation:
  // index = (N - 1) * percentile / 100.
  const double index =
    (static_cast<double>(sortedValues.size() - 1U) * percentile) /
      100.0;
  
  const auto lowerIndex =
    static_cast<std::size_t>(std::floor(index));
  const auto upperIndex =
    static_cast<std::size_t>(std::ceil(index));
  const double fraction = index - static_cast<double>(lowerIndex);
  
  const double lower = sortedValues[lowerIndex];
  const double upper = sortedValues[upperIndex];
  
  return lower + fraction * (upper - lower);
}

void validateChunkedMatteMIInputs(
    cv::Mat1b& validGlobal,
    std::size_t& validGlobalCount,
    const cv::Mat& fixed,
    const cv::Mat& moving,
    const cv::Mat& mask,
    ChunkSize chunkSize,
    int bins) {
  
  const int height = fixed.rows;
  const int width = fixed.cols;
  
  for (int y = 0; y < height; ++y) {
    const double* fixedRow = fixed.ptr<double>(y);
    const double* movingRow = moving.ptr<double>(y);
    const double* maskRow =
      mask.empty() ? nullptr : mask.ptr<double>(y);
    unsigned char* validRow = validGlobal.ptr<unsigned char>(y);
    
    for (int x = 0; x < width; ++x) {
      const bool insideMask =
        maskRow == nullptr || static_cast<bool>(maskRow[x]);
      
      const bool valid =
        insideMask &&
        std::isfinite(fixedRow[x]) &&
        std::isfinite(movingRow[x]);
      
      if (valid) {
        validRow[x] = 1U;
        ++validGlobalCount;
      }
    }
  }

  // stop if no pixels are valid  
  if (validGlobalCount == 0U) {
    throw std::invalid_argument(
        "The mask contains no valid pixels.");
  }

}

// ChunkedNmiMapResult chunkedMatteMIMap(const cv::Mat& fixed,
cv::Mat1d MatteMIMap(const cv::Mat& fixed,
                     const cv::Mat& moving,
                     const cv::Mat& mask,
                     int bins = 50) {
                          
  ChunkSize chunkSize = ChunkSize{};

  // validate chunks and return validation map of pixels
  // const int height = fixed.rows;
  // const int width = fixed.cols;
  // cv::Mat1b validGlobal(height, width, static_cast<unsigned char>(0));
  // std::size_t validGlobalCount = 0U;
  // validateChunkedMatteMIInputs(validGlobal,
  //                              validGlobalCount,
  //                              fixed, 
  //                              moving, 
  //                              mask, 
  //                              chunkSize, 
  //                              bins);
  
  
  // Do I need these to be cv_64f ? 
  cv::Mat fixed64;
  cv::Mat moving64;
  fixed.convertTo(fixed64, CV_64F);
  moving.convertTo(moving64, CV_64F);
  
  // Do I need these to be cv_64f ? 
  cv::Mat mask64;
  if (!mask.empty()) {
    mask.convertTo(mask64, CV_64F);
  }
  
  // validate
  const int height = fixed.rows;
  const int width = fixed.cols;
  
  cv::Mat1b validGlobal(height, width, static_cast<unsigned char>(0));
  std::size_t validGlobalCount = 0U;
  
  for (int y = 0; y < height; ++y) {
    const double* fixedRow = fixed64.ptr<double>(y);
    const double* movingRow = moving64.ptr<double>(y);
    const double* maskRow =
      mask64.empty() ? nullptr : mask64.ptr<double>(y);
    unsigned char* validRow = validGlobal.ptr<unsigned char>(y);
    
    for (int x = 0; x < width; ++x) {
      const bool insideMask =
        maskRow == nullptr || static_cast<bool>(maskRow[x]);
      
      const bool valid =
        insideMask &&
        std::isfinite(fixedRow[x]) &&
        std::isfinite(movingRow[x]);
      
      if (valid) {
        validRow[x] = 1U;
        ++validGlobalCount;
      }
    }
  }
  
  if (validGlobalCount == 0U) {
    throw std::invalid_argument(
        "The mask contains no valid pixels.");
  }
  
  // Calculate one fixed-image range and one moving-image range globally,
  // then reuse those ranges in every chunk. This makes chunk values
  // comparable across the image.
  std::vector<double> fixedGlobalValues;
  std::vector<double> movingGlobalValues;
  fixedGlobalValues.reserve(validGlobalCount);
  movingGlobalValues.reserve(validGlobalCount);
  
  for (int y = 0; y < height; ++y) {
    const double* fixedRow = fixed64.ptr<double>(y);
    const double* movingRow = moving64.ptr<double>(y);
    const double* maskRow =
      mask64.empty() ? nullptr : mask64.ptr<double>(y);
    const unsigned char* validRow =
      validGlobal.ptr<unsigned char>(y);
    
    for (int x = 0; x < width; ++x) {
      if (validRow[x] != 0U) {
        fixedGlobalValues.push_back(fixedRow[x]);
        movingGlobalValues.push_back(movingRow[x]);
      }
    }
  }
  
  std::sort(fixedGlobalValues.begin(), fixedGlobalValues.end());
  std::sort(movingGlobalValues.begin(), movingGlobalValues.end());
  
  constexpr double lowerPercentile = 0.5;
  constexpr double upperPercentile = 99.5;
  
  const IntensityRange fixedRange{
    linearPercentileFromSorted(fixedGlobalValues, lowerPercentile),
    linearPercentileFromSorted(fixedGlobalValues, upperPercentile)
  };
  
  // Intentional correction from the pasted Python: movingRange is derived
  // from moving-image values, not from fixed-image values.
  const IntensityRange movingRange{
    linearPercentileFromSorted(movingGlobalValues, lowerPercentile),
    linearPercentileFromSorted(movingGlobalValues, upperPercentile)
  };
  
  Rcout << fixedRange.min << " " << fixedRange.max << endl;
  if (!isValidRange(fixedRange)) {
    throw std::invalid_argument(
        "Invalid fixed intensity range.");
  }
  
  Rcout << movingRange.min << " " << movingRange.max << endl;
  if (!isValidRange(movingRange)) {
    throw std::invalid_argument(
        "Invalid moving intensity range.");
  }
  
  const int nRows = ceilDividePositive(height, chunkSize.height);
  const int nCols = ceilDividePositive(width, chunkSize.width);
  
  const double nan =
    std::numeric_limits<double>::quiet_NaN();
  
  // only set nmimap, why need bounds and centers
  cv::Mat1d NmiMap(nRows, nCols);
  NmiMap.setTo(cv::Scalar(nan));
  // result.nmiMap.setTo(cv::Scalar(nan));
  // result.bounds.setTo(cv::Scalar::all(0));
  // result.centers.setTo(cv::Scalar::all(0));
  
  constexpr std::size_t minValidPixels = 100U;
  constexpr double minValidFraction = 0.10;
  
  std::vector<double> fixedChunkValues;
  std::vector<double> movingChunkValues;
  
  for (int row = 0; row < nRows; ++row) {
    for (int col = 0; col < nCols; ++col) {
      const int y0 = row * chunkSize.height;
      const int x0 = col * chunkSize.width;
      
      // Written this way instead of y0 + chunkSize.height to avoid
      // signed integer overflow for extreme dimensions.
      const int y1 = y0 + std::min(chunkSize.height, height - y0);
      const int x1 = x0 + std::min(chunkSize.width, width - x0);
      
      const std::size_t totalPixels =
        static_cast<std::size_t>(y1 - y0) *
        static_cast<std::size_t>(x1 - x0);
      
      fixedChunkValues.clear();
      movingChunkValues.clear();
      fixedChunkValues.reserve(totalPixels);
      movingChunkValues.reserve(totalPixels);
      
      for (int y = y0; y < y1; ++y) {
        const double* fixedRow = fixed64.ptr<double>(y);
        const double* movingRow = moving64.ptr<double>(y);
        const unsigned char* validRow =
          validGlobal.ptr<unsigned char>(y);
        
        for (int x = x0; x < x1; ++x) {
          if (validRow[x] != 0U) {
            fixedChunkValues.push_back(fixedRow[x]);
            movingChunkValues.push_back(movingRow[x]);
          }
        }
      }
      
      const std::size_t validPixels = fixedChunkValues.size();
      
      if (validPixels < minValidPixels) {
        continue;
      }
      
      const double validFraction =
        static_cast<double>(validPixels) /
          static_cast<double>(totalPixels);
      
      if (validFraction < minValidFraction) {
        continue;
      }
      
      NmiMap(row, col) = mattesMiFromValues(
        fixedChunkValues.data(),
        movingChunkValues.data(),
        validPixels,
        static_cast<std::size_t>(bins),
        std::optional<IntensityRange>{fixedRange},
        std::optional<IntensityRange>{movingRange});
    }
  }
  
  return NmiMap;
}

const double* getRowAsDouble(
    const cv::Mat& image,
    int row,
    cv::Mat& scratch) {
  
  if (image.depth() == CV_64F) {
    return image.ptr<double>(row);
  }
  
  image.row(row).convertTo(scratch, CV_64F);
  return scratch.ptr<double>(0);
}

double MatteMI(
    const cv::Mat& fixed,
    const cv::Mat& moving,
    const cv::Mat& mask,
    int bins) {
  
  const double nan =
    std::numeric_limits<double>::quiet_NaN();
  
  if (bins < 4U) {
    throw std::invalid_argument(
        "bins must be >= 4 for cubic B-spline smoothing.");
  }
  
  const std::size_t binCount =
    static_cast<std::size_t>(bins);
  
  if (binCount >
        std::numeric_limits<std::size_t>::max() / binCount) {
    throw std::length_error(
        "Histogram dimensions are too large.");
  }
  
  /*
   * Convert the mask to a conventional CV_8U binary mask.
   *
   * Every nonzero mask value becomes 255.
   */
  cv::Mat1b mask8;
  
  if (!mask.empty()) {
    cv::compare(
      mask,
      cv::Scalar::all(0),
      mask8,
      cv::CMP_NE);
  }
  
  const int height = fixed.rows;
  const int width = fixed.cols;
  
  /*
   * Only row-sized double buffers are needed. This avoids converting
   * both complete images to CV_64F and avoids full-image value vectors.
   */
  cv::Mat fixedRowScratch;
  cv::Mat movingRowScratch;
  
  /*
   * Utility that visits every valid pixel pair.
   *
   * This function performs no sampling. Every valid pair is delivered
   * to the supplied visitor.
   */
  auto visitValidPairs = [&](auto&& visitor) {
    for (int y = 0; y < height; ++y) {
      const double* fixedRow =
        getRowAsDouble(
          fixed,
          y,
          fixedRowScratch);
      
      const double* movingRow =
        getRowAsDouble(
          moving,
          y,
          movingRowScratch);
      
      const unsigned char* maskRow =
        mask.empty()
        ? nullptr
          : mask8.ptr<unsigned char>(y);
      
      for (int x = 0; x < width; ++x) {
        if (maskRow != nullptr &&
            maskRow[x] == 0U) {
          continue;
        }
        
        const double fixedValue = fixedRow[x];
        const double movingValue = movingRow[x];
        
        if (!std::isfinite(fixedValue) ||
            !std::isfinite(movingValue)) {
            continue;
        }
        
        visitor(fixedValue, movingValue);
      }
    }
  };
  
  /*
   * First pass:
   * determine full-image intensity ranges from every valid pixel pair.
   *
   * This matches the default behavior of your original Python
   * _mattes_mi_from_values function when no explicit ranges are supplied.
   */
  std::size_t validCount = 0U;
  
  double fixedMin =
    std::numeric_limits<double>::infinity();
  
  double fixedMax =
    -std::numeric_limits<double>::infinity();
    
    double movingMin =
    std::numeric_limits<double>::infinity();
    
    double movingMax =
      -std::numeric_limits<double>::infinity();
      
      visitValidPairs(
        [&](double fixedValue, double movingValue) {
          ++validCount;
          
          fixedMin =
          std::min(fixedMin, fixedValue);
          
          fixedMax =
            std::max(fixedMax, fixedValue);
          
          movingMin =
            std::min(movingMin, movingValue);
          
          movingMax =
            std::max(movingMax, movingValue);
        });
      
      if (validCount < 2U) {
        return nan;
      }
      
      const IntensityRange fixedRange{
        fixedMin,
        fixedMax
      };
      
      const IntensityRange movingRange{
        movingMin,
        movingMax
      };
      
      /*
       * This follows the behavior of the supplied Python function:
       * a constant fixed or moving image has an invalid range and returns NaN.
       */
      if (!isValidRange(fixedRange) ||
          !isValidRange(movingRange)) {
          return nan;
      }
      
      /*
       * Joint histogram:
       *
       * rows    = fixed-image bins
       * columns = moving-image bins
       */
      std::vector<double> jointHistogram(
          binCount * binCount,
          0.0);
      
      /*
       * Second pass:
       * add every valid pixel pair to the Mattes joint histogram.
       */
      visitValidPairs(
        [&](double fixedValue, double movingValue) {
          /*
           * Fixed image:
           * nearest-bin assignment over [0, bins - 1].
           */
          const double fixedPosition =
            scaleToBinPosition(
              fixedValue,
              fixedRange,
              0.0,
              static_cast<double>(binCount - 1U));
          
          std::size_t fixedBin =
            roundToNearestEvenNonnegative(
              fixedPosition);
          
          fixedBin =
            std::min(
              fixedBin,
              binCount - 1U);
          
          /*
           * Moving image:
           * continuous coordinate over [1, bins - 2].
           *
           * The one-bin margin leaves room for the cubic B-spline
           * support at both histogram boundaries.
           */
          const double movingPosition =
            scaleToBinPosition(
              movingValue,
              movingRange,
              1.0,
              static_cast<double>(binCount - 2U));
          
          const std::ptrdiff_t baseBin =
            static_cast<std::ptrdiff_t>(
              std::floor(movingPosition));
          
          /*
           * Cubic B-spline support covers at most four bins.
           */
          for (int offset = -1;
               offset <= 2;
               ++offset) {
            
            const std::ptrdiff_t movingBin =
            baseBin +
            static_cast<std::ptrdiff_t>(offset);
            
            if (movingBin < 0 ||
                movingBin >=
                static_cast<std::ptrdiff_t>(
                  binCount)) {
              continue;
            }
            
            const double weight =
              cubicBSpline(
                movingPosition -
                  static_cast<double>(movingBin));
            
            if (weight <= 0.0) {
              continue;
            }
            
            const std::size_t histogramIndex =
              fixedBin * binCount +
              static_cast<std::size_t>(
                movingBin);
            
            jointHistogram[histogramIndex] +=
              weight;
          }
        });
      
      const double total =
        std::accumulate(
          jointHistogram.begin(),
          jointHistogram.end(),
          0.0);
      
      if (!(total > 0.0) ||
          !std::isfinite(total)) {
          return nan;
      }
      
      /*
       * Marginal probability distributions.
       */
      std::vector<double> px(
          binCount,
          0.0);
      
      std::vector<double> py(
          binCount,
          0.0);
      
      for (std::size_t fixedBin = 0;
           fixedBin < binCount;
           ++fixedBin) {
        
        for (std::size_t movingBin = 0;
             movingBin < binCount;
             ++movingBin) {
          
          const std::size_t index =
          fixedBin * binCount +
          movingBin;
          
          const double pxy =
            jointHistogram[index] / total;
          
          px[fixedBin] += pxy;
          py[movingBin] += pxy;
        }
      }
      
      /*
       * MI = sum p(x,y) log(p(x,y) / (p(x)p(y)))
       */
      double mi = 0.0;
      
      for (std::size_t fixedBin = 0;
           fixedBin < binCount;
           ++fixedBin) {
        
        for (std::size_t movingBin = 0;
             movingBin < binCount;
             ++movingBin) {
          
          const std::size_t index =
          fixedBin * binCount +
          movingBin;
          
          const double pxy =
            jointHistogram[index] / total;
          
          const double productOfMarginals =
            px[fixedBin] *
            py[movingBin];
          
          if (pxy > 0.0 &&
              productOfMarginals > 0.0) {
            
            mi +=
              pxy *
              std::log(
                pxy /
                  productOfMarginals);
          }
        }
      }
      
      // Natural logarithm: the result is in nats.
      return mi;
}