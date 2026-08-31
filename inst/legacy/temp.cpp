// ---------------------------------------------------------------------------
// tiled_metrics.cpp  --  single-pass tiled alignment metrics for VoltRon
//
// Computes histogram intersection, Bhattacharyya (Hellinger) distance and SSIM
// for every tile in one loop, then aggregates over tiles that overlap the mask.
//
// WHY ONE LOOP IS ENOUGH
//   For an 8-bit tile, a 256-bin histogram is a *complete* summary of that
//   tiles intensity distribution. So a single calcHist per image gives us:
//       - the valid pixel count   n  = sum(h)
//       - the mean                mu = sum(i*h)/n          (exact)
//       - the variance            var= sum(i^2*h)/n - mu^2 (exact)
//       - the histogram itself, which we rebin to nBins for the two
//         histogram metrics (summing groups of adjacent bins is exact)
//   The only thing a histogram cannot give us is the cross term E[XY], because
//   that depends on pixel correspondence. That costs one multiply + one mean.
//   Net cost per tile: 2 calcHist + 1 multiply + 1 mean. No redundant passes,
//   and no meanStdDev calls at all.
//
// NO CONVERSION OF THE INPUT IMAGES
//   Inputs stay CV_8U throughout. The one place a wider type is required is the
//   pixel-wise product, and cv::multiplys `dtype` argument handles that inline.
//   CV_32F is safe there because 8-bit products max out at 255*255 = 65025,
//   which float32 represents exactly (integers are exact below 2^24).
//   NEVER use `a - scalar` or plain `multiply` on CV_8U: both saturate, which
//   silently destroys the covariance term.
// ---------------------------------------------------------------------------

#include <Rcpp.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/check.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

// ---------------------------------------------------------------------------
// Exact n / mean / variance from a 256-bin histogram of an 8-bit tile.
// Bin i corresponds to intensity value i exactly, so nothing is approximated.
// The variance convention matches cv::meanStdDev (population, divide by n).
// ---------------------------------------------------------------------------
inline void histStats(const cv::Mat& h256, double& n, double& mu, double& var)
{
  double s0 = 0.0, s1 = 0.0, s2 = 0.0;
  const float* p = h256.ptr<float>(0);
  for (int i = 0; i < 256; ++i) {
    const double c = static_cast<double>(p[i]);
    s0 += c;
    s1 += static_cast<double>(i) * c;
    s2 += static_cast<double>(i) * static_cast<double>(i) * c;
  }
  n = s0;
  if (s0 <= 0.0) { mu = 0.0; var = 0.0; return; }
  mu  = s1 / s0;
  var = s2 / s0 - mu * mu;
  if (var < 0.0) var = 0.0;   // clamp the tiny negative that rounding can give
}

// ---------------------------------------------------------------------------
// Collapse a 256-bin histogram into nBins by summing adjacent groups, then
// L1-normalise so it is a probability distribution.
//
// WHY REBIN AT ALL: a 50x50 tile has only 2500 pixels. Spread over 256 bins
// that averages under 10 counts per bin, so most bins are near-empty and the
// comparison is dominated by sampling noise. 32 bins gives ~78 counts per bin.
// Fewer bins = more robust but coarser. nBins must divide 256.
//
// WHY L1-NORMALISE: intersection is scale-dependent -- without normalising it
// just reports pixel counts. After normalising, intersection lands in [0,1] and
// equals 1 - total-variation distance. It also simplifies OpenCV's
// Bhattacharyya, whose 1/sqrt(H1bar*H2bar*N^2) prefactor becomes exactly 1,
// leaving the clean Hellinger form sqrt(1 - sum(sqrt(h1*h2))).
// ---------------------------------------------------------------------------
inline cv::Mat rebinAndNormalise(const cv::Mat& h256, int nBins)
{
  const int group = 256 / nBins;
  cv::Mat out = cv::Mat::zeros(nBins, 1, CV_32F);
  const float* p = h256.ptr<float>(0);
  float*       q = out.ptr<float>(0);
  for (int i = 0; i < 256; ++i) q[i / group] += p[i];
  
  const double s = cv::sum(out)[0];
  if (s > 0.0) out /= s;
  return out;
}

} // namespace

// ---------------------------------------------------------------------------
// Result bundle. The three maps are coarse grids (one cell per tile) holding
// the per-tile value, or NaN where the tile was skipped. Keep them: a scalar
// tells you *whether* alignment is bad, the map tells you *where*.
// ---------------------------------------------------------------------------
struct TiledMetrics {
  double intersection  = kNaN;   // pixel-weighted mean, higher = better
  double bhattacharyya = kNaN;   // pixel-weighted mean, LOWER = better
  double ssim          = kNaN;   // pixel-weighted mean, higher = better
  double ssim_median   = kNaN;   // robust alternative, ignores outlier tiles
  
  int n_total          = 0;      // tiles in the grid
  int n_valid          = 0;      // tiles that passed both guards
  int n_drop_coverage  = 0;      // dropped: too little mask overlap
  int n_drop_flat      = 0;      // dropped: both tiles near-constant
  
  cv::Mat inter_map;             // CV_32F, nTilesY x nTilesX, NaN = skipped
  cv::Mat bhatta_map;
  cv::Mat ssim_map;
};

// ---------------------------------------------------------------------------
// tiledAlignmentMetrics
//
//   im1, im2     CV_8U, same size. Multi-channel input uses channel 0.
//   mask         CV_8U same size, non-zero = valid. Pass an empty Mat for all.
//   tile         tile edge in pixels (50 as requested)
//   nBins        histogram bins for the two histogram metrics (must divide 256)
//   minCoverage  fraction of a tile that must be inside the mask to count
//   varFactor    flat-tile rejection threshold, as a multiple of C2
// ---------------------------------------------------------------------------
TiledMetrics tiledAlignmentMetrics(const cv::Mat& im1,
                                   const cv::Mat& im2,
                                   const cv::Mat& mask,
                                   int    tile        = 50,
                                   int    nBins       = 32,
                                   double minCoverage = 0.5,
                                   double varFactor   = 4.0)
{
  // ---- 0. validate up front, with readable errors ------------------------
  // These fire as exceptions naming the actual types, which is what catches a
  // commented-out conversion or a caller handing you 16-bit microscopy data.
  CV_Assert(!im1.empty() && !im2.empty());
  CV_CheckEQ(im1.size(), im2.size(), "image pair must have identical size");
  CV_CheckDepthEQ(im1.depth(), CV_8U, "tiledAlignmentMetrics expects CV_8U");
  CV_CheckDepthEQ(im2.depth(), CV_8U, "tiledAlignmentMetrics expects CV_8U");
  CV_Assert(nBins > 0 && nBins <= 256 && (256 % nBins == 0));
  CV_Assert(tile > 1);
  if (!mask.empty()) {
    CV_CheckEQ(mask.size(), im1.size(), "mask must match image size");
    CV_CheckDepthEQ(mask.depth(), CV_8U, "mask must be CV_8U");
  }
  
  // ---- 1. reduce to single channel --------------------------------------
  // NOTE ON ALIASING: `g1 = im1` is a shallow copy -- the two Mats share one
  // pixel buffer. That is safe here only because we never write through g1.
  // Do not add in-place operations on g1/g2 without cloning first.
  cv::Mat g1, g2;
  if (im1.channels() > 1) cv::extractChannel(im1, g1, 0); else g1 = im1;
  if (im2.channels() > 1) cv::extractChannel(im2, g2, 0); else g2 = im2;
  
  // ---- 2. SSIM constants -------------------------------------------------
  // L is the dynamic range. For CV_8U it is 255 by definition, which is one
  // reason to keep the data 8-bit: you cannot get L wrong.
  const double L  = 255.0;
  const double C1 = (0.01 * L) * (0.01 * L);   //  6.5025
  const double C2 = (0.03 * L) * (0.03 * L);   // 58.5225
  
  // Flat-tile rejection threshold. Two INDEPENDENT near-constant patches score
  // SSIM ~= 1.0 and intersection ~= 1.0, because when the variances are small
  // next to C2 the constant dominates and forces the ratio toward 1. Empty
  // background tiles would then report perfect agreement while measuring
  // nothing, and in a tissue image most tiles are background. Requiring
  // max(varA,varB) >= varFactor*C2 (default 4*58.52 ~= 234) drops them.
  const double varThresh = varFactor * C2;
  
  // ---- 3. histogram setup, hoisted out of the loop ----------------------
  const int   histSize256 = 256;
  float       range[]     = {0.0f, 256.0f};   // NB upper bound is EXCLUSIVE
  const float* histRange  = {range};
  const int   channels[]  = {0};
  
  // ---- 4. tile grid ------------------------------------------------------
  const int nY = (im1.rows + tile - 1) / tile;   // ceil, so edge tiles included
  const int nX = (im1.cols + tile - 1) / tile;
  
  TiledMetrics R;
  R.n_total    = nY * nX;
  R.inter_map  = cv::Mat(nY, nX, CV_32F, cv::Scalar(std::nanf("")));
  R.bhatta_map = cv::Mat(nY, nX, CV_32F, cv::Scalar(std::nanf("")));
  R.ssim_map   = cv::Mat(nY, nX, CV_32F, cv::Scalar(std::nanf("")));
  
  // Weighted accumulators. We weight by valid pixel count rather than taking a
  // plain mean over tiles, because masked tiles contribute unequal areas and a
  // half-covered edge tile should not count as much as a full interior one.
  double wSum = 0.0, wInter = 0.0, wBhatta = 0.0, wSsim = 0.0;
  std::vector<double> ssimVals;
  ssimVals.reserve(R.n_total);
  
  cv::Mat h1, h2, prod;   // declared outside; OpenCV reuses the buffers
  
  // ---- 5. THE SINGLE LOOP -----------------------------------------------
  for (int ty = 0; ty < nY; ++ty) {
    for (int tx = 0; tx < nX; ++tx) {
      
      // Clip the ROI so edge tiles are simply smaller rather than skipped.
      const int y0 = ty * tile, x0 = tx * tile;
      const int th = std::min(tile, im1.rows - y0);
      const int tw = std::min(tile, im1.cols - x0);
      const cv::Rect roi(x0, y0, tw, th);
      
      // ROI views -- no pixel data is copied here.
      const cv::Mat ta = g1(roi);
      const cv::Mat tb = g2(roi);
      const cv::Mat tm = mask.empty() ? cv::Mat() : mask(roi);
      
      // --- GUARD 1: mask coverage. Cheap, so do it before any real work.
      const double area  = static_cast<double>(tw) * static_cast<double>(th);
      const double valid = tm.empty() ? area
      : static_cast<double>(cv::countNonZero(tm));
      if (valid < minCoverage * area) { ++R.n_drop_coverage; continue; }
      
      // --- ONE histogram per image. This is the whole trick: it yields the
      //     counts, the exact mean, the exact variance, and (after rebinning)
      //     the distributions for intersection and Bhattacharyya.
      cv::calcHist(&ta, 1, channels, tm, h1, 1, &histSize256, &histRange);
      cv::calcHist(&tb, 1, channels, tm, h2, 1, &histSize256, &histRange);
      
      double n1, mu1, var1, n2, mu2, var2;
      histStats(h1, n1, mu1, var1);
      histStats(h2, n2, mu2, var2);
      if (n1 <= 1.0 || n2 <= 1.0) { ++R.n_drop_coverage; continue; }
      
      // --- GUARD 2: flat tiles. Skip only when BOTH are near-constant. If one
      //     has texture and the other does not, that is a real mismatch and we
      //     want SSIM to report it as low, not to discard the tile.
      if (std::max(var1, var2) < varThresh) { ++R.n_drop_flat; continue; }
      
      // --- histogram metrics, off the rebinned pair (no extra image pass)
      const cv::Mat p1 = rebinAndNormalise(h1, nBins);
      const cv::Mat p2 = rebinAndNormalise(h2, nBins);
      const double inter  = cv::compareHist(p1, p2, cv::HISTCMP_INTERSECT);
      const double bhatta = cv::compareHist(p1, p2, cv::HISTCMP_BHATTACHARYYA);
      
      // --- the only thing histograms cannot give us: the cross term.
      //     dtype=CV_32F keeps the product out of 8-bit, where it would clamp
      //     at 255. cv::mean accumulates in double and honours the mask.
      cv::multiply(ta, tb, prod, 1.0, CV_32F);
      const double exy = cv::mean(prod, tm)[0];
      
      // cov = E[XY] - E[X]E[Y]. Both terms are doubles, so the subtraction is
      // done in double and the usual cancellation worry does not apply here.
      const double cov = exy - mu1 * mu2;
      
      // Standard SSIM with alpha=beta=gamma=1 and C3=C2/2, which collapses the
      // three components into two factors: luminance, then contrast*structure.
      const double lum = (2.0 * mu1 * mu2 + C1) / (mu1 * mu1 + mu2 * mu2 + C1);
      const double cs  = (2.0 * cov       + C2) / (var1      + var2      + C2);
      const double ssim = lum * cs;
      
      // --- record
      R.inter_map .at<float>(ty, tx) = static_cast<float>(inter);
      R.bhatta_map.at<float>(ty, tx) = static_cast<float>(bhatta);
      R.ssim_map  .at<float>(ty, tx) = static_cast<float>(ssim);
      
      const double w = valid;      // pixel-count weight
      wSum    += w;
      wInter  += w * inter;
      wBhatta += w * bhatta;
      wSsim   += w * ssim;
      ssimVals.push_back(ssim);
      ++R.n_valid;
    }
  }
  
  // ---- 6. aggregate ------------------------------------------------------
  if (wSum > 0.0) {
    R.intersection  = wInter  / wSum;
    R.bhattacharyya = wBhatta / wSum;
    R.ssim          = wSsim   / wSum;
    
    // Median is worth reporting alongside the mean: a handful of catastrophic
    // tiles (tissue tear, air bubble, stitching seam) drag the mean down while
    // the median still describes the bulk of the overlap.
    const size_t m = ssimVals.size() / 2;
    std::nth_element(ssimVals.begin(), ssimVals.begin() + m, ssimVals.end());
    R.ssim_median = ssimVals[m];
  }
  return R;
}

// ---------------------------------------------------------------------------
// Blow a coarse tile grid up to image resolution for overlay/visualisation.
// INTER_NEAREST keeps the tile blocks crisp -- do not interpolate, or you will
// invent gradients that were never measured.
// ---------------------------------------------------------------------------
cv::Mat upscaleMetricMap(const cv::Mat& tileMap, cv::Size dsize)
{
  cv::Mat out;
  cv::resize(tileMap, out, dsize, 0, 0, cv::INTER_NEAREST);
  return out;
}

// ---------------------------------------------------------------------------
// Hook for getAlignmentMetrics. Note the honest labelling: only the SSIM
// entries measure alignment. The two histogram entries measure whether the
// intensity distributions agree tile-by-tile, which is a different question --
// they are far less sensitive to misregistration than SSIM is, because a
// histogram discards pixel positions inside each tile.
// ---------------------------------------------------------------------------
void addTiledMetrics(std::map<std::string, double>& metrics,
                     const cv::Mat& im1, const cv::Mat& im2, const cv::Mat& mask,
                     int tile = 50)
{
  TiledMetrics t = tiledAlignmentMetrics(im1, im2, mask, tile);
  
  metrics["Tiled SSIM"]            = t.ssim;
  metrics["Tiled SSIM (median)"]   = t.ssim_median;
  metrics["Tiled Intersection"]    = t.intersection;
  metrics["Tiled Bhattacharyya"]   = t.bhattacharyya;
  metrics["Valid tiles"]           = static_cast<double>(t.n_valid);
  
  Rcout << "Tiled metrics (" << tile << "px tiles): " << std::endl;
  Rcout << "  SSIM:          " << t.ssim
        << "  (median " << t.ssim_median << ")" << std::endl;
  Rcout << "  Intersection:  " << t.intersection << std::endl;
  Rcout << "  Bhattacharyya: " << t.bhattacharyya << std::endl;
  Rcout << "  Tiles:         " << t.n_valid << "/" << t.n_total
        << " valid  (dropped " << t.n_drop_coverage << " low-coverage, "
        << t.n_drop_flat << " flat)" << std::endl;
  
  // A score computed from a handful of tiles is not comparable to one computed
  // from hundreds. Say so rather than letting the number stand alone.
  if (t.n_valid < 10) {
    Rcout << "  WARNING: very few valid tiles - metrics are unreliable."
          << std::endl;
  }
}

// ---------------------------------------------------------------------------
// Rcpp entry point. Returns the scalars plus the three tile maps so you can
// plot them in R (image(), or overlay after upscaling).
// Assumes you already have imageToMat / matToImage in auxiliary.h+image.h.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List tiled_alignment_metrics(Rcpp::NumericMatrix& ref_image,
                                   Rcpp::NumericMatrix& query_image,
                                   Rcpp::NumericMatrix& mask_mat,
                                   int tile = 50,
                                   int bins = 32,
                                   double min_coverage = 0.5,
                                   double var_factor = 4.0)
{
  cv::Mat im1  = numericMatrixToMat(ref_image);
  cv::Mat im2  = numericMatrixToMat(query_image);
  cv::Mat mask = numericMatrixToMat(mask_mat);
  
  // numericMatrixToMat gives doubles; the metrics want 8-bit. This is the one
  // permitted conversion: it is a type *fix* at the R boundary, not an
  // avoidable widening inside the hot loop.
  if (im1.depth()  != CV_8U) im1.convertTo(im1, CV_8U);
  if (im2.depth()  != CV_8U) im2.convertTo(im2, CV_8U);
  if (!mask.empty() && mask.depth() != CV_8U) mask.convertTo(mask, CV_8U);
  
  TiledMetrics t = tiledAlignmentMetrics(im1, im2, mask, tile, bins,
                                         min_coverage, var_factor);
  
  return Rcpp::List::create(
    Rcpp::Named("ssim")            = t.ssim,
    Rcpp::Named("ssim_median")     = t.ssim_median,
    Rcpp::Named("intersection")    = t.intersection,
    Rcpp::Named("bhattacharyya")   = t.bhattacharyya,
    Rcpp::Named("n_valid")         = t.n_valid,
    Rcpp::Named("n_total")         = t.n_total,
    Rcpp::Named("n_drop_coverage") = t.n_drop_coverage,
    Rcpp::Named("n_drop_flat")     = t.n_drop_flat,
    Rcpp::Named("ssim_map")        = matToImage(t.ssim_map),
    Rcpp::Named("inter_map")       = matToImage(t.inter_map),
    Rcpp::Named("bhatta_map")      = matToImage(t.bhatta_map)
  );
}
