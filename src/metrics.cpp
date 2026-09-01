#include <Rcpp.h>

// OpenCV
#include <opencv2/opencv.hpp>
#include "opencv2/features2d.hpp"
#include "opencv2/shape/shape_transformer.hpp"

// Internal functions
#include "auxiliary.h"
#include "image.h"
#include "matte_mi.h"
#include <algorithm>
#include <vector>
#include <cmath>

// Namespaces
using namespace Rcpp;
using namespace std;
using namespace cv;

////
// Quality Control
////

// check distribution of registered points
std::vector<double> checkMappedGridDistribution(Mat &im, Mat &h){
  
  // message
  std::string message;
  
  // get image shape
  int height = im.rows;
  int width = im.cols;
  int height_interval = height > 50 ? (double) height/50.0 : 1;
  int width_interval = width > 50 ? (double) width/50.0 : 1;

  // perspective transformation of grid points
  std::vector<cv::Point2f> gridpoints;
  for (double i = 0.0; i <= height; i += height_interval) {
    for (double j = 0.0; j <= width; j += width_interval) {
      gridpoints.push_back(cv::Point2f(j,i));
    }
  }

  // register grid points
  std::vector<cv::Point2f> gridpoints_reg;
  if (h.rows == 2){
    cv::transform(gridpoints, gridpoints_reg, h);
  } else if(h.rows == 3) {
    cv::perspectiveTransform(gridpoints, gridpoints_reg, h);
  } 

  // Compute the standard deviation of the transformed points
  std::vector<double> stds(2);
  stds[0] = cppSD(gridpoints_reg, 0);
  stds[1] = cppSD(gridpoints_reg, 1);
  return(stds);
}

bool checkMaskAbundance(Mat &mask){
  int j=0;
  for (int i = 0; i < mask.rows; i++) {
    if (mask.at<uchar>(i)) {
      j++;
    }
  }
  return j > 6;
}

// compare the distance between two sets of match points
double medianMappingDistance(std::vector<cv::Point2f> &keypoints1, std::vector<cv::Point2f> &keypoints2, Mat &h) {
  std::vector<cv::Point2f> keypoints1_warped;
  if(keypoints1.size() > 0){
    if (h.rows == 2){
      cv::transform(keypoints1, keypoints1_warped, h);
    } else {
      cv::perspectiveTransform(keypoints1, keypoints1_warped, h);
    }
  }
  
  return medianDistances(keypoints1_warped, keypoints2);
}

// calculate inlier percentage
int checkInlierPercentage(Mat &mask){
  int j=0;
  for (int i = 0; i < mask.rows; i++) {
    if (mask.at<uchar>(i)) {
      j++;
    }
  }
  double ratio = (double) j/mask.rows; 
  double perc = round(100.0 * ratio);
  return (int) perc;
}

void maskKeypoints(std::vector<cv::KeyPoint> &keypoints1_good, std::vector<cv::KeyPoint> &keypoints2_good,
                   std::vector<cv::KeyPoint> &keypoints1_masked, std::vector<cv::KeyPoint> &keypoints2_masked, 
                   std::vector<cv::DMatch> &top_matches, Mat &mask)
{
  int j=0;
  for (int i = 0; i < mask.rows; i++) {
    if (mask.at<uchar>(i)) {
      keypoints1_masked.push_back(keypoints1_good[i]);
      keypoints2_masked.push_back(keypoints2_good[i]);
      top_matches.push_back(cv::DMatch(static_cast<int>(j), static_cast<int>(j), 0));
      j++;
    }
  }
}

// check if keypoints are degenerate
bool checkDegenerate(double pts1, double pts2, double pts3, double pts4) {
  
  // get warning message
  bool is_degenerate_x = FALSE;
  if(pts1 < 1.0 || pts2 < 1.0){
    is_degenerate_x = TRUE;
  } 
  bool is_degenerate_y = FALSE;
  if(pts3 < 1.0 || pts4 < 1.0){
    is_degenerate_y = TRUE;
  } 
  if(is_degenerate_y || is_degenerate_x){
    Rcout << "WARNING: points may be in a degenerate configuration." << endl;
  }

  return (is_degenerate_y || is_degenerate_x);
}

cv::Mat generateOverlapMask(cv::Size dsize,
                            cv::Mat& h, 
                            cv::Size ssize)
{
  // generate mask
  cv::Mat mask = cv::Mat::ones(ssize, CV_8UC1) * 255;
  cv::Mat warped;
  
  // Keep masks crisp: nearest-neighbor only.
  const int interp = cv::INTER_NEAREST;
  const int borderMode = cv::BORDER_CONSTANT;
  const cv::Scalar borderValue(0);

  // warp mask
  if (h.rows == 2){
    cv::warpAffine(mask, warped, h, dsize, 
                   interp, borderMode, borderValue);
  } else {
    cv::warpPerspective(mask, warped, h, dsize, 
                        interp, borderMode, borderValue);
  }

  // Force binary mask again.
  cv::threshold(warped, warped, 0, 255, cv::THRESH_BINARY);
  return warped;
}

cv::Mat generateOverlapMask(cv::Mat& ref_image, 
                            Ptr<ThinPlateSplineShapeTransformer>& tps,
                            cv::Size ssize)
{
  // generate mask
  cv::Mat mask = cv::Mat::ones(ssize, CV_8UC1) * 255;

  // Keep masks crisp: nearest-neighbor only.
  const int interp = cv::INTER_NEAREST;
  
  mask = warpTPSImage(ref_image, mask, tps, 
                      ref_image.rows, ref_image.cols, interp);
  
  // Force binary mask again.
  cv::threshold(mask, mask, 0, 255, cv::THRESH_BINARY);
  return mask;
}

// [[Rcpp::export]]
Rcpp::IntegerVector generateOverlapMask(Rcpp::NumericVector& dsize, 
                                        Rcpp::NumericMatrix& trans_mat,
                                        Rcpp::NumericVector& ssize){
  cv::Mat h = numericMatrixToMat(trans_mat);
  cv::Mat mask = generateOverlapMask(cv::Size((int) dsize[0], (int) dsize[1]),
                                     h,
                                     cv::Size((int) ssize[0], (int) ssize[1]));
  return matToMask(mask);
  // return matToImage(mask);
}

double Entropy(cv::Mat& im1, cv::Mat& overlapMask, int bins = 256) {
  
  // Histogram settings
  int histSize = 256;
  float range[] = {0.0, 256.0};
  const float* histRange = {range};
  int channels[] = {0};
  
  // Compute histograms
  cv::Mat hist;
  cv::calcHist(&im1, 1, channels, overlapMask, 
               hist, 1, &histSize, &histRange);
  
  // Normalize histograms
  cv::normalize(hist, hist, 0, 1, cv::NORM_MINMAX);
  
  // Convert counts to probabilities
  hist /= cv::sum(hist)[0];
  
  double entropy = 0.0;
  for (int r = 0; r < hist.rows; ++r)
  {
    const float* ptr = hist.ptr<float>(r);
    
    for (int c = 0; c < hist.cols; ++c)
    {
      double p = ptr[c];
      
      if (p > 0.0)
        entropy -= p * std::log(p);
    }
  }
  
  return entropy;
}

double jointEntropy(cv::Mat& im1, cv::Mat& im2, 
                    cv::Mat& overlapMask, int bins = 256) {
  
  // 2D histogram parameters
  int histSize[] = {bins, bins};
  float range[] = {0.f, 256.f};
  const float* ranges[] = {range, range};
  int channels[] = {0, 1};
  
  // calculate histogram
  cv::Mat images[] = {im1, im2};
  cv::Mat hist;
  cv::calcHist(images,
               2,
               channels,
               overlapMask,
               hist,
               2,
               histSize,
               ranges,
               true,
               false);
  cv::normalize(hist, hist, 0, 1, cv::NORM_MINMAX);
  
  // Convert counts to probabilities
  hist /= cv::sum(hist)[0];
  
  double entropy = 0.0;
  for (int r = 0; r < hist.rows; ++r)
  {
    const float* ptr = hist.ptr<float>(r);
    
    for (int c = 0; c < hist.cols; ++c)
    {
      double p = ptr[c];
      
      if (p > 0.0)
        entropy -= p * std::log(p);
    }
  }
  
  return entropy;
}

double MutualInfo(cv::Mat& im1, cv::Mat& im2, 
                  cv::Mat& overlapMask, int bins = 256) {
  double ent1=Entropy(im1, overlapMask, bins);
  double ent2=Entropy(im2, overlapMask, bins);
  double ent12=jointEntropy(im1, im2, overlapMask, bins);
  return ent1+ent2-ent12;
}

double NormalizedMutualInfo(cv::Mat& im1, cv::Mat& im2, 
                  cv::Mat& overlapMask, int bins = 256) {
  double ent1=Entropy(im1, overlapMask, bins);
  double ent2=Entropy(im2, overlapMask, bins);
  double ent12=jointEntropy(im1, im2, overlapMask, bins);
  return (ent1+ent2)/ent12;
}

const double kNaN = std::numeric_limits<double>::quiet_NaN();

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
  
  cv::Mat1d inter_map;           // CV_64FC1, nTilesY x nTilesX, NaN = skipped
  cv::Mat1d bhatta_map;
  cv::Mat1d ssim_map;
};

// ---------------------------------------------------------------------------
// Exact n / mean / variance from a 256-bin histogram of an 8-bit tile.
// Bin i corresponds to intensity value i exactly, so nothing is approximated.
// The variance convention matches cv::meanStdDev (population, divide by n).
// ---------------------------------------------------------------------------
void histStats(const cv::Mat& h256, double& n, double& mu, double& var) {
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
cv::Mat rebinAndNormalise(const cv::Mat& h256, int nBins) {
  const int group = 256 / nBins;
  cv::Mat out = cv::Mat::zeros(nBins, 1, CV_32F);
  const float* p = h256.ptr<float>(0);
  float*       q = out.ptr<float>(0);
  for (int i = 0; i < 256; ++i) q[i / group] += p[i];
  
  const double s = cv::sum(out)[0];
  if (s > 0.0) out /= s;
  return out;
}

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
  R.inter_map  = cv::Mat1d(nY, nX, kNaN);
  R.bhatta_map = cv::Mat1d(nY, nX, kNaN);
  R.ssim_map   = cv::Mat1d(nY, nX, kNaN);
  
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
      // if (std::max(var1, var2) < varThresh) { ++R.n_drop_flat; continue; }
      
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
      // R.inter_map .at<float>(ty, tx) = static_cast<float>(inter);
      // R.bhatta_map.at<float>(ty, tx) = static_cast<float>(bhatta);
      // R.ssim_map.at<float>(ty, tx) = static_cast<float>(ssim);
      // --- record
      R.inter_map (ty, tx) = inter;    // Mat1d: no template arg, no cast
      R.bhatta_map(ty, tx) = bhatta;
      R.ssim_map  (ty, tx) = ssim;
      
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

cv::Mat1d getTiledAlignmentMetrics(Mat &im1, Mat &im2, Mat &mask){
  TiledMetrics t = tiledAlignmentMetrics(im1, im2, mask);
  return t.ssim_map;
}

std::map<std::string, double> getAlignmentMetrics(Mat &im1, Mat &im2, 
                                                  Mat &mask, std::string type){
  
  // Metrics
  std::map<std::string, double> metrics;
  
  // // Compute histograms
  // int histSize = 256;
  // float range[] = {0.0, 256.0};
  // const float* histRange = {range};
  // int channels[] = {0};
  // cv::Mat hist1, hist2;
  // cv::calcHist(&im1, 1, channels, mask, 
  //              hist1, 1, &histSize, &histRange);
  // cv::calcHist(&im2, 1, channels, mask, 
  //              hist2, 1, &histSize, &histRange);
  
  // // Normalize histograms
  // // cv::normalize(hist1, hist1, 0, 1, cv::NORM_MINMAX);
  // // cv::normalize(hist2, hist2, 0, 1, cv::NORM_MINMAX);
  // hist1 /= cv::sum(hist1)[0];
  // hist2 /= cv::sum(hist2)[0];
  
  // tiled metrics:
  TiledMetrics t = tiledAlignmentMetrics(im1, im2, mask);
  
  // Summary
  Rcout << "Alignment Accuracy (" << type << "): " << endl;
  // metrics["Intersection"] = cv::compareHist(hist1, hist2, cv::HISTCMP_INTERSECT);
  // metrics["Bhattacharyya"] = cv::compareHist(hist1, hist2, cv::HISTCMP_BHATTACHARYYA);
  // metrics["SSIM"] = tiledSSIM(im2, im1, mask, 50);
  metrics["Matte's MI"] = MatteMI(im2, im1, mask, 50);
  metrics["SSIM"] = t.ssim;
  metrics["Intersection"] = t.intersection;
  metrics["Bhattacharyya"] = t.bhattacharyya;
  
  Rcout << "  Matte's MI:    " << metrics["Matte's MI"] << std::endl;
  Rcout << "  SSIM:          " << metrics["SSIM"] << std::endl;
  Rcout << "  Intersection:  " << metrics["Intersection"] << std::endl;
  Rcout << "  Bhattacharyya: " << metrics["Bhattacharyya"] << std::endl;
  
  return metrics;
}

// do overall checks on keypoints and images
std::map<std::string, double> getKeypointMetrics(std::vector<cv::Point2f> &points1, 
                                       std::vector<cv::Point2f> &points2, 
                                       Mat &im1, Mat &im2, 
                                       Mat &h, Mat &mask) {
  
  // metrics list
  std::map<std::string, double> metrics;
  
  // Alignment report
  Rcout << "Keypoint Report: " << endl;
  
  // Report final keypoints
  Rcout << "  Calculated transformation matrix with " << points1.size() << " keypoints" << endl;
  metrics["#Matches"] = points1.size();
  
  // get inlier percentages
  double ratio = checkInlierPercentage(mask);
  Rcout << "  Inlier Percentage: " << ratio << endl;
  metrics["Inlier Perc."] = ratio;
  
  // points stand. dev.
  double points1_sd = cppSD(points1);
  double points1_sd_y = cppSD(points1, 1);
  metrics["sd query(x) (>1?)"] = points1_sd;
  metrics["sd query(y) (>1?)"] = points1_sd_y;
  Rcout << "  Std dev of query points: x=" << points1_sd << " y="  << points1_sd_y << endl;
  double points2_sd = cppSD(points2);
  double points2_sd_y = cppSD(points2, 1);
  metrics["sd ref(x) (>1?)"] = points2_sd;
  metrics["sd ref(y) (>1?)"] = points2_sd_y;
  Rcout << "  Std dev of ref points: x=" << points2_sd << " y="  << points2_sd_y << endl;
  
  // degenerate ?
  bool degenerate_points = checkDegenerate(points1_sd, points2_sd, 
                                           points1_sd_y, points2_sd_y);
  metrics["Degenerate"] = (double) degenerate_points;
  
  // check distribution of points
  std::vector<double> stddev = checkMappedGridDistribution(im1, h);
  Rcout << "  Std dev of registered grid points: x=" << stddev[0] << " y="  << stddev[1] << endl;
  if(stddev[0] < 1.0 || stddev[0] > im2.cols ||
     stddev[1] < 1.0 || stddev[1] > im2.rows ){
    Rcout << "  WARNING: Transformation may be poor - transformed points grid seem to be concentrated!" << endl;
    metrics["Degenerate"] = 1.0;
  }
  metrics["sd grid (x) [w,h]?"] = stddev[0];
  metrics["sd grid (y) [w,h]?"] = stddev[1];
  
  // warp keypoints and check median distances 
  double md = medianMappingDistance(points1, points2, h);
  Rcout << "  Median distance between points: " << md << endl;
  metrics["Median distance"] = md;

  // report degenerate
  if((bool) metrics["Degenerate"]){
    Rcout << "  WARNING: Registration is degenerate!" << endl;
  } 
  
  // return is_degenerate;
  return metrics;
}