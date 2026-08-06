#include <Rcpp.h>

// OpenCV
#include <opencv2/opencv.hpp>
#include "opencv2/features2d.hpp"
#include "opencv2/shape/shape_transformer.hpp"

// Internal functions
#include "auxiliary.h"
#include "image.h"
#include "matte_mi.h"

// Namespaces
using namespace Rcpp;
using namespace std;
using namespace cv;

////
// Quality Control
////

// check distribution of registered points
double checkMappedGridDistribution(Mat &im, Mat &h){
  
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
  return cppSD(gridpoints_reg);
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
bool checkDegenerate(double pts1, double pts2) {
  
  // get warning message
  bool is_degenerate = FALSE;
  if(pts1 < 1.0 | pts2 < 1.0){
    is_degenerate = TRUE;
    Rcout << "WARNING: points may be in a degenerate configuration." << endl;
  } 
  
  return is_degenerate;
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
  // const int borderMode = cv::BORDER_CONSTANT;
  // const cv::Scalar borderValue(0);
  
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

std::map<std::string, double> getAlignmentMetrics(Mat &im1, Mat &im2, 
                                                  Mat &mask, std::string type){
  
  // Metrics
  std::map<std::string, double> metrics;
  
  // Compute histograms
  int histSize = 256;
  float range[] = {0.0, 256.0};
  const float* histRange = {range};
  int channels[] = {0};
  cv::Mat hist1, hist2;
  cv::calcHist(&im1, 1, channels, mask, 
               hist1, 1, &histSize, &histRange);
  cv::calcHist(&im2, 1, channels, mask, 
               hist2, 1, &histSize, &histRange);
  
  // Normalize histograms
  cv::normalize(hist1, hist1, 0, 1, cv::NORM_MINMAX);
  cv::normalize(hist2, hist2, 0, 1, cv::NORM_MINMAX);
  hist1 /= cv::sum(hist1)[0];
  hist2 /= cv::sum(hist2)[0];
  
  // Summary
  Rcout << "Alignment Accuracy (" << type << "): " << endl;
  metrics["Intersection"] = cv::compareHist(hist1, hist2, cv::HISTCMP_INTERSECT);
  metrics["Bhattacharyya"] = cv::compareHist(hist1, hist2, cv::HISTCMP_BHATTACHARYYA);
  metrics["Matte's MI"] = MatteMI(im2, im1, mask, 50);
  
  Rcout << "  Intersection:  " << metrics["Intersection"] << std::endl;
  Rcout << "  Bhattacharyya: " << metrics["Bhattacharyya"] << std::endl;
  Rcout << "  Matte's MI:    " << metrics["Matte's MI"] << std::endl;
  
  // old metrics, keep for comparison
  //metrics.push_back(cv::compareHist(hist1, hist2, cv::HISTCMP_CHISQR));
  // metrics.push_back(jointEntropy(im1, im2, mask, histSize));
  // metrics.push_back(MutualInfo(im1, im2, mask, histSize));
  // metrics.push_back(NormalizedMutualInfo(im1, im2, mask, histSize));
  
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
  metrics["#Keypoints"] = points1.size();
  
  // get inlier percentages
  double ratio = checkInlierPercentage(mask);
  Rcout << "  Inlier Percentage: " << ratio << endl;
  metrics["Inlier Perc."] = ratio;
  
  // points stand. dev.
  double points1_sd = cppSD(points1);
  double points2_sd = cppSD(points2);
  Rcout << "  Std dev of points: x=" << points1_sd << " y="  << points2_sd << endl;
  metrics["sd query kpts (>1?)"] = points1_sd;
  metrics["sd ref. kpts (>1?)"] = points2_sd;
  
  // degenerate ?
  bool degenerate_points = checkDegenerate(points1_sd, points2_sd);
  metrics["Degenerate"] = (double) degenerate_points;
  
  // check distribution of points
  double stddev = checkMappedGridDistribution(im1, h);
  Rcout << "  Std dev of registered points: " << stddev << endl;
  if(stddev < 1.0 | stddev > max(im2.rows, im2.cols)){
    Rcout << "  WARNING: Transformation may be poor - transformed points grid seem to be concentrated!" << endl;
    metrics["Degenerate"] = 1.0;
  }
  metrics["sd grid (in [w,h]?)"] = stddev;
  
  // warp keypoints and check median distances 
  double md = medianMappingDistance(points1, points2, h);
  Rcout << "  Median distance between points: " << md << endl;
  if(md > 3){
    Rcout << "  WARNING: Transformation may be poor - mean euclidean distance of mapped source and destination key points is high!" << endl;
  } 
  metrics["Median distance"] = md;

  // report degenerate
  if((bool) metrics["Degenerate"]){
    Rcout << "  WARNING: Registration is degenerate!" << endl;
  } 
  
  // return is_degenerate;
  return metrics;
}