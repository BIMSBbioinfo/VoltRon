#include <Rcpp.h>

// OpenCV
#include <opencv2/opencv.hpp>
#include "opencv2/features2d.hpp"
#include "opencv2/shape/shape_transformer.hpp"
// #include <opencv2/imgproc.hpp>

// Internal functions
#include "auxiliary.h"
#include "image.h"

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
  double gridpoints_reg_sd = cppSD(gridpoints_reg);
  
  // get warning message
  if(gridpoints_reg_sd < 1.0 | gridpoints_reg_sd > max(height, width)){
    Rcout << "  WARNING: Transformation may be poor - transformed points grid seem to be concentrated!" << endl;
  }
  
  return gridpoints_reg_sd;
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

// void maskKeypoints(std::vector<cv::KeyPoint> &keypoints1_good, std::vector<cv::KeyPoint> &keypoints2_good,
//                     std::vector<cv::KeyPoint> &keypoints1_masked, std::vector<cv::KeyPoint> &keypoints2_masked, 
//                     Mat &mask)
// {
//   int j=0;
//   for (int i = 0; i < mask.rows; i++) {
//     if (mask.at<uchar>(i)) {
//       keypoints1_masked.push_back(keypoints1_good[i]);
//       keypoints2_masked.push_back(keypoints2_good[i]);
//       j++;
//     }
//   }
// }

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

cv::Mat generateOverlapMask(cv::Mat& im, cv::Mat& h, 
                            cv::Size dsize, cv::Size ssize)
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

std::vector<double> getAlignmentMetrics(Mat &im1, Mat &im2, Mat &h, 
                                        cv::Size ssize){
  
  // Histogram settings
  int histSize = 256;
  float range[] = {0.0, 256.0};
  const float* histRange = {range};
  int channels[] = {0};
  
  // get overlap mask
  cv::Mat alignmentMask = generateOverlapMask(im1, h, im2.size(), ssize);
  // imwrite("mask.tif", alignmentMask);
  
  // Compute histograms
  cv::Mat hist1, hist2;
  cv::calcHist(&im1, 1, channels, alignmentMask, 
               hist1, 1, &histSize, &histRange);
  cv::calcHist(&im2, 1, channels, alignmentMask, 
               hist2, 1, &histSize, &histRange);
  
  // Normalize histograms
  cv::normalize(hist1, hist1, 0, 1, cv::NORM_MINMAX);
  cv::normalize(hist2, hist2, 0, 1, cv::NORM_MINMAX);
  
  // Summary
  Rcout << "Alignment Report: " << endl;
  std::vector<double> metrics;
  metrics.push_back(cv::compareHist(hist1, hist2, cv::HISTCMP_INTERSECT));
  metrics.push_back(cv::compareHist(hist1, hist2, cv::HISTCMP_BHATTACHARYYA));
  //metrics.push_back(cv::compareHist(hist1, hist2, cv::HISTCMP_CHISQR));
  // metrics.push_back(jointEntropy(im1, im2, alignmentMask, histSize));
  // metrics.push_back(MutualInfo(im1, im2, alignmentMask, histSize));
  // metrics.push_back(NormalizedMutualInfo(im1, im2, alignmentMask, histSize));
  
  Rcout << "  Intersection:     " << metrics[0] << std::endl;
  Rcout << "  Bhattacharyya:    " << metrics[1] << std::endl;
  // Rcout << "  Chi-Square:       " << metrics[0] << std::endl;
  // Rcout << "  Joint Entropy:    " << metrics[3] << std::endl;
  // Rcout << "  MutualInfo:       " << metrics[4] << std::endl;
  // Rcout << "  NormalizedMutualInfo: " << metrics[5] << std::endl;
  return metrics;
}

// do overall checks on keypoints and images
std::vector<double> getKeypointMetrics(std::vector<cv::Point2f> &points1, 
                                       std::vector<cv::Point2f> &points2, 
                                       Mat &im1, Mat &im2, 
                                       Mat &h, Mat &mask) {
  
  // metrics list
  std::vector<double> metrics_list;
  
  // Alignment report
  Rcout << "Keypoint Report: " << endl;
  
  // Report final keypoints
  Rcout << "  Calculated transformation matrix with " << points1.size() << " keypoints" << endl;
  
  // points stand. dev.
  double points1_sd = cppSD(points1);
  double points2_sd = cppSD(points2);
  if(points1_sd < 1.0 | points2_sd < 1.0){
    Rcout << "  WARNING: points may be in a degenerate configuration." << endl;
  } 
  Rcout << "  Std dev of points: x=" << points1_sd << " y="  << points2_sd << endl;
  metrics_list.push_back(checkDegenerate(points1_sd, points2_sd));
  metrics_list.push_back(points1_sd);
  metrics_list.push_back(points2_sd);
  
  // check distribution of points
  double stddev = checkMappedGridDistribution(im2, h);
  Rcout << "  Std dev of registered points: " << stddev << endl;
  metrics_list.push_back(stddev);
  
  // warp keypoints and compare 
  double md = medianMappingDistance(points1, points2, h);
  Rcout << "  Median distance between points: " << md << endl;
  if(md > 3){
    Rcout << "  WARNING: Transformation may be poor - mean euclidean distance of mapped source and destination key points is high!" << endl;
  } 
  
  // get inlier percentages
  double ratio = checkInlierPercentage(mask);
  Rcout << "  Inlier Percentage: " << ratio << endl;
  
  // degenerate ?
  Rcout << "Registration is " << (metrics_list[0] ? "degenerate!" : "not degenerate!") << endl;
  
  // return is_degenerate;
  return metrics_list;
}