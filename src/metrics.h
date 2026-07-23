#include "Rcpp.h"
#include <opencv2/opencv.hpp>
#include "opencv2/shape/shape_transformer.hpp"

// Namespaces
using namespace Rcpp;
using namespace std;
using namespace cv;

#ifndef METRICS_H
#define METRICS_H

// check distribution of registered points
double checkMappedGridDistribution(cv::Mat &im, cv::Mat &h);

bool checkMaskAbundance(cv::Mat &mask);

// compare the distance between two sets of match points
double medianMappingDistance(std::vector<cv::Point2f> &keypoints1, std::vector<cv::Point2f> &keypoints2, cv::Mat &h);

// calculate inlier percentage
int checkInlierPercentage(cv::Mat &mask);

void maskKeypoints(std::vector<cv::KeyPoint> &keypoints1_good, std::vector<cv::KeyPoint> &keypoints2_good,
                   std::vector<cv::KeyPoint> &keypoints1_masked, std::vector<cv::KeyPoint> &keypoints2_masked, 
                   std::vector<cv::DMatch> &top_matches, cv::Mat &mask);

// check if keypoints are degenerate
bool checkDegenerate(double pts1, double pts2);

// generate overlap mask for alignment
cv::Mat generateOverlapMask(cv::Size dsize,
                            cv::Mat& h, 
                            cv::Size ssize);

cv::Mat generateOverlapMask(cv::Mat& ref_image, 
                            Ptr<ThinPlateSplineShapeTransformer>& tps, 
                            cv::Size ssize);

// cv::Mat generateOverlapMask(Rcpp::NumericVector dsize,
//                             Rcpp::NumericMatrix trans_mat,
//                             Rcpp::NumericVector ssize);

// get alignment metrics
std::map<std::string, double> getAlignmentMetrics(cv::Mat &im1, 
                                                  cv::Mat &im2, 
                                                  cv::Mat &mask, 
                                                  std::string type);

// do overall checks on keypoints and metrics
std::map<std::string, double> getKeypointMetrics(std::vector<cv::Point2f> &points1, 
                                       std::vector<cv::Point2f> &points2, 
                                       cv::Mat &im1, cv::Mat &im2, 
                                       cv::Mat &h, cv::Mat &mask);
  
#endif