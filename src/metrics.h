#include "Rcpp.h"
#include <opencv2/opencv.hpp>

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

cv::Mat generateOverlapMask(cv::Mat& im, cv::Mat& h, cv::Size dsize);

std::vector<double> getAlignmentMetrics(cv::Mat &im1, cv::Mat &im2, cv::Mat &h);

// do overall checks on keypoints and images
std::vector<double> getKeypointMetrics(std::vector<cv::Point2f> &points1, 
                                       std::vector<cv::Point2f> &points2, 
                                       cv::Mat &im1, cv::Mat &im2, 
                                       cv::Mat &h, cv::Mat &mask);
  
#endif