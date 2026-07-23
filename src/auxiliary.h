#include "Rcpp.h"
#include <opencv2/opencv.hpp>
#include "opencv2/shape/shape_transformer.hpp"

#ifndef AUXILIARY_H
#define AUXILIARY_H

////
// replace NAs in matrices
////

Rcpp::NumericMatrix replaceNaMatrix(Rcpp::NumericMatrix mat, int replace);

////
// Conversion
////

// cv::Mat vs Rcpp::RawVector(Image)
Rcpp::RawVector matToImage(const cv::Mat &mat);
cv::Mat imageToMat(Rcpp::RawVector &image_data, int width, int height);

// cv::Mat vs Rcpp::RawVector(Image) with 2 dim (mostly for masks)
Rcpp::IntegerVector matToMask(const cv::Mat &mat);

// cv::Mat vs Rcpp::NumericMatrix 
cv::Mat numericMatrixToMat(Rcpp::NumericMatrix nm);
Rcpp::NumericMatrix matToNumericMatrix(cv::Mat m);

// std::vector<cv::Point2f> vs Rcpp::NumericMatrix
std::vector<cv::Point2f> numericMatrixToPoint2f(Rcpp::NumericMatrix mat);
Rcpp::NumericMatrix point2fToNumericMatrix(std::vector<cv::Point2f> &points);
  
// std::vector<cv::Point2f> vs cv::Mat
std::vector<cv::Point2f> matToPoint2f(cv::Mat &mat);
cv::Mat point2fToMat(std::vector<cv::Point2f> &points);

// std::vector<uint8_t> vs cv::Mat
cv::Mat IntVectorToMat(std::vector<uint8_t> &points);

// std::vector<double> vs std::vector<cv::Point2f>
std::vector<double> KeyPointToDoubleVector(std::vector<cv::KeyPoint> &points);
std::vector<double> Point2fToDoubleVector(std::vector<cv::Point2f> &points);

// std::vector<cv::KeyPoint> vs std::vector<cv::Point2f>
std::vector<cv::Point2f> KeyPointToPoint2f(std::vector<cv::KeyPoint> &keypoints);
  
////
// stats
////
  
// standard deviation
double cppSD(std::vector<cv::KeyPoint> &points);
double cppSD(std::vector<cv::Point2f> &points);

// mean distance between points
double meanDistances(std::vector<cv::Point2f> &pts1, std::vector<cv::Point2f> &pts2);
double medianDistances(std::vector<cv::Point2f> &pts1, std::vector<cv::Point2f> &pts2);

#endif