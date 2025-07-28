#include "Rcpp.h"
#include <opencv2/opencv.hpp>
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/features2d.hpp"
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
Rcpp::RawVector matToImage(cv::Mat mat);
cv::Mat imageToMat(Rcpp::RawVector &image_data, int width, int height);

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
  
////
// validation
////
  
// standard deviation
double cppSD(std::vector<cv::KeyPoint> &points);
double cppSD(std::vector<cv::Point2f> &points);
  
#endif