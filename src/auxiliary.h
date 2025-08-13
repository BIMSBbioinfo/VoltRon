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
Rcpp::RawVector matToImage(const cv::Mat &mat);
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
// stats
////
  
// standard deviation
double cppSD(std::vector<cv::KeyPoint> &points);
double cppSD(std::vector<cv::Point2f> &points);

////
// memory
////

// void log_mem_usage(const std::string& label);
// void log_mem_macos(const std::string& label);
// double object_size_long(long bsize);
// double object_size_double(double bsize);
// double get_resident_bytes();
// double bytes_to_gb(double bytes);
// 
// struct MemProfiler {
//   size_t start;
//   std::string label;
// 
//   MemProfiler(const std::string& lbl) : label(lbl) {
//     start = get_resident_bytes();
//   }
// 
//   ~MemProfiler() {
//     size_t end = get_resident_bytes();
//     double diff = (double) end - (double) start;
//     if(diff < 0.0) diff = 0.0;
//     double diff_gb = bytes_to_gb(diff);
//     Rcpp::Rcout << "[MEM] " << label << ": +" << diff_gb << " GB" << std::endl;
//   }
// };

#endif