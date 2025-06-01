#include "Rcpp.h"
#include <opencv2/opencv.hpp>
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/shape/shape_transformer.hpp"

#ifndef IMAGE_H
#define IMAGE_H

////
// Processing
////

cv::Mat preprocessImage(cv::Mat &im, const bool invert, const char* flipflop, const char* rotate);
cv::Mat reversepreprocessImage(cv::Mat &im, const char* flipflop, const char* rotate);
cv::Mat resize_image(cv::Mat im, int width);
std::vector<cv::KeyPoint> resize_keypoints(std::vector<cv::KeyPoint> keypoints, 
                                           cv::Mat im,
                                           int width);
void scaledDrawMatches(cv::Mat im1, std::vector<cv::KeyPoint> keypoints1, 
                       cv::Mat im2, std::vector<cv::KeyPoint> keypoints2,
                       std::vector<cv::DMatch> top_matches,
                       cv::Mat &imMatches);

////
// Warping
////
  
Rcpp::RawVector warpImage(Rcpp::RawVector ref_image, Rcpp::RawVector query_image, 
                          Rcpp::List mapping,
                          const int width1, const int height1,
                          const int width2, const int height2);

Rcpp::RawVector warpImageAuto(Rcpp::RawVector ref_image, Rcpp::RawVector query_image, 
                              Rcpp::List mapping,
                              const int width1, const int height1,
                              const int width2, const int height2);

Rcpp::RawVector warpImageManual(Rcpp::RawVector ref_image, Rcpp::RawVector query_image, 
                                Rcpp::List mapping,
                                const int width1, const int height1,
                                const int width2, const int height2);

#endif