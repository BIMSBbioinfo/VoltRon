#include <Rcpp.h>

// OpenCV
#include <opencv2/opencv.hpp>

// Library
#include "auxiliary.h"
#include "image.h"
#include "metrics.h"
#include "matte_mi.h"

// Namespaces
using namespace Rcpp;
using namespace std;
using namespace cv;

// [[Rcpp::export]]
Rcpp::List accuracy_rawvector(Rcpp::RawVector& ref_image, 
                              Rcpp::RawVector& query_image,
                              Rcpp::RawVector& mask,
                              const int width, 
                              const int height,
                              std::string type,
                              bool overlay_images = true,
                              const bool compute_matte_map = true) {
  // results
  Rcpp::List out(3);
  
  // Read images
  cv::Mat imReference = imageToMat(ref_image, width, height);
  cv::Mat imReg = imageToMat(query_image, width, height);
  cv::Mat maskReg = imageToMat(mask, width, height);

  // process 
  Mat im1Proc, im2Proc;
  cvtColor(imReg, im1Proc, cv::COLOR_BGR2GRAY);
  cvtColor(imReference, im2Proc, cv::COLOR_BGR2GRAY);
  cvtColor(maskReg, maskReg, cv::COLOR_BGR2GRAY);
  
  // get metrics
  std::map<std::string, double> accuracy;
  accuracy = getAlignmentMetrics(im1Proc, im2Proc, maskReg, type);
  out[0] = accuracy;
  
  // get matte map
  Mat1d accuracyMatte;
  if(compute_matte_map){
    // accuracyMatte = MatteMIMap(im2Proc, im1Proc, maskReg, 50);
    accuracyMatte = getTiledAlignmentMetrics(im2Proc, im1Proc, maskReg);
    out[1] = matToNumericMatrix(accuracyMatte); // Matte MI metric 
  } else {
    out[1] = R_NilValue;
  }
  
  // image overlay
  if(overlay_images){
    cv::addWeighted(im2Proc, 0.7, im1Proc, 0.3, 0, im1Proc);
    cvtColor(im1Proc, im1Proc, cv::COLOR_GRAY2BGR);
    im1Proc = resize_image(im1Proc, 500);
    out[2] = matToImage(im1Proc); 
  } else {
    out[2] = R_NilValue;
  }
  
  // return
  return out;
}