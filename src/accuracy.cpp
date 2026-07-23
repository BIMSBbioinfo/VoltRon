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
                              Rcpp::NumericMatrix trans_mat,
                              const int width1, 
                              const int height1,
                              const int width2, 
                              const int height2,                    
                              const bool invert_query, 
                              const bool invert_ref)
{
  // results
  Rcpp::List out(2);
  
  // Read images
  cv::Mat imReference = imageToMat(ref_image, width1, height1);
  cv::Mat imReg = imageToMat(query_image, width1, height1);
  cv::Mat h = numericMatrixToMat(trans_mat);
  
  // generate mask
  cv::Mat maskReg = generateOverlapMask(imReference.size(), 
                                        h,
                                        cv::Size(width2, height2));
  
  // process 
  Mat im1Proc, im2Proc;
  cvtColor(imReg, im1Proc, cv::COLOR_BGR2GRAY);
  cvtColor(imReference, im2Proc, cv::COLOR_BGR2GRAY);
  im1Proc = preprocessImage(im1Proc, invert_query, "None", "0");
  im2Proc = preprocessImage(im2Proc, invert_ref, "None", "0");
  
  // get metrics
  std::map<std::string, double> accuracy;
  accuracy = getAlignmentMetrics(im1Proc, im2Proc, maskReg, "Course");
  out[0] = accuracy;
  
  // get matte map
  Mat1d accuracyMatte;
  accuracyMatte = MatteMIMap(im2Proc, im1Proc, maskReg, 50);
  out[1] = matToNumericMatrix(accuracyMatte); // Matte MI metric
  
  // return
  return out;
}