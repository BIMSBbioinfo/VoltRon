#include <Rcpp.h>

// OpenCV
#include <opencv2/opencv.hpp>
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/shape/shape_transformer.hpp"

// Internal functions
#include "auxiliary.h"

using namespace Rcpp;
using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

// [[Rcpp::export]]
Rcpp::NumericMatrix applyMapping(Rcpp::NumericMatrix coords, Rcpp::List mapping)
{
  // Get coordinates as Point2f
  std::vector<cv::Point2f> coords_mat = numericMatrixToPoint2f(coords);
  std::vector<cv::Point2f> coords_temp;
  
  // seed
  cv::setRNGSeed(0);
  RNG rng(12345);
  Scalar value;
  
  // list 
  int n = mapping.size();
  
  // iterate over the list 
  for(int i=0; i<n; i++) {
    
    // get the current mapping
    Rcpp::List cur_mapping = mapping[i];
    
    // Get transformation matrix
    cv::Mat h = numericMatrixToMat(cur_mapping[0]);
    if(h.cols > 0){
      
      // transform coordinates
      cv::perspectiveTransform(coords_mat, coords_temp, h);
      coords_mat = coords_temp;
      
    } 
    
    // non-rigid warping
    if(cur_mapping[1] != R_NilValue){
      
      // Get landmarks as Point2f 
      Rcpp::List keypoints = cur_mapping[1];
      std::vector<cv::Point2f> ref_mat = numericMatrixToPoint2f(keypoints[0]);
      std::vector<cv::Point2f> query_mat = numericMatrixToPoint2f(keypoints[1]);
      
      // get matches
      std::vector<cv::DMatch> matches;
      for (unsigned int i = 0; i < ref_mat.size(); i++)
        matches.push_back(cv::DMatch(i, i, 0));
      
      // calculate transformation
      Ptr<ThinPlateSplineShapeTransformer> tps = cv::createThinPlateSplineShapeTransformer(0);
      tps->estimateTransformation(query_mat, ref_mat, matches);
      
      // apply transformation to coordinates
      tps->applyTransformation(coords_mat, coords_temp);
      
    } else {
      
      coords_temp = coords_mat;
    }
    
    // set initial coordinates and registered one
    coords_mat = coords_temp;
  }
  
  // return registered coordinates as numeric matrix
  Rcpp::NumericMatrix coords_regToMat = point2fToNumericMatrix(coords_mat);
  
  // return
  return coords_regToMat;
}