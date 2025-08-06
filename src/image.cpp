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

////
// processing
////

// process images before keypoint detection
cv::Mat preprocessImage(Mat &im, const bool invert, const char* flipflop, const char* rotate)
{
  // normalize
  Mat imNorm;
  cv::normalize(im, imNorm, 0, 255, cv::NORM_MINMAX, CV_8UC1);
  
  // rotate image
  Mat imRotate;
  if(atoi(rotate) > 0){
    cv::rotate(imNorm, imRotate, (atoi(rotate)/90)-1);
  } else {
    imRotate = imNorm;
  }
  
  // Flipflop image
  Mat imFlipFlop;
  if(strcmp(flipflop, "Flip") == 0){
    cv::flip(imRotate, imFlipFlop, 0);
  } else if(strcmp(flipflop, "Flop") == 0){
    cv::flip(imRotate, imFlipFlop, 1);
  } else if(strcmp(flipflop, "None") == 0){
    imFlipFlop = imRotate;
  }
  
  // invert/negate image and full processed image
  Mat imProcess;
  if(invert) {
    cv::bitwise_not(imFlipFlop, imProcess);
  } else {
    imProcess = imFlipFlop;
  }
  
  // return
  return imProcess;
}

// revert the processing on registrated images using the reference image
cv::Mat reversepreprocessImage(Mat &im, const char* flipflop, const char* rotate)
{
  
  // Flipflop image
  Mat imFlipFlop;
  if(strcmp(flipflop, "Flip") == 0){
    cv::flip(im, imFlipFlop, 0);
  } else if(strcmp(flipflop, "Flop") == 0){
    cv::flip(im, imFlipFlop, 1);
  } else if(strcmp(flipflop, "None") == 0){
    imFlipFlop = im;
  }
  
  // rotate image
  Mat imRotate;
  if(atoi(rotate) > 0){
    cv::rotate(imFlipFlop, imRotate, ((360-atoi(rotate))/90)-1);
  } else {
    imRotate = imFlipFlop;
  }
  
  // return
  return imRotate;
}

cv::Mat resize_image(cv::Mat &im, int width)
{
  // calculate resize factor
  int resized_cols = width;
  int resized_rows = (int) round(im.rows * ((double)width / im.cols));
  
  // resize image
  Mat im_resized;
  cv::resize(im, im_resized, cv::Size(resized_cols,resized_rows),0,0,cv::INTER_AREA);
  return im_resized;
}

std::vector<cv::KeyPoint> resize_keypoints(std::vector<cv::KeyPoint> &keypoints, 
                                           cv::Mat &im,
                                           int width)
{
  // scale factor
  float factor = (float) width / im.cols;
  
  // scale keypoints
  std::vector<cv::KeyPoint> new_keypoints;
  for (int i = 0; i < keypoints.size(); i++) {
    cv::KeyPoint kp = keypoints[i];
    kp.pt.x *= factor;
    kp.pt.y *= factor;
    new_keypoints.push_back(kp);
  }
  
  // return
  return new_keypoints;
}

void scaledDrawMatches(cv::Mat im1, std::vector<cv::KeyPoint> &keypoints1, 
                       cv::Mat im2, std::vector<cv::KeyPoint> &keypoints2,
                       std::vector<cv::DMatch> &top_matches,
                       cv::Mat &imMatches)
{
  // resize keypoints
  keypoints1 = resize_keypoints(keypoints1, im1, 500);
  keypoints2 = resize_keypoints(keypoints2, im2, 500);

  // resize images
  im1 = resize_image(im1, 500);
  im2 = resize_image(im2, 500);
  
  // draw matches
  drawMatches(im1, keypoints1, im2, keypoints2, top_matches, imMatches);
}
  
// [[Rcpp::export]]
Rcpp::RawVector warpImage(Rcpp::RawVector ref_image, Rcpp::RawVector query_image, 
                          Rcpp::List mapping,
                          const int width1, const int height1,
                          const int width2, const int height2)
{
  // Read reference image
  cv::Mat imReference = imageToMat(ref_image, width1, height1);
  
  // Read image to be aligned
  cv::Mat im = imageToMat(query_image, width2, height2);
  cv::Mat im_temp;
  
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
      if(h.rows == 2){
        cv::warpAffine(im, im_temp, h, imReference.size()); 
      } else {
        cv::warpPerspective(im, im_temp, h, imReference.size());    
      }
      im = im_temp;
    } 
    
    // non-rigid warping
    if(cur_mapping[1] != R_NilValue){
      
      // get landmarks
      Rcpp::List keypoints = cur_mapping[1];
      std::vector<cv::Point2f> ref_mat = numericMatrixToPoint2f(keypoints[0]);
      std::vector<cv::Point2f> query_mat = numericMatrixToPoint2f(keypoints[1]);
      
      // get matches
      std::vector<cv::DMatch> matches;
      for (unsigned int i = 0; i < ref_mat.size(); i++)
        matches.push_back(cv::DMatch(i, i, 0));
      
      // calculate transformation
      Ptr<ThinPlateSplineShapeTransformer> tps = cv::createThinPlateSplineShapeTransformer(0);
      tps->estimateTransformation(ref_mat, query_mat, matches);
      
      // determine extension limits for both images
      int y_max = max(im.rows, imReference.rows);
      int x_max = max(im.cols, imReference.cols);
      
      // extend images
      cv::copyMakeBorder(im, im, 0.0, (int) (y_max - im.rows), 0.0, (x_max - im.cols), cv::BORDER_CONSTANT, Scalar(0, 0, 0));
      
      // transform image
      tps->warpImage(im, im_temp);
      
      // resize image
      cv::Mat im_temp_cropped  = im_temp(cv::Range(0,imReference.size().height), cv::Range(0,imReference.size().width));
      im_temp = im_temp_cropped.clone();
      
    } else {
      
      // pass registered object
      im_temp = im;
    }
    
    im = im_temp;
  }
  
  // return
  return matToImage(im);
} 


// [[Rcpp::export]]
Rcpp::RawVector warpImageAuto(Rcpp::RawVector ref_image, Rcpp::RawVector query_image, 
                              Rcpp::List mapping,
                              const int width1, const int height1,
                              const int width2, const int height2)
{
  // Get transformation matrix
  cv::Mat h = numericMatrixToMat(mapping[0]);
  
  // Read reference image
  cv::Mat imReference = imageToMat(ref_image, width1, height1);
  
  // Read image to be aligned
  cv::Mat im = imageToMat(query_image, width2, height2);
  
  // transform coordinates
  Mat imWarp;
  cv::warpPerspective(im, imWarp, h, imReference.size()); 
  
  // non-rigid warping
  cv::Mat imReg;
  if(mapping[1] != R_NilValue){
    
    // get landmarks
    Rcpp::List keypoints = mapping[1];
    std::vector<cv::Point2f> ref_mat = numericMatrixToPoint2f(keypoints[0]);
    std::vector<cv::Point2f> query_mat = numericMatrixToPoint2f(keypoints[1]);
    
    // get matches
    std::vector<cv::DMatch> matches;
    for (unsigned int i = 0; i < ref_mat.size(); i++)
      matches.push_back(cv::DMatch(i, i, 0));
    
    // calculate transformation
    Ptr<ThinPlateSplineShapeTransformer> tps = cv::createThinPlateSplineShapeTransformer(0);
    tps->estimateTransformation(ref_mat, query_mat, matches);
    
    // determine extension limits for both images
    int y_max = max(imWarp.rows, imReference.rows);
    int x_max = max(imWarp.cols, imReference.cols);
    
    // extend images
    cv::copyMakeBorder(imWarp, imWarp, 0.0, (int) (y_max - imWarp.rows), 0.0, (x_max - imWarp.cols), cv::BORDER_CONSTANT, Scalar(0, 0, 0));
    
    // transform image
    tps->warpImage(imWarp, imReg);
    
    // resize image
    // cv::resize(im1Reg, im1Reg, im2.size());
    cv::Mat imReg_cropped  = imReg(cv::Range(0,imReference.size().height), cv::Range(0,imReference.size().width));
    imReg = imReg_cropped.clone();
    
  } else {
    
    // pass registered object
    imReg = imWarp;
  }
  
  // return
  return matToImage(imReg);
}


// [[Rcpp::export]]
Rcpp::RawVector warpImageManual(Rcpp::RawVector ref_image, Rcpp::RawVector query_image, 
                                Rcpp::List mapping,
                                const int width1, const int height1,
                                const int width2, const int height2)
{
  // Get landmarks as Point2f and transformation matrix
  cv::Mat h = numericMatrixToMat(mapping[0]);
  Rcpp::List keypoints = mapping[1];
  std::vector<cv::Point2f> ref_mat = numericMatrixToPoint2f(keypoints[0]);
  std::vector<cv::Point2f> query_mat = numericMatrixToPoint2f(keypoints[1]);
  
  // Read reference image
  cv::Mat imReference = imageToMat(ref_image, width1, height1);
  
  // Read image to be aligned
  cv::Mat im = imageToMat(query_image, width2, height2);
  
  // transform coordinates
  Mat imWarp;
  if(h.cols > 0){
    cv::warpPerspective(im, imWarp, h, imReference.size()); 
  } else {
    imWarp = im;
  }
  
  // get matches
  std::vector<cv::DMatch> matches;
  for (unsigned int i = 0; i < ref_mat.size(); i++)
    matches.push_back(cv::DMatch(i, i, 0));
  
  // calculate transformation
  Ptr<ThinPlateSplineShapeTransformer> tps = cv::createThinPlateSplineShapeTransformer(0);
  tps->estimateTransformation(ref_mat, query_mat, matches);
  
  // determine extension limits for both images
  int y_max = max(imWarp.rows, imReference.rows);
  int x_max = max(imWarp.cols, imReference.cols);
  
  // extend images
  cv::copyMakeBorder(imWarp, imWarp, 0.0, (int) (y_max - imWarp.rows), 0.0, (x_max - imWarp.cols), cv::BORDER_CONSTANT, Scalar(0, 0, 0));
  
  // transform image
  cv::Mat imReg;
  tps->warpImage(imWarp, imReg);
  
  // resize image
  // cv::resize(im1Reg, im1Reg, im2.size());
  cv::Mat imReg_cropped  = imReg(cv::Range(0,imReference.size().height), cv::Range(0,imReference.size().width));
  imReg = imReg_cropped.clone();
  
  // return
  return matToImage(imReg);;
}