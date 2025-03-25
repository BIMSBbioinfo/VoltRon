#include <Rcpp.h>

// OpenCV
#include <opencv2/opencv.hpp>
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/shape/shape_transformer.hpp"
// #include <opencv2/imgproc.hpp>

// Auxiliary
#include "auxiliary.h"

// Namespaces
using namespace Rcpp;
using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

// align images with TPS algorithm
void alignImagesTPS(Mat &im1, Mat &im2, Mat &im1Reg, Rcpp::List &keypoints, 
                    Rcpp::NumericMatrix query_landmark, Rcpp::NumericMatrix reference_landmark)
{

  // seed
  cv::setRNGSeed(0);
  RNG rng(12345);
  Scalar value;

  // Get landmarks as Point2f
  std::vector<cv::Point2f> query_mat = numericMatrixToPoint2f(query_landmark);
  std::vector<cv::Point2f> ref_mat = numericMatrixToPoint2f(reference_landmark);

  // get matches
  std::vector<cv::DMatch> matches;
  for (unsigned int i = 0; i < ref_mat.size(); i++)
    matches.push_back(cv::DMatch(i, i, 0));

  // calculate transformation
  Ptr<ThinPlateSplineShapeTransformer> tps = cv::createThinPlateSplineShapeTransformer(0);
  tps->estimateTransformation(ref_mat, query_mat, matches);

  // save keypoints 
  keypoints[0] = point2fToNumericMatrix(ref_mat);
  keypoints[1] = point2fToNumericMatrix(query_mat); 
  
  // determine extension limits for both images
  int y_max = max(im1.rows, im2.rows);
  int x_max = max(im1.cols, im2.cols);

  // extend images
  cv::copyMakeBorder(im1, im1, 0.0, (int) (y_max - im1.rows), 0.0, (x_max - im1.cols), cv::BORDER_CONSTANT, Scalar(0, 0, 0));

  // transform image
  tps->warpImage(im1, im1Reg);

  // resize image
  cv::Mat im1Reg_cropped  = im1Reg(cv::Range(0,im2.size().height), cv::Range(0,im2.size().width));
  im1Reg = im1Reg_cropped.clone();
}

// align images with affine transformation with TPS algorithm
void alignImagesAffineTPS(Mat &im1, Mat &im2, Mat &im1Reg, Mat &h, Rcpp::List &keypoints,
                          Rcpp::NumericMatrix query_landmark, Rcpp::NumericMatrix reference_landmark, 
                          const bool run_Affine)
{
  
  // seed
  cv::setRNGSeed(0);
  RNG rng(12345);
  Scalar value;
  
  // Get landmarks as Point2f
  std::vector<cv::Point2f> query_mat = numericMatrixToPoint2f(query_landmark);
  std::vector<cv::Point2f> ref_mat = numericMatrixToPoint2f(reference_landmark);
  
  // get matches
  std::vector<cv::DMatch> matches;
  for (unsigned int i = 0; i < ref_mat.size(); i++)
    matches.push_back(cv::DMatch(i, i, 0));
  
  // calculate homography transformation
  Mat im1Affine;
  std::vector<cv::Point2f> query_reg;
  if(run_Affine){
    h = estimateAffine2D(query_mat, ref_mat);
    Rcout << h << endl;
    cv::warpAffine(im1, im1Affine, h, im2.size());
    cv::transform(query_mat, query_reg, h);
  } else {
    h = findHomography(query_mat, ref_mat);
    Rcout << h << endl;
    cv::warpPerspective(im1, im1Affine, h, im2.size());
    cv::perspectiveTransform(query_mat, query_reg, h);
  }
  
  // calculate TPS transformation
  Ptr<ThinPlateSplineShapeTransformer> tps = cv::createThinPlateSplineShapeTransformer(0);
  tps->estimateTransformation(ref_mat, query_reg, matches);
  
  // save keypoints 
  keypoints[0] = point2fToNumericMatrix(ref_mat);
  keypoints[1] = point2fToNumericMatrix(query_reg); 
  
  // determine extension limits for both images
  int y_max = max(im1Affine.rows, im2.rows);
  int x_max = max(im1Affine.cols, im2.cols);

  // extend images
  cv::copyMakeBorder(im1Affine, im1Affine, 0.0, (int) (y_max - im1Affine.rows), 0.0, (x_max - im1Affine.cols), cv::BORDER_CONSTANT, Scalar(0, 0, 0));

  // transform image
  tps->warpImage(im1Affine, im1Reg);
  
  // resize image
  cv::Mat im1Reg_cropped  = im1Reg(cv::Range(0,im2.size().height), cv::Range(0,im2.size().width));
  im1Reg = im1Reg_cropped.clone();
}

// [[Rcpp::export]]
Rcpp::List manual_registeration_rawvector(Rcpp::RawVector ref_image, Rcpp::RawVector query_image,
                                          Rcpp::NumericMatrix reference_landmark, Rcpp::NumericMatrix query_landmark,
                                          const int width1, const int height1,
                                          const int width2, const int height2,
                                          Rcpp::String method)
{
  // Return data
  Rcpp::List out(2);
  Rcpp::List out_trans(2);
  Rcpp::List keypoints(2);
  Mat imReg, h;
  
  // Read reference image
  cv::Mat imReference = imageToMat(ref_image, width1, height1);

  // Read image to be aligned
  cv::Mat im = imageToMat(query_image, width2, height2);


  // Homography + Non-rigid (TPS)
  if(strcmp(method.get_cstring(), "Homography + Non-Rigid") == 0 || strcmp(method.get_cstring(), "Affine + Non-Rigid") == 0){
    const bool run_Affine = (strcmp(method.get_cstring(), "Affine + Non-Rigid") == 0);
    alignImagesAffineTPS(im, imReference, imReg, 
                         h, keypoints,
                         query_landmark, reference_landmark, run_Affine);
  }
  
  // Non-rigid (TPS) only
  if(strcmp(method.get_cstring(), "Non-Rigid") == 0){
    alignImagesTPS(im, imReference, imReg, 
                   keypoints, 
                   query_landmark, reference_landmark);
  }

  // transformation matrix, can be either a matrix, set of keypoints or both
  out_trans[0] = matToNumericMatrix(h.clone());
  out_trans[1] = keypoints;
  out[0] = out_trans;
  
  // registered image
  out[1] = matToImage(imReg.clone());

  return out;
}