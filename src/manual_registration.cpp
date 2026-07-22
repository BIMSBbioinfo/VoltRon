#include <Rcpp.h>

// OpenCV
#include <opencv2/opencv.hpp>
#include "opencv2/shape/shape_transformer.hpp"

// Library
#include "auxiliary.h"
#include "image.h"
#include "metrics.h"
#include "matte_mi.h"

// Namespaces
using namespace Rcpp;
using namespace std;
using namespace cv;

// align images with TPS algorithm
void alignImagesTPS(Mat &im1, Mat &im2, Mat &im1Reg, Rcpp::List &keypoints, 
                    Rcpp::NumericMatrix query_landmark, Rcpp::NumericMatrix reference_landmark,
                    const bool invert_query, const bool invert_ref,
                    Mat1d &accuracyMatte, 
                    std::map<std::string, double> &accuracy)
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

  // message
  Rcout << "Running Course Alignment (Thin-Plate-Spline)" << endl;

  // calculate transformation
  Ptr<ThinPlateSplineShapeTransformer> tps = cv::createThinPlateSplineShapeTransformer(0);
  tps->estimateTransformation(ref_mat, query_mat, matches);

  // save keypoints 
  keypoints[0] = point2fToNumericMatrix(ref_mat);
  keypoints[1] = point2fToNumericMatrix(query_mat); 
  
  // transform image using trained tps
  im1Reg = warpTPSImage(im2, im1, tps, 
                        im2.rows, im2.cols, cv::INTER_LINEAR);
  
  // // determine extension limits for both images
  // int y_max = max(im1.rows, im2.rows);
  // int x_max = max(im1.cols, im2.cols);
  // 
  // // extend images
  // cv::copyMakeBorder(im1, im1, 
  //                    0.0, (int) (y_max - im1.rows), 
  //                    0.0, (x_max - im1.cols), 
  //                    cv::BORDER_CONSTANT, 
  //                    Scalar(0, 0, 0));
  // 
  // // transform image
  // tps->warpImage(im1, im1Reg);
  // 
  // 
  // // resize image
  // cv::Mat im1Reg_cropped  = im1Reg(cv::Range(0,im2.size().height), 
  //                                  cv::Range(0,im2.size().width));
  // im1Reg = im1Reg_cropped.clone();
  
  // process 
  Mat im1Proc, im2Proc;
  cvtColor(im1Reg, im1Proc, cv::COLOR_BGR2GRAY);
  cvtColor(im2, im2Proc, cv::COLOR_BGR2GRAY);
  im1Proc = preprocessImage(im1Proc, invert_query, "None", "0");
  im2Proc = preprocessImage(im2Proc, invert_ref, "None", "0");
  
  // get alignment mask
  cv::Mat alignmentMask = generateOverlapMask(im2Proc,
                                              tps, 
                                              im1Proc.size());
  
  // get alignment metrics
  accuracy = getAlignmentMetrics(im1Proc, im2Proc, alignmentMask, "Course");
  accuracyMatte = MatteMIMap(im2Proc, im1Proc, alignmentMask, 50);
}

// align images with TPS algorithm
void alignImagesTPS_points(Rcpp::NumericMatrix &query_data,
                           Rcpp::NumericMatrix &dataReg,
                           Rcpp::List &keypoints, 
                           Rcpp::NumericMatrix query_landmark, 
                           Rcpp::NumericMatrix reference_landmark)
{
  // seed
  cv::setRNGSeed(0);
  RNG rng(12345);
  Scalar value;
  
  // Get landmarks as Point2f
  std::vector<cv::Point2f> query_mat = numericMatrixToPoint2f(query_landmark);
  std::vector<cv::Point2f> ref_mat = numericMatrixToPoint2f(reference_landmark);
  
  // Get data as Point2f
  std::vector<cv::Point2f> query_data_mat = numericMatrixToPoint2f(query_data);

  // get matches
  std::vector<cv::DMatch> matches;
  for (unsigned int i = 0; i < ref_mat.size(); i++)
    matches.push_back(cv::DMatch(i, i, 0));
  
  // calculate transformation
  Ptr<ThinPlateSplineShapeTransformer> tps = cv::createThinPlateSplineShapeTransformer(0);
  tps->estimateTransformation(ref_mat, query_mat, matches);
  
  // apply transformation to coordinates
  std::vector<cv::Point2f> query_data_reg;
  tps->applyTransformation(query_data_mat, query_data_reg);
  
  // save keypoints 
  keypoints[0] = point2fToNumericMatrix(ref_mat);
  keypoints[1] = point2fToNumericMatrix(query_mat); 
  
  // transform points
  dataReg = point2fToNumericMatrix(query_data_reg);
}

// align images with FLANN algorithm
void alignImagesAffineTPS(Mat &im1, Mat &im2, Mat &im1Reg, Mat &h, Rcpp::List &keypoints,
                          Rcpp::NumericMatrix query_landmark, Rcpp::NumericMatrix reference_landmark,
                          const bool invert_query, const bool invert_ref,
                          const bool run_Affine, const bool run_TPS,
                          Mat1d &accuracyMatte, 
                          std::map<std::string, double> &accuracy)
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
  Rcout << "Calculating" << (run_Affine ? " (Affine) " : " (Homography) ") << "Transformation Matrix" << endl;

  Mat im1Affine;
  std::vector<cv::Point2f> query_reg;
  if(run_Affine){
    h = estimateAffine2D(query_mat, ref_mat);
    cv::warpAffine(im1, im1Affine, h, im2.size());
    cv::transform(query_mat, query_reg, h);
  } else {
    h = findHomography(query_mat, ref_mat);
    cv::warpPerspective(im1, im1Affine, h, im2.size());
    cv::perspectiveTransform(query_mat, query_reg, h);
  }
  
  // get alignment metrics for course registration
  cv::Mat alignmentMask = generateOverlapMask(im2.size(), 
                                              h, 
                                              im1.size());

  // get matte metric, process image before
  Mat im1Proc, im2Proc;
  cvtColor(im1Affine, im1Proc, cv::COLOR_BGR2GRAY);
  cvtColor(im2, im2Proc, cv::COLOR_BGR2GRAY);
  im1Proc = preprocessImage(im1Proc, invert_query, "None", "0");
  im2Proc = preprocessImage(im2Proc, invert_ref, "None", "0");
  accuracy = getAlignmentMetrics(im1Proc, im2Proc, alignmentMask, "Course");
  accuracyMatte = MatteMIMap(im2Proc, im1Proc, alignmentMask, 50);
  
  if(!run_TPS){
    
    // clone and exit
    im1Reg = im1Affine.clone();
    
  } else {
    
    // message
    Rcout << "Running Fine Alignment (Thin-Plate-Spline)" << endl;

    // calculate TPS transformation
    Ptr<ThinPlateSplineShapeTransformer> tps = cv::createThinPlateSplineShapeTransformer(0);
    tps->estimateTransformation(ref_mat, query_reg, matches);
    
    // save keypoints 
    keypoints[0] = point2fToNumericMatrix(ref_mat);
    keypoints[1] = point2fToNumericMatrix(query_reg); 
    
    // transform image using trained tps
    im1Reg = warpTPSImage(im2, im1Affine, tps, 
                          im2.rows, im2.cols, cv::INTER_LINEAR);
    
    // // determine extension limits for both images
    // int y_max = max(im1Affine.rows, im2.rows);
    // int x_max = max(im1Affine.cols, im2.cols);
    // 
    // // extend images
    // cv::copyMakeBorder(im1Affine, im1Affine, 0.0, (int) (y_max - im1Affine.rows), 0.0, (x_max - im1Affine.cols), cv::BORDER_CONSTANT, Scalar(0, 0, 0));
    // 
    // // transform image
    // tps->warpImage(im1Affine, im1Reg);
    // 
    // // resize image
    // cv::Mat im1Reg_cropped  = im1Reg(cv::Range(0,im2.size().height), cv::Range(0,im2.size().width));
    // im1Reg = im1Reg_cropped.clone();
    
    // get alignment metrics
    cv::Mat alignmentMask = generateOverlapMask(im2, 
                                                tps, 
                                                im1Affine.size());
    
    // get matte metric, process 
    Mat im1Proc, im2Proc;
    cvtColor(im1Reg, im1Proc, cv::COLOR_BGR2GRAY);
    cvtColor(im2, im2Proc, cv::COLOR_BGR2GRAY);
    im1Proc = preprocessImage(im1Proc, invert_query, "None", "0");
    im2Proc = preprocessImage(im2Proc, invert_ref, "None", "0");
    accuracy = getAlignmentMetrics(im1Proc, im2Proc, alignmentMask, "Fine");
    accuracyMatte = MatteMIMap(im2Proc, im1Proc, alignmentMask, 50);
  }
}

// align images with FLANN algorithm
void alignImagesAffineTPS_points(Rcpp::NumericMatrix &query_data,
                                 Rcpp::NumericMatrix &dataReg,
                                 Mat &h, Rcpp::List &keypoints,
                                 Rcpp::NumericMatrix query_landmark, 
                                 Rcpp::NumericMatrix reference_landmark,
                                 const bool run_Affine, 
                                 const bool run_TPS)
{
  // seed
  cv::setRNGSeed(0);
  RNG rng(12345);
  Scalar value;
  
  // Get landmarks as Point2f
  std::vector<cv::Point2f> query_mat = numericMatrixToPoint2f(query_landmark);
  std::vector<cv::Point2f> ref_mat = numericMatrixToPoint2f(reference_landmark);
  
  // Get data as Point2f
  std::vector<cv::Point2f> query_data_mat = numericMatrixToPoint2f(query_data);

  // get matches
  std::vector<cv::DMatch> matches;
  for (unsigned int i = 0; i < ref_mat.size(); i++)
    matches.push_back(cv::DMatch(i, i, 0));
  
  // calculate homography transformation
  Rcout << "Calculating" << (run_Affine ? " (Affine) " : " (Homography) ") << "Transformation Matrix" << endl;

  std::vector<cv::Point2f> query_reg;
  std::vector<cv::Point2f> query_data_reg;
  if(run_Affine){
    h = estimateAffine2D(query_mat, ref_mat);
    cv::transform(query_mat, query_reg, h);
    cv::transform(query_data_mat, query_data_reg, h);
  } else {
    h = findHomography(query_mat, ref_mat);
    cv::perspectiveTransform(query_mat, query_reg, h);
    cv::perspectiveTransform(query_data_mat, query_data_reg, h);
  }

  if(run_TPS){
    
    // message
    Rcout << "Running Fine Alignment (Thin-Plate-Spline)" << endl;

    // calculate TPS transformation
    Ptr<ThinPlateSplineShapeTransformer> tps = cv::createThinPlateSplineShapeTransformer(0);
    tps->estimateTransformation(ref_mat, query_reg, matches);
      
    // apply transformation to coordinates
    tps->applyTransformation(query_data_reg, query_data_reg);
    
    // save keypoints 
    keypoints[0] = point2fToNumericMatrix(ref_mat);
    keypoints[1] = point2fToNumericMatrix(query_reg); 
  
  } 
  
  // save data  
  dataReg = point2fToNumericMatrix(query_data_reg);
}

// [[Rcpp::export]]
Rcpp::List manual_registeration_rawvector(Rcpp::RawVector ref_image, 
                                          Rcpp::RawVector query_image,
                                          Rcpp::NumericMatrix reference_landmark, 
                                          Rcpp::NumericMatrix query_landmark,
                                          const int width1, 
                                          const int height1,
                                          const int width2, 
                                          const int height2,
                                          const bool invert_query, 
                                          const bool invert_ref,
                                          Rcpp::String method, 
                                          Rcpp::String nonrigid)
{
  // Return data
  Rcpp::List out(4);
  Rcpp::List out_trans(2);
  Rcpp::List keypoints(2);
  Mat imReg, h;
  Mat1d accuracyMatte;
  std::map<std::string, double> accuracy;
  
  // get params
  const bool run_TPS = (strcmp(method.get_cstring(), "Homography + Non-Rigid") == 0 || 
                        strcmp(method.get_cstring(), "Affine + Non-Rigid") == 0) && 
                        strcmp(nonrigid.get_cstring(), "TPS (OpenCV)") == 0;
  const bool run_Affine = (strcmp(method.get_cstring(), "Affine") == 0 || 
                           strcmp(method.get_cstring(), "Affine + Non-Rigid") == 0);
  
  // Read reference and query images
  cv::Mat imReference = imageToMat(ref_image, width1, height1);
  cv::Mat im = imageToMat(query_image, width2, height2);
  
  // AffineHomography + Non-rigid (TPS)
  if(strcmp(method.get_cstring(), "Non-Rigid") != 0){
    alignImagesAffineTPS(im, imReference, imReg, 
                         h, keypoints,
                         query_landmark, reference_landmark, 
                         invert_query, 
                         invert_ref,
                         run_Affine, run_TPS, 
                         accuracyMatte, 
                         accuracy);
  }
  
  // Non-rigid (TPS) only
  if(strcmp(method.get_cstring(), "Non-Rigid") == 0){
    alignImagesTPS(im, imReference, imReg, 
                   keypoints, 
                   query_landmark, 
                   reference_landmark,
                   invert_query, 
                   invert_ref,
                   accuracyMatte, 
                   accuracy);
  }
  
  // transformation matrix, can be either a matrix, set of keypoints or both
  out_trans[0] = matToNumericMatrix(h.clone());
  out_trans[1] = keypoints;
  out[0] = out_trans;
  
  // registered image if exists
  out[1] = matToImage(imReg.clone()); 
  out[2] = matToNumericMatrix(accuracyMatte); // Matte MI metric
  out[3] = accuracy;
  
  return out;
}

// [[Rcpp::export]]
Rcpp::List manual_registeration_matrix(Rcpp::NumericMatrix query_data,
                                       Rcpp::NumericMatrix reference_landmark, 
                                       Rcpp::NumericMatrix query_landmark,
                                       Rcpp::String method, 
                                       Rcpp::String nonrigid)
{
  // Return data
  Rcpp::List out(2);
  Rcpp::List out_trans(2);
  Rcpp::List keypoints(2);
  Mat h;
  Rcpp::NumericMatrix dataReg;
  
  // get params
  const bool run_TPS = (strcmp(method.get_cstring(), "Homography + Non-Rigid") == 0 || 
                        strcmp(method.get_cstring(), "Affine + Non-Rigid") == 0) && 
                        strcmp(nonrigid.get_cstring(), "TPS (OpenCV)") == 0;
  const bool run_Affine = (strcmp(method.get_cstring(), "Affine") == 0 || 
                           strcmp(method.get_cstring(), "Affine + Non-Rigid") == 0);
 
 // AffineHomography + Non-rigid (TPS)
 if(strcmp(method.get_cstring(), "Non-Rigid") != 0){
   alignImagesAffineTPS_points(query_data, dataReg,
                               h, keypoints,
                               query_landmark, reference_landmark,
                               run_Affine, run_TPS);
   out_trans[0] = matToNumericMatrix(h.clone());
 }
 
 if(strcmp(method.get_cstring(), "Non-Rigid") == 0){
   alignImagesTPS_points(query_data, dataReg,
                         keypoints,
                         query_landmark, 
                         reference_landmark);
   keypoints[0] = keypoints[0];
   keypoints[1] = keypoints[1];
 } 
 
 // transformation matrix, can be either a matrix, set of keypoints or both
 out_trans[1] = keypoints;
 out[0] = out_trans;
 out[1] = dataReg;
 
 return out;
}
  