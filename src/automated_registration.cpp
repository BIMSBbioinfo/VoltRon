#include <Rcpp.h>

// OpenCV
#include <opencv2/opencv.hpp>
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/shape/shape_transformer.hpp"
// #include <opencv2/imgproc.hpp>

// Internal functions
#include "auxiliary.h"
#include "image.h"

// Namespaces
using namespace Rcpp;
using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

// Parameters
struct Parameters
{ 
  const int sift_nfeatures=20000;
  const int tile_size=5000;
  const int tile_overlap=50;
};

// check if keypoints are degenerate
bool check_degenerate(std::vector<cv::KeyPoint> keypoints1, std::vector<cv::KeyPoint> keypoints2) {

  // get sd
  double keypoints1_sd = cppSD(keypoints1);
  double keypoints2_sd = cppSD(keypoints2);
  
  // get warning message
  bool is_degenerate = FALSE;
  if(keypoints1_sd < 1.0 | keypoints2_sd < 1.0){
    is_degenerate = TRUE;
    Rcout << "WARNING: points may be in a degenerate configuration." << endl;
  } 

  return is_degenerate;
}

// check distribution of registered points
std::string check_transformation_by_point_distribution(Mat &im, Mat &h){

  // get image shape
  int height = im.rows;
  int width = im.cols;
  int height_interval = (double) height/50.0;
  int width_interval = (double) width/50.0;
  
  // perspective transformation of grid points
  std::vector<cv::Point2f> gridpoints;
  for (double i = 0.0; i <= height; i += height_interval) {
    for (double j = 0.0; j <= width; j += width_interval) {
      gridpoints.push_back(cv::Point2f(j,i));
    }
  }
  
  // register grid points
  std::vector<cv::Point2f> gridpoints_reg;
  cv::perspectiveTransform(gridpoints, gridpoints_reg, h);

  // Compute the standard deviation of the transformed points
  double gridpoints_reg_sd = cppSD(gridpoints_reg);

  // get warning message
  std::string message;
  if(gridpoints_reg_sd < 1.0 | gridpoints_reg_sd > max(height, width)){
    message = "large distribution";
    Rcout << "WARNING: Transformation may be poor - transformed points grid seem to be concentrated!" << endl;
  } else {
    message = "small distribution";
  }

  return message;
}


// do overall checks on keypoints and images
bool check_transformation_metrics(std::vector<cv::KeyPoint> keypoints1, std::vector<cv::KeyPoint> keypoints2, Mat &im1, Mat &im2, Mat &h, Mat &mask) {
  
  // check keypoint standard deviation
  bool is_degenerate = check_degenerate(keypoints1, keypoints2);
  
  //  check transformation
  // std::string transformation;
  // transformation = check_transformation_by_pts_mean_sqrt(keypoints1, keypoints2, h, mask);
  
  //  check distribution
  std::string distribution;
  distribution = check_transformation_by_point_distribution(im2, h);
    
  return is_degenerate;
}

// get good matching keypoints
void getGoodMatches(std::vector<std::vector<DMatch>> matches, std::vector<DMatch> &good_matches, const float lowe_ratio = 0.8)
{
  for (size_t i = 0; i < matches.size(); i++) {
    if (matches[i][0].distance < lowe_ratio * matches[i][1].distance) {
      good_matches.push_back(matches[i][0]);
    }
  }
}

// remove duplicate keypoints for TPS
void removeCloseMatches(std::vector<cv::Point2f>& points1, std::vector<cv::Point2f>& points2, float threshold = std::numeric_limits<float>::epsilon()) {
  
  // Create a vector to store filtered points
  std::vector<cv::Point2f> filtered_points1;
  std::vector<cv::Point2f> filtered_points2;
  
  // Iterate through the points
  for (size_t i = 0; i < points1.size(); ++i) {
    bool is_close = false;
    
    // Compare the current point with all already filtered points
    for (size_t j = 0; j < filtered_points1.size(); ++j) {
      float dist_squared1 = (points1[i].x - filtered_points1[j].x) * (points1[i].x - filtered_points1[j].x) +
        (points1[i].y - filtered_points1[j].y) * (points1[i].y - filtered_points1[j].y);
      float dist_squared2 = (points2[i].x - filtered_points2[j].x) * (points2[i].x - filtered_points2[j].x) +
        (points2[i].y - filtered_points2[j].y) * (points2[i].y - filtered_points2[j].y);
      
      // If the distance is less than the threshold (considered close), break
      if (dist_squared1 < threshold | dist_squared2 < threshold) {
        is_close = true;
        break;
      }
    }
    
    // If no close point was found, add the current point to the filtered list
    if (!is_close) {
      filtered_points1.push_back(points1[i]);
      filtered_points2.push_back(points2[i]);
    }
  }
  
  // Update the original point set with the filtered points
  points1 = filtered_points1;
  points2 = filtered_points2;
}

// align images with BRUTE FORCE algorithm
void alignImagesBRUTE(Mat &im1, Mat &im2, Mat &im1Reg, Mat &im1Overlay, Mat &imMatches, Mat &h, const float GOOD_MATCH_PERCENT, const int MAX_FEATURES)
{
  // Convert images to grayscale
  Mat im1Gray, im2Gray;
  cvtColor(im1, im1Gray, cv::COLOR_BGR2GRAY);
  cvtColor(im2, im2Gray, cv::COLOR_BGR2GRAY);

  // Variables to store keypoints and descriptors
  std::vector<KeyPoint> keypoints1, keypoints2;
  Mat descriptors1, descriptors2;

  // Detect ORB features and compute descriptors.
  Ptr<Feature2D> orb = ORB::create(MAX_FEATURES);
  orb->detectAndCompute(im1Gray, Mat(), keypoints1, descriptors1);
  orb->detectAndCompute(im2Gray, Mat(), keypoints2, descriptors2);
  Rcout << "DONE: orb based key-points detection and descriptors computation" << endl;

  // Match features.
  std::vector<DMatch> matches;
  Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("BruteForce-Hamming");
  matcher->match(descriptors1, descriptors2, matches, Mat());
  Rcout << "DONE: BruteForce-Hamming - descriptor matching" << endl;

  // Sort matches by score
  std::sort(matches.begin(), matches.end());

  // Remove not so good matches
  const int numGoodMatches = matches.size() * GOOD_MATCH_PERCENT;
  matches.erase(matches.begin()+numGoodMatches, matches.end());
  Rcout << "DONE: get good matches by distance thresholding" << endl;

  // Extract location of good matches
  std::vector<Point2f> points1, points2;
  for( size_t i = 0; i < matches.size(); i++ )
  {
    points1.push_back( keypoints1[ matches[i].queryIdx ].pt );
    points2.push_back( keypoints2[ matches[i].trainIdx ].pt );
  }

  // Find homography
  h = findHomography(points1, points2, RANSAC);
  Rcout << "DONE: calculated homography matrix" << endl;

  // Extract location of good matches in terms of keypoints
  std::vector<KeyPoint> keypoints1_best, keypoints2_best;
  std::vector<cv::DMatch> goodMatches;
  for( size_t i = 0; i < matches.size(); i++ )
  {
    keypoints1_best.push_back(keypoints1[matches[i].queryIdx]);
    keypoints2_best.push_back(keypoints2[matches[i].trainIdx]);
    goodMatches.push_back(cv::DMatch(static_cast<int>(i), static_cast<int>(i), 0));
  }

  // Draw top matches and good ones only
  // Mat imMatches;
  drawMatches(im1, keypoints1_best, im2, keypoints2_best, goodMatches, imMatches);

  // Use homography to warp image
  warpPerspective(im1, im1Reg, h, im2.size());
  warpPerspective(im1, im1Overlay, h, im2.size());
}

void computeSIFTTiles(Mat &im, std::vector<KeyPoint> &keypoints, Mat &descriptors, Ptr<Feature2D> &sift,
                      Parameters params){
  
  // image size
  int width = im.size().width;
  int height = im.size().height;
  
  // tile overlaps
  int step_x = params.tile_size - params.tile_overlap;
  int step_y = params.tile_size - params.tile_overlap;
  
  // compare size and tiles
  if(width <= params.tile_size or height <= params.tile_size){
    
    sift->detectAndCompute(im, Mat(), keypoints, descriptors);
    return;
    
  } else {
    
    // Loop over tile positions
    for (int y = 0; y < height - params.tile_overlap; y += step_y) {
      for (int x = 0; x < width - params.tile_overlap; x += step_x) {
        
        // Define the region of interest (tile)
        int tile_x_end = std::min(x + params.tile_size, width);
        int tile_y_end = std::min(y + params.tile_size, height);
        cv::Rect tile_rect(x, y, tile_x_end - x, tile_y_end - y);
        
        // 
        cv::Mat tile = im(tile_rect);
        
        // Detect keypoints and compute descriptors in the tile
        std::vector<cv::KeyPoint> tile_keypoints;
        cv::Mat tile_descriptors;
        sift->detectAndCompute(tile, Mat(), tile_keypoints, tile_descriptors);
        
        // Adjust keypoint coordinates
        for (auto& kp : tile_keypoints) {
          kp.pt.x += x;
          kp.pt.y += y;
          keypoints.push_back(kp);
        }
        
        // Append tile descriptors to the global descriptors matrix
        if (!tile_descriptors.empty()) {
          if (descriptors.empty()) {
            descriptors = tile_descriptors.clone();
          } else {
            cv::vconcat(descriptors, tile_descriptors, descriptors);
          }
        }
        
      }
    }
    
  }
}

void getSIFTTransformationMatrix(Mat im1Proc, Mat im1NormalProc, Mat im2Proc, Mat im1, Mat im2, Mat &h, Mat &mask, 
                                 Mat &imMatches, std::vector<Point2f> &points1, std::vector<Point2f> &points2,
                                 Parameters params, bool &is_faulty){
  
  //////////////////////
  /// Compute SIFT /////
  //////////////////////
  
  // Variables to store keypoints and descriptors
  std::vector<KeyPoint> keypoints1, keypoints2;
  Mat descriptors1, descriptors2;
  
  // Detect SIFT features
  Ptr<Feature2D> sift = cv::SIFT::create(params.sift_nfeatures);
  // sift->detectAndCompute(im1Proc, Mat(), keypoints1, descriptors1);
  // sift->detectAndCompute(im2Proc, Mat(), keypoints2, descriptors2);
  computeSIFTTiles(im1Proc, keypoints1, descriptors1, sift, params);
  computeSIFTTiles(im2Proc, keypoints2, descriptors2, sift, params);
  Rcout << "DONE: sift based key-points detection and descriptors computation" << endl;
  
  ///////////////////////
  /// Compute FLANN /////
  ///////////////////////
  
  // Match features using FLANN matching
  std::vector<std::vector<DMatch>> matches;
  cv::FlannBasedMatcher custom_matcher = cv::FlannBasedMatcher(cv::makePtr<cv::flann::KDTreeIndexParams>(5), cv::makePtr<cv::flann::SearchParams>(50, 0, TRUE));
  cv::Ptr<cv::FlannBasedMatcher> matcher = custom_matcher.create();
  matcher->knnMatch(descriptors1, descriptors2, matches, 2);
  Rcout << "DONE: FLANN - Fast Library for Approximate Nearest Neighbors - descriptor matching" << endl;
  
  // Find good matches
  std::vector<DMatch> good_matches;
  getGoodMatches(matches, good_matches);
  Rcout << "DONE: get good matches by distance thresholding" << endl;
  
  ///////////////////////
  /// Find Homography ///
  ///////////////////////
  
  // Extract location of good matches
  //std::vector<Point2f> points1, points2;
  for( size_t i = 0; i < good_matches.size(); i++ )
  {
    points1.push_back(keypoints1[good_matches[i].queryIdx].pt);
    points2.push_back(keypoints2[good_matches[i].trainIdx].pt);
  }
  
  // Find homography
  h = findHomography(points1, points2, RANSAC, 5, mask);
  Rcout << "DONE: calculated homography matrix" << endl;
  
  ///////////////////////
  /// Find Matches ///
  ///////////////////////
  
  // Draw top matches and good ones only
  std::vector<cv::DMatch> top_matches;
  std::vector<cv::KeyPoint> keypoints1_best, keypoints2_best;
  for(size_t i = 0; i < good_matches.size(); i++ )
  {
    keypoints1_best.push_back(keypoints1[good_matches[i].queryIdx]);
    keypoints2_best.push_back(keypoints2[good_matches[i].trainIdx]);
    top_matches.push_back(cv::DMatch(static_cast<int>(i), static_cast<int>(i), 0));
  }
  drawMatches(im1Proc, keypoints1_best, im2Proc, keypoints2_best, top_matches, imMatches);
  
  // check keypoints
  is_faulty = check_transformation_metrics(keypoints1_best, keypoints2_best, im1, im2, h, mask);
}


// align images with FLANN algorithm
void alignImagesFLANN_old(Mat &im1, Mat &im2, Mat &im1Reg, Mat &im1Overlay, Mat &imMatches, Mat &h,
                 const bool invert_query, const bool invert_ref,
                 const char* flipflop_query, const char* flipflop_ref,
                 const char* rotate_query, const char* rotate_ref)
{

  // seed
  cv::setRNGSeed(0);

  // Convert images to grayscale
  Mat im1Gray, im2Gray;
  cvtColor(im1, im1Gray, cv::COLOR_BGR2GRAY);
  cvtColor(im2, im2Gray, cv::COLOR_BGR2GRAY);

  // Process images
  Mat im1Proc, im2Proc, im1NormalProc;
  im1Proc = preprocessImage(im1Gray, invert_query, flipflop_query, rotate_query);
  im1NormalProc = preprocessImage(im1, FALSE, flipflop_query, rotate_query);
  im2Proc = preprocessImage(im2Gray, invert_ref, flipflop_ref, rotate_ref);

  // Variables to store keypoints and descriptors
  std::vector<KeyPoint> keypoints1, keypoints2;
  Mat descriptors1, descriptors2;

  // Detect SIFT features
  Ptr<Feature2D> sift = cv::SIFT::create();
  sift->detectAndCompute(im1Proc, Mat(), keypoints1, descriptors1);
  sift->detectAndCompute(im2Proc, Mat(), keypoints2, descriptors2);
  Rcout << "DONE: sift based key-points detection and descriptors computation" << endl;

  // Match features using FLANN matching
  std::vector<std::vector<DMatch>> matches;
  cv::FlannBasedMatcher custom_matcher = cv::FlannBasedMatcher(cv::makePtr<cv::flann::KDTreeIndexParams>(5), cv::makePtr<cv::flann::SearchParams>(50, 0, TRUE));
  cv::Ptr<cv::FlannBasedMatcher> matcher = custom_matcher.create();
  matcher->knnMatch(descriptors1, descriptors2, matches, 2);
  Rcout << "DONE: FLANN - Fast Library for Approximate Nearest Neighbors - descriptor matching" << endl;

  // Find good matches
  // goodMatches = get_good_matches(matches)
  std::vector<DMatch> good_matches;
  getGoodMatches(matches, good_matches);
  Rcout << "DONE: get good matches by distance thresholding" << endl;

  // Extract location of good matches
  std::vector<Point2f> points1, points2;
  for( size_t i = 0; i < good_matches.size(); i++ )
  {
    points1.push_back(keypoints1[good_matches[i].queryIdx].pt);
    points2.push_back(keypoints2[good_matches[i].trainIdx].pt);
  }

  // Find homography
  cv::Mat mask;
  h = findHomography(points1, points2, RANSAC, 5, mask);
  Rcpp::Rcout << "DONE: calculated homography matrix" << endl;

  // Draw top matches and good ones only
  std::vector<cv::DMatch> top_matches;
  std::vector<cv::KeyPoint> keypoints1_best, keypoints2_best;
  for(size_t i = 0; i < good_matches.size(); i++ )
  {
    keypoints1_best.push_back(keypoints1[good_matches[i].queryIdx]);
    keypoints2_best.push_back(keypoints2[good_matches[i].trainIdx]);
    top_matches.push_back(cv::DMatch(static_cast<int>(i), static_cast<int>(i), 0));
  }
  drawMatches(im1Proc, keypoints1_best, im2Proc, keypoints2_best, top_matches, imMatches);
  
  // check keypoints
  bool is_faulty = check_transformation_metrics(keypoints1_best, keypoints2_best, im1, im2, h, mask);

  // Use homography to warp image
  Mat im1Warp, im1NormalWarp;
  warpPerspective(im1Proc, im1Warp, h, im2Proc.size());
  warpPerspective(im1NormalProc, im1NormalWarp, h, im2Proc.size());
  im1Reg = reversepreprocessImage(im1NormalWarp, flipflop_ref, rotate_ref);
  Rcout << "DONE: warped query image" << endl;

  // change color map
  Mat im1Combine;
  cv::addWeighted(im2Proc, 0.7, im1Warp, 0.3, 0, im1Combine);

  // return as rgb
  cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);
  cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);
}

// align images with FLANN algorithm
void alignImagesFLANN(Mat &im1, Mat &im2, Mat &im1Reg, Mat &im1Overlay, 
                      Mat &imMatches, Mat &h, Rcpp::List &keypoints,
                      const bool invert_query, const bool invert_ref,
                      const char* flipflop_query, const char* flipflop_ref,
                      const char* rotate_query, const char* rotate_ref,
                      const bool run_TPS)
{
  
  // seed
  cv::setRNGSeed(0);
  
  ////////////////////
  /// Parameters /////
  ////////////////////
  
  Parameters params;

  //////////////////////
  /// Process Images ///
  //////////////////////
  
  // Convert images to grayscale
  Mat im1Gray, im2Gray;
  cvtColor(im1, im1Gray, cv::COLOR_BGR2GRAY);
  cvtColor(im2, im2Gray, cv::COLOR_BGR2GRAY);
  
  // Process images
  Mat im1Proc, im2Proc, im1NormalProc;
  im1Proc = preprocessImage(im1Gray, invert_query, flipflop_query, rotate_query);
  im1NormalProc = preprocessImage(im1, FALSE, flipflop_query, rotate_query);
  im2Proc = preprocessImage(im2Gray, invert_ref, flipflop_ref, rotate_ref);
  
  // ////////////////////////////////////
  // /// Compute SIFT+FLANN+Homograpy ///
  // ////////////////////////////////////
  
  // RUN SIFT+FLANN+Homography
  bool is_faulty = FALSE;
  cv::Mat mask;
  std::vector<Point2f> points1, points2;
  getSIFTTransformationMatrix(im1Proc, im1NormalProc, im2Proc, im1, im2, h, mask, imMatches, 
                              points1, points2, params, is_faulty);
  
  // Use homography to warp image
  Mat im1Warp, im1NormalWarp;
  warpPerspective(im1Proc, im1Warp, h, im2Proc.size());
  warpPerspective(im1NormalProc, im1NormalWarp, h, im2Proc.size());
  Rcout << "DONE: warped query image" << endl;
  
  ///////////////////////
  /// Find Homography ///
  ///////////////////////
  
  // continue with TPS or do FLANN only
  Mat im1Reg_Warp_nonrigid;
  Mat im1Reg_NormalWarp_nonrigid;
  if(is_faulty || run_TPS){
    
    // change color map
    Mat im1Combine;
    cv::addWeighted(im2Proc, 0.7, im1Warp, 0.3, 0, im1Combine);
    
    // Reverse process
    im1Reg = reversepreprocessImage(im1NormalWarp, flipflop_ref, rotate_ref);
    
    // return as rgb
    cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);
    cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);

  // TPS is requested (only if FLANN succeeded)
  } else {
    
    // Filtered points (inliers) based on the mask
    std::vector<cv::Point2f> filtered_points1;
    std::vector<cv::Point2f> filtered_points2;
    for (int i = 0; i < mask.rows; i++) {
      if (mask.at<uchar>(i)) {
        filtered_points1.push_back(points1[i]);
        filtered_points2.push_back(points2[i]);
      }
    }
    removeCloseMatches(filtered_points1, filtered_points2);
    
    // transform query
    std::vector<cv::Point2f> filtered_points1_reg;
    cv::perspectiveTransform(filtered_points1, filtered_points1_reg, h);
    
    // get TPS matches
    std::vector<cv::DMatch> matches;
    for (unsigned int i = 0; i < filtered_points2.size(); i++)
      matches.push_back(cv::DMatch(i, i, 0));
    
    // calculate TPS transformation
    Ptr<ThinPlateSplineShapeTransformer> tps = cv::createThinPlateSplineShapeTransformer(0);
    tps->estimateTransformation(filtered_points2, filtered_points1_reg, matches);
    
    // save keypoints 
    keypoints[0] = point2fToNumericMatrix(filtered_points2);
    keypoints[1] = point2fToNumericMatrix(filtered_points1_reg); 
    
    // determine extension limits for both images
    int y_max = max(im1Warp.rows, im2.rows);
    int x_max = max(im1Warp.cols, im2.cols);
    
    // extend images
    cv::copyMakeBorder(im1Warp, im1Warp, 0.0, (int) (y_max - im1Warp.rows), 0.0, (x_max - im1Warp.cols), cv::BORDER_CONSTANT, Scalar(0, 0, 0));
    cv::copyMakeBorder(im1NormalWarp, im1NormalWarp, 0.0, (int) (y_max - im1NormalWarp.rows), 0.0, (x_max - im1NormalWarp.cols), cv::BORDER_CONSTANT, Scalar(0, 0, 0));

    // transform image
    Mat im1Reg_Warp_nonrigid;
    Mat im1Reg_NormalWarp_nonrigid;
    tps->warpImage(im1Warp, im1Reg_Warp_nonrigid);
    tps->warpImage(im1NormalWarp, im1Reg_NormalWarp_nonrigid);
    
    // resize image
    cv::Mat im1Reg_NormalWarp_nonrigid_cropped  = im1Reg_NormalWarp_nonrigid(cv::Range(0,im2Proc.size().height), cv::Range(0,im2Proc.size().width));
    im1Reg_NormalWarp_nonrigid = im1Reg_NormalWarp_nonrigid_cropped.clone();
    
    cv::Mat im1Reg_Warp_nonrigid_cropped  = im1Reg_Warp_nonrigid(cv::Range(0,im2Proc.size().height), cv::Range(0,im2Proc.size().width));
    im1Reg_Warp_nonrigid = im1Reg_Warp_nonrigid_cropped.clone();
    
    // change color map
    Mat im1Combine;
    cv::addWeighted(im2Proc, 0.7, im1Reg_Warp_nonrigid, 0.3, 0, im1Combine);

    // Reverse process
    im1Reg = reversepreprocessImage(im1Reg_NormalWarp_nonrigid, flipflop_ref, rotate_ref);
    
    // return as rgb
    cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);
    cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);
  }
}

// [[Rcpp::export]]
Rcpp::List automated_registeration_rawvector(Rcpp::RawVector ref_image, Rcpp::RawVector query_image,
                                             const int width1, const int height1,
                                             const int width2, const int height2,
                                             const float GOOD_MATCH_PERCENT, const int MAX_FEATURES,
                                             const bool invert_query, const bool invert_ref,
                                             Rcpp::String flipflop_query, Rcpp::String flipflop_ref,
                                             Rcpp::String rotate_query, Rcpp::String rotate_ref,
                                             Rcpp::String matcher, Rcpp::String method)
{
  // Return data
  Rcpp::List out(5);
  Rcpp::List out_trans(2);
  Rcpp::List keypoints(2);
  Mat imOverlay, imReg, h, imMatches;
  
  // Read reference image
  cv::Mat imReference = imageToMat(ref_image, width1, height1);

  // Read image to be aligned
  cv::Mat im = imageToMat(query_image, width2, height2);
  
  // Homography (with FLANN)
  if(strcmp(matcher.get_cstring(), "FLANN") == 0 && (strcmp(method.get_cstring(), "Homography") == 0 || strcmp(method.get_cstring(), "Homography + Non-Rigid") == 0)){
    const bool run_TPS = strcmp(method.get_cstring(), "Homography + Non-Rigid") == 0;
    Rcout << "Running SIFT+FLANN Alignment" << ((run_TPS) ? " with TPS" : "") << endl;
    alignImagesFLANN(im, imReference, imReg, imOverlay, imMatches, 
                     h, keypoints, 
                     invert_query, invert_ref,
                     flipflop_query.get_cstring(), flipflop_ref.get_cstring(), 
                     rotate_query.get_cstring(), rotate_ref.get_cstring(), run_TPS);
  }
  
  // Homography (with Brute-Force matching)
  if(strcmp(matcher.get_cstring(), "BRUTE-FORCE") == 0 && strcmp(method.get_cstring(), "Homography") == 0){
    Rcout << "Running BRUTE-FORCE Alignment" << endl;
    alignImagesBRUTE(im, imReference, imReg, imOverlay, imMatches, 
                     h, GOOD_MATCH_PERCENT, MAX_FEATURES);
  }

  // transformation matrix, can be either a matrix, set of keypoints or both
  out_trans[0] = matToNumericMatrix(h.clone());
  out_trans[1] = keypoints;
  out[0] = out_trans;
  
  // destination image, registered image, keypoint matching image
  out[1] = matToImage(imReference.clone());
  
  // registered image
  out[2] = matToImage(imReg.clone());
  
  // keypoint matching image
  out[3] = matToImage(imMatches.clone());
  
  // overlay image
  out[4] = matToImage(imOverlay.clone());
  
  // return
  return out;
}