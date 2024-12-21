#include <Rcpp.h>
#include <opencv2/opencv.hpp>
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/shape/shape_transformer.hpp"
// #include <opencv2/imgproc.hpp>

using namespace Rcpp;
using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

// Function to convert a RawVector for magick images to a cv::Mat object
cv::Mat imageToMat(Rcpp::RawVector image_data, int width, int height) {

  // Create cv::Mat object
  cv::Mat mat(height, width, CV_8UC3, image_data.begin());

  // Convert from RGBA to BGRA
  cv::cvtColor(mat, mat, cv::COLOR_RGBA2BGR);

  return mat;
}

// Function to convert a cv::Mat object to a RawVector for magick images
Rcpp::RawVector matToImage(cv::Mat mat) {

  // Create RawVector object
  Rcpp::RawVector rawvec(mat.total() * mat.elemSize());
  rawvec.attr("dim") = Rcpp::Dimension(3, mat.cols, mat.rows);

  // Copy Mat data to RawVector
  std::memcpy(rawvec.begin(), mat.data, rawvec.size());

  return rawvec;
}

// Function to convert a NumericMatrix object to a cv::Mat
cv::Mat numericMatrixToMat(Rcpp::NumericMatrix nm) {
  cv::Mat m(nm.rows(), nm.cols(), CV_64F);
  for (int i = 0; i < nm.rows(); ++i) {
    for (int j = 0; j < nm.cols(); ++j) {
      m.at<double>(i, j) = nm(i, j);
    }
  }
  return m;
}

// Function to convert a NumericMatrix object to a cv::Point2f
std::vector<cv::Point2f> numericMatrixToPoint2f(Rcpp::NumericMatrix mat) {
  std::vector<cv::Point2f> points;
  for (int i = 0; i < mat.nrow(); i++) {
    points.push_back(cv::Point2f(mat(i, 0), mat(i, 1)));
  }
  return points;
}

// Function to convert a cv::Mat object to a NumericMatrix
Rcpp::NumericMatrix matToNumericMatrix(cv::Mat m) {
  Rcpp::NumericMatrix nm(m.rows, m.cols);
  for (int i = 0; i < m.rows; ++i) {
    for (int j = 0; j < m.cols; ++j) {
      nm(i, j) = m.at<double>(i, j);
    }
  }
  return nm;
}

// Function to convert a cv::Point2f object to a NumericMatrix
Rcpp::NumericMatrix point2fToNumericMatrix(std::vector<cv::Point2f> points) {
  int n = points.size();
  Rcpp::NumericMatrix mat(n, 2);
  for (int i = 0; i < n; i++) {
    mat(i, 0) = points[i].x;
    mat(i, 1) = points[i].y;
  }
  return mat;
}

// Function to convert a cv::Keypoint object to a std::vector<double>
std::vector<double> KeyPointToDoubleVector(std::vector<cv::KeyPoint> points) {
  int n = points.size();
  std::vector<double> vec(n);
  for (int i = 0; i < n; i++) {
    vec[i] = (double) points[i].pt.x;
  }
  return vec;
}

// Function to convert a cv::Point2f object to a std::vector<double>
std::vector<double> Point2fToDoubleVector(std::vector<cv::Point2f> points) {
  int n = points.size();
  std::vector<double> vec(n);
  for (int i = 0; i < n; i++) {
    vec[i] = (double) points[i].x;
  }
  return vec;
}

// Function to convert a cv::Point2f object to a cv::Mat
std::vector<cv::Point2f> matToPoint2f(cv::Mat mat) {
  std::vector<cv::Point2f> points;

  // Assuming the matrix has 2 columns (x and y coordinates)
  if (mat.cols != 2) {
    // cerr << "Input matrix must have exactly 2 columns for x and y coordinates." << endl;
    return points;
  }

  // Iterate over the rows of the matrix
  for (int i = 0; i < mat.rows; ++i) {
    // Extract x and y coordinates from the matrix
    float x = mat.at<float>(i, 0);
    float y = mat.at<float>(i, 1);

    // Create Point2f object and add it to the vector
    points.push_back(Point2f(x, y));
  }

  return points;
}

cv::Mat point2fToMat(std::vector<cv::Point2f> points) {
  cv::Mat mat(points.size(), 2, CV_32F);

  // Iterate over the vector of Point2f
  for (size_t i = 0; i < points.size(); ++i) {

    // Assign x and y coordinates to the matrix
    mat.at<float>(i, 0) = points[i].x;
    mat.at<float>(i, 1) = points[i].y;
  }

  return mat;
}

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
    // cout << "Image is inverted!" << endl;
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

// calculate standard deviation of a vector
double cppSD(std::vector<cv::KeyPoint> points)
{
  std::vector<double> inVec = KeyPointToDoubleVector(points);
  int n = inVec.size();
  double sum = std::accumulate(inVec.begin(), inVec.end(), 0.0);
  double mean = sum / inVec.size();

  for(std::vector<double>::iterator iter = inVec.begin();
      iter != inVec.end(); ++iter){
    double temp;
    temp= (*iter - mean)*(*iter - mean);
    *iter = temp;
  }

  double sd = std::accumulate(inVec.begin(), inVec.end(), 0.0);
  return std::sqrt( sd / (n-1) );
}

double cppSD(std::vector<cv::Point2f> points)
{
  std::vector<double> inVec = Point2fToDoubleVector(points);
  int n = inVec.size();
  double sum = std::accumulate(inVec.begin(), inVec.end(), 0.0);
  double mean = sum / inVec.size();
  
  for(std::vector<double>::iterator iter = inVec.begin();
      iter != inVec.end(); ++iter){
    double temp;
    temp= (*iter - mean)*(*iter - mean);
    *iter = temp;
  }
  
  double sd = std::accumulate(inVec.begin(), inVec.end(), 0.0);
  return std::sqrt( sd / (n-1) );
}


// check if keypoints are degenerate
std::string check_degenerate(std::vector<cv::KeyPoint> keypoints1, std::vector<cv::KeyPoint> keypoints2) {

  // get sd
  double keypoints1_sd = cppSD(keypoints1);
  double keypoints2_sd = cppSD(keypoints2);
  
  // get warning message
  std::string message;
  if(keypoints1_sd < 1.0 | keypoints2_sd < 1.0){
    message = "degenarate";
    // cout << "WARNING: points may be in a degenerate configuration." << endl;
  } else {
    message = "not degenarate";
  }

  return message;
}

// // check transformation of points
// std::string check_transformation_by_pts_mean_sqrt(std::vector<cv::KeyPoint> keypoints1, std::vector<cv::KeyPoint> keypoints2, Mat &h, Mat &mask){
//   
//   // perspective transformation of query keypoints
//   std::vector<cv::Point2f> keypoints1_reg;
//   cv::perspectiveTransform(keypoints1, keypoints1_reg, h);
//   
// }

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
    // cout << "WARNING: Transformation may be poor - transformed points grid seem to be concentrated!" << endl;
  } else {
    message = "small distribution";
  }

  return message;
}


// do overall checks on keypoints and images
std::string check_transformation_metrics(std::vector<cv::KeyPoint> keypoints1, std::vector<cv::KeyPoint> keypoints2, Mat &im1, Mat &im2, Mat &h, Mat &mask) {
  
  // check keypoint standard deviation
  std::string degenerate;
  degenerate = check_degenerate(keypoints1, keypoints2);
  
  // //  check transformation
  // std::string transformation;
  // transformation = check_transformation_by_pts_mean_sqrt(keypoints1, keypoints2, h, mask);
  
  //  check distribution
  std::string distribution;
  distribution = check_transformation_by_point_distribution(im2, h);
    
  return degenerate;
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
  // cout << "DONE: orb based key-points detection and descriptors computation" << endl;

  // Match features.
  std::vector<DMatch> matches;
  Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("BruteForce-Hamming");
  matcher->match(descriptors1, descriptors2, matches, Mat());
  // cout << "DONE: BruteForce-Hamming - descriptor matching" << endl;

  // Sort matches by score
  std::sort(matches.begin(), matches.end());

  // Remove not so good matches
  const int numGoodMatches = matches.size() * GOOD_MATCH_PERCENT;
  matches.erase(matches.begin()+numGoodMatches, matches.end());
  // cout << "DONE: get good matches by distance thresholding" << endl;

  // Extract location of good matches
  std::vector<Point2f> points1, points2;
  for( size_t i = 0; i < matches.size(); i++ )
  {
    points1.push_back( keypoints1[ matches[i].queryIdx ].pt );
    points2.push_back( keypoints2[ matches[i].trainIdx ].pt );
  }

  // Find homography
  h = findHomography(points1, points2, RANSAC);
  // cout << "DONE: calculated homography matrix" << endl;

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
  // im1Reg = im1Overlay;
  // cout << "DONE: warped query image" << endl;

  // change color map
  // Mat im1Combine;
  // cv::addWeighted(im2Gray, 0.7, im1Reg, 0.3, 0, im1Overlay);
  //
  // // return as rgb
  // cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);

  // cv::imwrite("dest.jpg", im2);
  // cv::imwrite("source.jpg", im1Overlay);
}

// align images with FLANN algorithm
void alignImagesFLANN(Mat &im1, Mat &im2, Mat &im1Reg, Mat &im1Overlay, Mat &imMatches, Mat &h,
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
  // cout << "DONE: sift based key-points detection and descriptors computation" << endl;

  // Match features using FLANN matching
  std::vector<std::vector<DMatch>> matches;
  cv::FlannBasedMatcher custom_matcher = cv::FlannBasedMatcher(cv::makePtr<cv::flann::KDTreeIndexParams>(5), cv::makePtr<cv::flann::SearchParams>(50, 0, TRUE));
  cv::Ptr<cv::FlannBasedMatcher> matcher = custom_matcher.create();
  matcher->knnMatch(descriptors1, descriptors2, matches, 2);
  // cout << "DONE: FLANN - Fast Library for Approximate Nearest Neighbors - descriptor matching" << endl;

  // Find good matches
  // goodMatches = get_good_matches(matches)
  std::vector<DMatch> good_matches;
  getGoodMatches(matches, good_matches);
  // cout << "DONE: get good matches by distance thresholding" << endl;

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
  // cout << "DONE: calculated homography matrix" << endl;

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
  
  // Visualize matches, use for demonstration purposes
  // std::vector<cv::DMatch> top_matches_vis;
  // std::vector<cv::KeyPoint> keypoints1_best_vis, keypoints2_best_vis;
  // for(size_t i = 0; i < good_matches.size(); i++)
  // {
  //   keypoints1_best_vis.push_back(keypoints1[good_matches[i].queryIdx]);
  //   keypoints2_best_vis.push_back(keypoints2[good_matches[i].trainIdx]);
  //   top_matches_vis.push_back(cv::DMatch(static_cast<int>(i), static_cast<int>(i), 0));
  // }
  // Mat im1Proc_vis, im2Proc_vis, im1NormalProc_vis, imMatches_vis;
  // im1Proc_vis = preprocessImage(im1, invert_query, flipflop_query, rotate_query);
  // im2Proc_vis = preprocessImage(im2, invert_ref, flipflop_ref, rotate_ref);
  // drawMatches(im1Proc_vis, keypoints1_best, im2Proc_vis, keypoints2_best, top_matches, imMatches_vis);
  // cv::imwrite("matches.jpg", imMatches_vis);
  
  // check keypoints
  std::string is_faulty = check_transformation_metrics(keypoints1_best, keypoints2_best, im1, im2, h, mask);

  // Use homography to warp image
  Mat im1Warp, im1NormalWarp;
  warpPerspective(im1Proc, im1Warp, h, im2Proc.size());
  warpPerspective(im1NormalProc, im1NormalWarp, h, im2Proc.size());
  im1Reg = reversepreprocessImage(im1NormalWarp, flipflop_ref, rotate_ref);
  // cout << "DONE: warped query image" << endl;

  // // overlay image
  // Mat im1ColorMap;
  // cv::applyColorMap(im1Warp, im1ColorMap, cv::COLORMAP_HOT);

  // change color map
  Mat im1Combine;
  cv::addWeighted(im2Proc, 0.7, im1Warp, 0.3, 0, im1Combine);

  // return as rgb
  // cvtColor(im1Combine, im1Reg, cv::COLOR_GRAY2BGR);
  cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);
  cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);
}

// align images with FLANN algorithm
void alignImagesFLANNTPS(Mat &im1, Mat &im2, Mat &im1Reg, Mat &im1Overlay, 
                         Mat &imMatches, Mat &h, Rcpp::List &keypoints,
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
  // cout << "DONE: sift based key-points detection and descriptors computation" << endl;
  
  // Match features using FLANN matching
  std::vector<std::vector<DMatch>> matches;
  cv::FlannBasedMatcher custom_matcher = cv::FlannBasedMatcher(cv::makePtr<cv::flann::KDTreeIndexParams>(5), cv::makePtr<cv::flann::SearchParams>(50, 0, TRUE));
  cv::Ptr<cv::FlannBasedMatcher> matcher = custom_matcher.create();
  matcher->knnMatch(descriptors1, descriptors2, matches, 2);
  // cout << "DONE: FLANN - Fast Library for Approximate Nearest Neighbors - descriptor matching" << endl;
  
  // Find good matches
  std::vector<DMatch> good_matches;
  getGoodMatches(matches, good_matches);
  // cout << "DONE: get good matches by distance thresholding" << endl;
  
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
  // cout << "DONE: calculated homography matrix" << endl;

  // Use homography to warp image
  Mat im1Warp, im1NormalWarp;
  warpPerspective(im1Proc, im1Warp, h, im2Proc.size());
  warpPerspective(im1NormalProc, im1NormalWarp, h, im2Proc.size());
  // im1Reg = reversepreprocessImage(im1NormalWarp, flipflop_ref, rotate_ref);
  // cout << "DONE: warped query image" << endl;
  
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
  std::string is_faulty = check_transformation_metrics(keypoints1_best, keypoints2_best, im1, im2, h, mask);
  
  // continue with TPS or do FLANN only
  Mat im1Reg_Warp_nonrigid;
  Mat im1Reg_NormalWarp_nonrigid;
  if(is_faulty == "degenerate"){
    
    // change color map
    Mat im1Combine;
    cv::addWeighted(im2Proc, 0.7, im1Warp, 0.3, 0, im1Combine);
    
    // Reverse process
    im1Reg = reversepreprocessImage(im1NormalWarp, flipflop_ref, rotate_ref);
    
    // return as rgb
    cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);
    cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);

  // if FLANN succeeds   
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
  // auto tps = cv::createThinPlateSplineShapeTransformer();
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
  
  // calculate homography transformation
  Mat im1Affine;
  // h = estimateAffine2D(query_mat, ref_mat);
  h = findHomography(query_mat, ref_mat);
  // std::cout << h.size() << "\n";
  // std::cout << h << "\n";
  cv::warpPerspective(im1, im1Affine, h, im2.size());
  // cv::warpAffine(im1, im1Affine, h, im2.size());
  std::vector<cv::Point2f> query_reg;
  cv::perspectiveTransform(query_mat, query_reg, h);
  
  // calculate TPS transformation
  // auto tps = cv::createThinPlateSplineShapeTransformer();
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
  // tps->warpImage(im1, im1Reg,  cv::INTER_LINEAR, cv::WARP_FILL_OUTLIERS);
  tps->warpImage(im1Affine, im1Reg);
  
  // resize image
  // cv::resize(im1Reg, im1Reg, im2.size());
  cv::Mat im1Reg_cropped  = im1Reg(cv::Range(0,im2.size().height), cv::Range(0,im2.size().width));
  im1Reg = im1Reg_cropped.clone();
}

// [[Rcpp::export]]
Rcpp::NumericMatrix applyTransform(Rcpp::NumericMatrix coords, Rcpp::List mapping)
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
      cv::warpPerspective(im, im_temp, h, imReference.size()); 
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
      // cv::resize(im1Reg, im1Reg, im2.size());
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
  
  // Homography (with FLANN) + Non-rigid (TPS)
  if(strcmp(matcher.get_cstring(), "FLANN") == 0 && strcmp(method.get_cstring(), "Homography + Non-Rigid") == 0){
    alignImagesFLANNTPS(im, imReference, imReg, imOverlay, imMatches, 
                        h, keypoints, 
                        invert_query, invert_ref,
                        flipflop_query.get_cstring(), flipflop_ref.get_cstring(), 
                        rotate_query.get_cstring(), rotate_ref.get_cstring());
  }
  
  // Homography (with FLANN)
  if(strcmp(matcher.get_cstring(), "FLANN") == 0 && strcmp(method.get_cstring(), "Homography") == 0){
  alignImagesFLANN(im, imReference, imReg, imOverlay, imMatches, 
                   h, invert_query, invert_ref,
                   flipflop_query.get_cstring(), flipflop_ref.get_cstring(), 
                   rotate_query.get_cstring(), rotate_ref.get_cstring());
  }
  
  // Homography (with Brute-Force matching)
  if(strcmp(matcher.get_cstring(), "BRUTE-FORCE") == 0 && strcmp(method.get_cstring(), "Homography") == 0){
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
  if(strcmp(method.get_cstring(), "Homography + Non-Rigid") == 0){
    alignImagesAffineTPS(im, imReference, imReg, 
                         h, keypoints,
                         query_landmark, reference_landmark);
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