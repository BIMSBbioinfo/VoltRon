#include <Rcpp.h>
#include <opencv2/opencv.hpp>
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/features2d.hpp"

using namespace Rcpp;
using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

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

// void computeSiftTiles(Mat image, std::vector<KeyPoint> &keypoints, Mat &descriptors, const int tile_size_tuple[2], const int tile_overlap, Ptr<Feature2D> sift){
//
//   // image shape
//   int height = image.rows;
//   int width = image.cols;
//
//   // Extend the image so that it can be divided into equal size tiles
//   float h_extend = ((height - 1) / (tile_size_tuple[1] - tile_overlap) + 1) * (tile_size_tuple[1] - tile_overlap) + tile_overlap;
//   float w_extend = ((width - 1) / (tile_size_tuple[0] - tile_overlap) + 1) * (tile_size_tuple[0] - tile_overlap) + tile_overlap;
//   cv::Mat image_extended = cv::Mat::zeros(height, width, image.type());
//
// }

void getGoodMatches(std::vector<std::vector<DMatch>> matches, std::vector<DMatch> &good_matches, const float lowe_ratio = 0.8)
{
  for (size_t i = 0; i < matches.size(); i++)
  {
    if (matches[i][0].distance < lowe_ratio * matches[i][1].distance)
    {
      good_matches.push_back(matches[i][0]);
    }
  }
}

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
    cout << "Image is inverted!" << endl;
  } else {
    imProcess = imFlipFlop;
  }

  // return
  return imProcess;
}

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
  cout << "DONE: orb based key-points detection and descriptors computation" << endl;

  // Match features.
  std::vector<DMatch> matches;
  Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("BruteForce-Hamming");
  matcher->match(descriptors1, descriptors2, matches, Mat());
  cout << "DONE: BruteForce-Hamming - descriptor matching" << endl;

  // Sort matches by score
  std::sort(matches.begin(), matches.end());

  // Remove not so good matches
  const int numGoodMatches = matches.size() * GOOD_MATCH_PERCENT;
  matches.erase(matches.begin()+numGoodMatches, matches.end());
  cout << "DONE: get good matches by distance thresholding" << endl;

  // Extract location of good matches
  std::vector<Point2f> points1, points2;
  for( size_t i = 0; i < matches.size(); i++ )
  {
    points1.push_back( keypoints1[ matches[i].queryIdx ].pt );
    points2.push_back( keypoints2[ matches[i].trainIdx ].pt );
  }

  // Find homography
  h = findHomography(points1, points2, RANSAC);
  cout << "DONE: calculated homography matrix" << endl;

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
  im1Reg = im1Overlay;
  cout << "DONE: warped query image" << endl;
}

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
  cout << "DONE: sift based key-points detection and descriptors computation" << endl;

  // Match features using FLANN matching
  std::vector<std::vector<DMatch>> matches;
  cv::FlannBasedMatcher custom_matcher = cv::FlannBasedMatcher(cv::makePtr<cv::flann::KDTreeIndexParams>(5), cv::makePtr<cv::flann::SearchParams>(50, 0, TRUE));
  cv::Ptr<cv::FlannBasedMatcher> matcher = custom_matcher.create();
  matcher->knnMatch(descriptors1, descriptors2, matches, 2);
  cout << "DONE: FLANN - Fast Library for Approximate Nearest Neighbors - descriptor matching" << endl;

  // Find good matches
  // goodMatches = get_good_matches(matches)
  std::vector<DMatch> good_matches;
  getGoodMatches(matches, good_matches);
  cout << "DONE: get good matches by distance thresholding" << endl;

  // Extract location of good matches
  std::vector<Point2f> points1, points2;
  for( size_t i = 0; i < good_matches.size(); i++ )
  {
    points1.push_back(keypoints1[good_matches[i].queryIdx].pt);
    points2.push_back(keypoints2[good_matches[i].trainIdx].pt);
  }

  // Find homography
  h = findHomography(points1, points2, RANSAC, 5);
  cout << "DONE: calculated homography matrix" << endl;

  // Draw top matches and good ones only
  std::vector<cv::DMatch> top_matches;
  std::vector<KeyPoint> keypoints1_best, keypoints2_best;
  for( size_t i = 0; i < good_matches.size(); i++ )
  {
    keypoints1_best.push_back(keypoints1[good_matches[i].queryIdx]);
    keypoints2_best.push_back(keypoints2[good_matches[i].trainIdx]);
    top_matches.push_back(cv::DMatch(static_cast<int>(i), static_cast<int>(i), 0));
  }
  drawMatches(im1Proc, keypoints1_best, im2Proc, keypoints2_best, top_matches, imMatches);

  // Use homography to warp image
  Mat im1Warp, im1NormalWarp;
  warpPerspective(im1Proc, im1Warp, h, im2Proc.size());
  warpPerspective(im1NormalProc, im1NormalWarp, h, im2Proc.size());
  im1Reg = reversepreprocessImage(im1NormalWarp, flipflop_ref, rotate_ref);
  cout << "DONE: warped query image" << endl;

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

  //destination and overlay
  cv::imwrite("dest.jpg", im2);
  cv::imwrite("source.jpg", im1Overlay);
}

// [[Rcpp::export]]
Rcpp::List automated_registeration_rawvector(Rcpp::RawVector ref_image, Rcpp::RawVector query_image,
                                             const int width1, const int height1,
                                             const int width2, const int height2,
                                             const float GOOD_MATCH_PERCENT, const int MAX_FEATURES,
                                             const bool invert_query, const bool invert_ref,
                                             Rcpp::String flipflop_query, Rcpp::String flipflop_ref,
                                             Rcpp::String rotate_query, Rcpp::String rotate_ref,
                                             Rcpp::String method)
{
  // define return data, 1 = transformation matrix, 2 = aligned image
  Rcpp::List out(5);

  // Read reference image
  cv::Mat imReference = imageToMat(ref_image, width1, height1);

  // Read image to be aligned
  cv::Mat im = imageToMat(query_image, width2, height2);

  // Registered image will be resotred in imReg.
  // The estimated homography will be stored in h.
  // The matching illustration of both images with be given in imMatches.
  Mat imOverlay, imReg, h, imMatches;

  // Align images
  if(strcmp(method.get_cstring(), "FLANN") == 0){
    cout << "Fast Library for Approximate Nearest Neighbors (FLANN) - descriptor matching" << endl;
    alignImagesFLANN(im, imReference, imReg, imOverlay, imMatches, h, invert_query, invert_ref,
                     flipflop_query.get_cstring(), flipflop_ref.get_cstring(), rotate_query.get_cstring(), rotate_ref.get_cstring());
  }
  if(strcmp(method.get_cstring(), "BRUTE-FORCE") == 0){
    cout << "BruteForce-Hamming - descriptor matching" << endl;
    alignImagesBRUTE(im, imReference, imReg, imOverlay, imMatches, h, GOOD_MATCH_PERCENT, MAX_FEATURES);
  }

  // return transformation matrix, destinated image, registered image, and keypoint matching image
  out[0] = matToNumericMatrix(h.clone());
  out[1] = matToImage(imReference.clone());
  out[2] = matToImage(imReg.clone());
  out[3] = matToImage(imMatches.clone());
  out[4] = matToImage(imOverlay.clone());
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix perspectiveTransform(Rcpp::NumericMatrix coords, Rcpp::NumericMatrix hmatrix)
{
  // Get coordinates as cv::Mat
  std::vector<cv::Point2f> coords_mat = numericMatrixToPoint2f(coords);
  std::vector<cv::Point2f> coords_reg;
  cv::Mat hmatrix_mat = numericMatrixToMat(hmatrix);

  // transform coordinates
  cv::perspectiveTransform(coords_mat, coords_reg, hmatrix_mat);

  // return transformation matrix and alignment images
  Rcpp::NumericMatrix coords_regToMat = point2fToNumericMatrix(coords_reg);

  // return
  return coords_regToMat;
}
