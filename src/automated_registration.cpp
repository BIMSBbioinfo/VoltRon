#include <Rcpp.h>
#include <sys/resource.h>

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

// SIFT Parameters
struct SIFTParameters
{ 
  const int sift_nfeatures=20000;
  const int max_nfeatures=10000;
  const int tile_size=5000;
  const int tile_overlap=50;
  const int ransac_pixel_threshold = 5;
  const float ransac_confidence = 0.95;
  const int ransac_maxIters=2000;
};

// memory check
void log_mem_usage(const std::string& label = "") {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  long rss_b = usage.ru_maxrss;
  
  double rss_kb = rss_b / 1024.0;
  double rss_mb = rss_kb / 1024.0;
  double rss_gb = rss_mb / 1024.0;
  
  Rcpp::Rcout << "Used Memory [" << label << "]: " << rss_gb << " GB" << std::endl;
}

double object_size(long bsize) {
  
  double rss_kb = bsize / 1024.0;
  double rss_mb = rss_kb / 1024.0;
  double rss_gb = rss_mb / 1024.0;
  
  return rss_gb;
}


// check if keypoints are degenerate
bool check_degenerate(std::vector<cv::Point2f> &points1, std::vector<cv::Point2f> &points2) {

  // get sd
  double points1_sd = cppSD(points1);
  double points2_sd = cppSD(points2);
  
  // get warning message
  bool is_degenerate = FALSE;
  if(points1_sd < 1.0 | points2_sd < 1.0){
    is_degenerate = TRUE;
    Rcout << "WARNING: points may be in a degenerate configuration." << endl;
  } 

  return is_degenerate;
}

// check distribution of registered points
std::string check_transformation_by_point_distribution(Mat &im, Mat &h){

  // message
  std::string message;
  
  // get image shape
  int height = im.rows;
  int width = im.cols;
  int height_interval = height > 50 ? (double) height/50.0 : 1;
  int width_interval = width > 50 ? (double) width/50.0 : 1;

  // perspective transformation of grid points
  std::vector<cv::Point2f> gridpoints;
  for (double i = 0.0; i <= height; i += height_interval) {
    for (double j = 0.0; j <= width; j += width_interval) {
      gridpoints.push_back(cv::Point2f(j,i));
    }
  }

  // register grid points
  std::vector<cv::Point2f> gridpoints_reg;
  if (h.rows == 2){
    cv::transform(gridpoints, gridpoints_reg, h);
  } else if(h.rows == 3) {
    cv::perspectiveTransform(gridpoints, gridpoints_reg, h);
  } else {
    message = "no distribution";
    return message;
  }
  // std::vector<cv::Point2f>().swap(gridpoints);

  // Compute the standard deviation of the transformed points
  double gridpoints_reg_sd = cppSD(gridpoints_reg);
  // std::vector<cv::Point2f>().swap(gridpoints_reg);
  
  // get warning message
  if(gridpoints_reg_sd < 1.0 | gridpoints_reg_sd > max(height, width)){
    message = "large distribution";
    Rcout << "WARNING: Transformation may be poor - transformed points grid seem to be concentrated!" << endl;
  } else {
    message = "small distribution";
  }

  return message;
}

bool check_matches(Mat &mask){
  int j=0;
  for (int i = 0; i < mask.rows; i++) {
    if (mask.at<uchar>(i)) {
      j++;
    }
  }
  return j > 6;
}

// do overall checks on keypoints and images
bool check_transformation_metrics(std::vector<cv::Point2f> &points1, std::vector<cv::Point2f> &points2, Mat &im2, Mat &h, Mat &mask) {
  
  // check keypoint standard deviation
  bool is_degenerate = check_degenerate(points1, points2);

  // TODO: check transformation
  // make keypoints from points
  // std::string transformation;
  // transformation = check_transformation_by_pts_mean_sqrt(keypoints1, keypoints2, h, mask);

  //  check distribution
  std::string distribution;
  distribution = check_transformation_by_point_distribution(im2, h);

  return is_degenerate;
}

// get good matching keypoints
void getGoodMatches(std::vector<std::vector<DMatch>> &matches12,std::vector<std::vector<DMatch>> &matches21,
                    std::vector<DMatch> &good_matches, const float lowe_ratio = 0.8)
{
  // direction wise good matches
  std::vector<DMatch> good_matches12, good_matches21;
  for (size_t i = 0; i < matches12.size(); i++) {
    if (matches12[i][0].distance < lowe_ratio * matches12[i][1].distance) {
      good_matches12.push_back(matches12[i][0]);
    }
  }
  for (size_t i = 0; i < matches21.size(); i++) {
    if (matches21[i][0].distance < lowe_ratio * matches21[i][1].distance) {
      good_matches21.push_back(matches21[i][0]);
    }
  }

  // get good matches as dictionaries
  std::map<std::pair<int, int>, std::vector<DMatch>> matches12_map;
  std::map<std::pair<int, int>, std::vector<DMatch>> matches21_map;
  for (const auto& match : good_matches12) {
    auto rounded_pt = std::make_pair(match.queryIdx, match.trainIdx);
    matches12_map[rounded_pt].push_back(match);
  }
  for (const auto& match : good_matches21) {
    auto rounded_pt = std::make_pair(match.trainIdx, match.queryIdx);
    matches21_map[rounded_pt].push_back(match);
  }

  // get mutual matches
  for (const auto& item : matches12_map) {
    auto query_train_pair = item.first;
    const cv::DMatch& m12 = item.second[0];

    auto it = matches21_map.find(query_train_pair);
    if (it != matches21_map.end()) {
      const cv::DMatch& m21 = it->second[0];

      // Calculate average distance
      float avg_distance = (m12.distance + m21.distance) / 2.0f;

      // Create a new match object with averaged distance
      cv::DMatch mutual_match(m12.queryIdx, m12.trainIdx, avg_distance);
      good_matches.push_back(mutual_match);
    }
  }

  // TODO: can I release there now ?
  // std::vector<DMatch>().swap(good_matches12);
  // std::vector<DMatch>().swap(good_matches21);
  // std::map<std::pair<int, int>, std::vector<DMatch>>().swap(matches12_map);
  // std::map<std::pair<int, int>, std::vector<DMatch>>().swap(matches21_map);
}

void getGoodMatches_temp(std::vector<std::vector<DMatch>> matches, std::vector<DMatch> &good_matches, const float lowe_ratio = 0.8)
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
  // std::vector<Point2f>().swap(points1);
  // std::vector<Point2f>().swap(points2);
  // points1 = std::move(filtered_points1);
  // points2 = std::move(filtered_points2);
  points1 = filtered_points1;
  points2 = filtered_points2;
}

// get good matching keypoints
void getFLANNMatches(cv::Mat &descriptors1, 
                     cv::Mat &descriptors2,
                     std::vector<std::vector<DMatch>> &matches12, 
                     std::vector<std::vector<DMatch>> &matches21)
{
  cv::FlannBasedMatcher custom_matcher = cv::FlannBasedMatcher(cv::makePtr<cv::flann::KDTreeIndexParams>(5), 
                                                               cv::makePtr<cv::flann::SearchParams>(50, 0, TRUE));
  cv::Ptr<cv::FlannBasedMatcher> matcher = custom_matcher.create();
  matcher->knnMatch(descriptors1, descriptors2, matches12, 2);
  matcher->knnMatch(descriptors2, descriptors1, matches21, 2);
}

void filterDuplicateKeypoints(std::vector<cv::KeyPoint> &keypoints, cv::Mat &descriptors) {
  
  // Map to group keypoints by their location (rounded for integer coordinates)
  std::map<std::pair<int, int>, std::vector<int>> point_to_indices;
  
  for (size_t i = 0; i < keypoints.size(); ++i) {
    auto rounded_pt = std::make_pair(static_cast<int>(keypoints[i].pt.x), static_cast<int>(keypoints[i].pt.y));
    point_to_indices[rounded_pt].push_back(i);
  }
  
  std::vector<int> unique_indices;
  
  for (const auto& [_, indices] : point_to_indices) {
    if (indices.size() > 1) {
      // Find the index with the highest response among duplicates
      int best_idx = *std::max_element(indices.begin(), indices.end(), [&](int a, int b) {
        return keypoints[a].response < keypoints[b].response;
      });
      unique_indices.push_back(best_idx);
    } else {
      unique_indices.push_back(indices[0]);
    }
  }
  
  // Create new keypoints and descriptors with unique indices
  std::vector<cv::KeyPoint> unique_keypoints;
  cv::Mat unique_descriptors;
  
  for (int idx : unique_indices) {
    unique_keypoints.push_back(keypoints[idx]);
    // unique_descriptors.push_back(descriptors.row(idx));
    unique_descriptors.push_back(descriptors.row(idx).clone());
  }
  
  keypoints = std::move(unique_keypoints);
  descriptors = std::move(unique_descriptors);
}

void keepTopKeypoints(std::vector<KeyPoint> &keypoints, Mat &descriptors, SIFTParameters params){
  
  // Sort keypoints and descriptors by response
  std::vector<std::pair<cv::KeyPoint, cv::Mat>> kp_desc_pairs;
  for (size_t i = 0; i < keypoints.size(); ++i) {
    kp_desc_pairs.emplace_back(keypoints[i], descriptors.row(i).clone());
  }
  
  std::sort(kp_desc_pairs.begin(), kp_desc_pairs.end(), [](const auto& a, const auto& b) {
    return a.first.response > b.first.response;
  });
  
  // Keep only the top 100000 keypoints and descriptors
  if (kp_desc_pairs.size() > params.max_nfeatures) {
    kp_desc_pairs.resize(params.max_nfeatures);
  }
  
  // TODO: we had problems here, were getting heap memory allocation
  // issues
  keypoints.clear();
  // descriptors.release();
  // descriptors = cv::Mat();
  descriptors.create(0, descriptors.cols, descriptors.type());

  for (const auto& pair : kp_desc_pairs) {
    keypoints.push_back(pair.first);
    descriptors.push_back(pair.second);
  }
}

void computeSIFTTiles(Mat &im, std::vector<KeyPoint> &keypoints, Mat &descriptors, Ptr<Feature2D> &sift,
                      SIFTParameters params){
  
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
            // tile_descriptors.release();
          }
        }
      }
    }
  }
}

bool getSIFTTransformationMatrixSingle(
    Mat &im1Proc, Mat &im2Proc, Mat &h, Mat &mask, 
    Mat &imMatches, 
    std::vector<Point2f> &points1, std::vector<Point2f> &points2, 
    // std::vector<DMatch> &good_matches, 
    const bool &run_Affine, SIFTParameters params, bool &is_faulty){
  
  //////////////////////
  /// Compute SIFT /////
  //////////////////////
  
  // Variables to store keypoints and descriptors
  std::vector<KeyPoint> keypoints1, keypoints2;
  Mat descriptors1, descriptors2;
  
  // Detect SIFT features
  Ptr<Feature2D> sift = cv::SIFT::create(params.sift_nfeatures);
  computeSIFTTiles(im1Proc, keypoints1, descriptors1, sift, params);
  computeSIFTTiles(im2Proc, keypoints2, descriptors2, sift, params);
  Rcout << "MESSAGE: Generated " << keypoints1.size() << " and " << keypoints2.size() << " keypoints"  << endl;
  Rcout << "DONE: SIFT based key-points detection and descriptors computation" << endl;
  log_mem_usage("sift landmarks");
  
  // filter duplicates
  filterDuplicateKeypoints(keypoints1, descriptors1);
  filterDuplicateKeypoints(keypoints2, descriptors2);
  
  // get top key points
  keepTopKeypoints(keypoints1, descriptors1, params);
  keepTopKeypoints(keypoints2, descriptors2, params);
  Rcout << "MESSAGE: After filtering " << keypoints1.size() << " and " << keypoints2.size() << " keypoints" << endl;
  log_mem_usage("filter landmarks");
  
  ///////////////////////
  /// Compute FLANN /////
  ///////////////////////
  
  // Match features using FLANN matching
  std::vector<std::vector<DMatch>> matches12, matches21;
  getFLANNMatches(descriptors1, descriptors2, matches12, matches21);
  Rcout << "DONE: FLANN - Fast Library for Approximate Nearest Neighbors - descriptor matching" << endl;
  log_mem_usage("flann matching");
  
  // TODO: can I release there now ?
  descriptors1.release();
  descriptors2.release();

  // Find good matches
  std::vector<DMatch> good_matches;
  getGoodMatches(matches12, matches21, good_matches);
  Rcout << "DONE: get good mutual matches by distance thresholding" << endl;
  log_mem_usage("good matches");
  
  // TODO: can I release there now ? 
  std::vector<std::vector<DMatch>>().swap(matches12);
  std::vector<std::vector<DMatch>>().swap(matches21);
  
  ///////////////////////
  /// Find Homography ///
  ///////////////////////
  
  // Extract location of good matches
  for( size_t i = 0; i < good_matches.size(); i++ )
  {
    points1.push_back(keypoints1[good_matches[i].queryIdx].pt);
    points2.push_back(keypoints2[good_matches[i].trainIdx].pt);
  }

  // check variable
  Rcout << "MESSAGE: Calculating" << (run_Affine ? " (Affine) " : " (Homography) ") << "Transformation Matrix" << endl;

  // Find transformation matrix
  if(points1.size() > 0){
    if(run_Affine){
      std::vector<uint8_t> match_mask;
      h = estimateAffine2D(points1,
                           points2,
                           match_mask,
                           cv::RANSAC,
                           params.ransac_pixel_threshold,
                           params.ransac_maxIters,
                           params.ransac_confidence);
      mask = IntVectorToMat(match_mask);
    } else {
      h = findHomography(points1,
                         points2,
                         cv::RANSAC,
                         params.ransac_pixel_threshold,
                         mask,
                         params.ransac_maxIters,
                         params.ransac_confidence);
    } 
  } else {
    Rcout <<  "WARNING: Found no matches!" << endl;
    return false;
  }
  log_mem_usage("find transformation");
  
  // Draw top matches and good ones only
  std::vector<cv::DMatch> top_matches;
  std::vector<cv::KeyPoint> keypoints1_best, keypoints2_best;
  for(size_t i = 0; i < good_matches.size(); i++ )
  {
    keypoints1_best.push_back(keypoints1[good_matches[i].queryIdx]);
    keypoints2_best.push_back(keypoints2[good_matches[i].trainIdx]);
  }
  std::vector<cv::KeyPoint> keypoints1_best2, keypoints2_best2;
  int j=0;
  for (int i = 0; i < mask.rows; i++) {
    if (mask.at<uchar>(i)) {
      keypoints1_best2.push_back(keypoints1_best[i]);
      keypoints2_best2.push_back(keypoints2_best[i]);
      top_matches.push_back(cv::DMatch(static_cast<int>(j), static_cast<int>(j), 0));
      j++;
    }
  }
  scaledDrawMatches(im1Proc, keypoints1_best2, im2Proc, keypoints2_best2, top_matches, imMatches);
  log_mem_usage("draw matches");
  
  // TODO: can I release there now ? 
  // std::vector<cv::KeyPoint>().swap(keypoints1_best);
  // std::vector<cv::KeyPoint>().swap(keypoints2_best);
  // std::vector<cv::KeyPoint>().swap(keypoints1_best2);
  // std::vector<cv::KeyPoint>().swap(keypoints2_best2);
  // std::vector<cv::DMatch>().swap(top_matches);
  
  // check number of matches
  return check_matches(mask);
}

void getSIFTTransformationMatrix(
    Mat &im1Proc, Mat &im2Proc, Mat &h, Mat &mask, 
    Mat &imMatches, std::vector<Point2f> &points1, std::vector<Point2f> &points2, 
    const bool &run_Affine, bool &is_faulty){
  
  // parameters
  SIFTParameters params;
  
  //////////////////////////
  /// Run multiple SIFT ////
  //////////////////////////
  
  // images
  Mat im1Proc_eq, im1Proc_eq2;
  Mat im2Proc_eq, im2Proc_eq2;
  
  // check variable
  bool check;
  Rcout << "MESSAGE: Calculating" << (run_Affine ? " (Affine) " : " (Homography) ") << "Transformation Matrix" << endl;
  
  // find matches and points
  check = getSIFTTransformationMatrixSingle(im1Proc, im2Proc, h, mask,
                                            imMatches,
                                            points1, points2,
                                            run_Affine, params, is_faulty);
  Rcout << "DONE: calculated homography matrix with " << points1.size() << " points" << endl;
  
  // equalize first image if fails
  if(!check){

    Mat im1Proc_eq;
    cv::equalizeHist(im1Proc, im1Proc_eq);
    Rcout << "MESSAGE: Calculating Transformation Matrix with histogram equalization (1)" << endl;

    check = getSIFTTransformationMatrixSingle(im1Proc_eq, im2Proc, h, mask,
                                              imMatches,
                                              points1, points2,
                                              run_Affine, params, is_faulty);
    Rcout << "DONE: calculated homography matrix with " << points1.size() << " points" << endl;
  } else {
    return;
  }

  // equalize second image if fails
  if(!check){

    cv::equalizeHist(im2Proc, im2Proc_eq);
    Rcout << "MESSAGE: Calculating Transformation Matrix with histogram equalization (2)" << endl;

    check = getSIFTTransformationMatrixSingle(im1Proc, im2Proc_eq, h, mask,
                                              imMatches,
                                              points1, points2,
                                              run_Affine, params, is_faulty);
    Rcout << "DONE: calculated homography matrix with " << points1.size() << " points" << endl;
  } else {
    return;
  }

  // last try with both equalized images
  if(!check){

    cv::equalizeHist(im1Proc, im1Proc_eq2);
    cv::equalizeHist(im2Proc, im2Proc_eq2);
    Rcout << "MESSAGE: Calculating Transformation Matrix with histogram equalization (3)" << endl;

    check = getSIFTTransformationMatrixSingle(im1Proc_eq2, im2Proc_eq2, h, mask,
                                              imMatches,
                                              points1, points2,
                                              run_Affine, params, is_faulty);
    Rcout << "DONE: calculated homography matrix with " << points1.size() << " points" << endl;
  } else {
    return;
  }
  
  // TODO: release ?
  // im1Proc_eq.release();
  // im1Proc_eq2.release();
  // im2Proc_eq.release();
  // im2Proc_eq2.release();
}

bool getORBTransformationMatrix(
    Mat &im1Proc, Mat &im2Proc, Mat &h, Mat &mask, 
    Mat &imMatches, std::vector<Point2f> &points1, std::vector<Point2f> &points2, 
    const bool &run_Affine, const float GOOD_MATCH_PERCENT, const int MAX_FEATURES, bool &is_faulty){
  
  //////////////////////
  /// Compute ORB /////
  //////////////////////
  
  // Variables to store keypoints and descriptors
  std::vector<KeyPoint> keypoints1, keypoints2;
  Mat descriptors1, descriptors2;
  
  // Detect ORB features and compute descriptors.
  Ptr<Feature2D> orb = ORB::create(MAX_FEATURES);
  orb->detectAndCompute(im1Proc, Mat(), keypoints1, descriptors1);
  orb->detectAndCompute(im2Proc, Mat(), keypoints2, descriptors2);
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
  // std::vector<Point2f> points1, points2;
  for( size_t i = 0; i < matches.size(); i++ )
  {
    points1.push_back(keypoints1[ matches[i].queryIdx ].pt );
    points2.push_back(keypoints2[ matches[i].trainIdx ].pt );
  }
  
  // check variable
  Rcout << "MESSAGE: Calculating" << (run_Affine ? " (Affine) " : " (Homography) ") << "Transformation Matrix" << endl;
    
  // Find transformation matrix
  if(points1.size() > 0){
    if(run_Affine){
      std::vector<uint8_t> match_mask;
      h = estimateAffine2D(points1,
                           points2,
                           match_mask,
                           cv::RANSAC);
      mask = IntVectorToMat(match_mask);
    } else {
      h = findHomography(points1,
                         points2,
                         cv::RANSAC,
                         5,
                         mask);
    }
  } else {
    Rcout <<  "Found no matches!" << endl;
    return false;
  }
  
  // Draw top matches and good ones only
  std::vector<cv::DMatch> top_matches;
  std::vector<cv::KeyPoint> keypoints1_best, keypoints2_best;
  for(size_t i = 0; i < matches.size(); i++ )
  {
    keypoints1_best.push_back(keypoints1[matches[i].queryIdx]);
    keypoints2_best.push_back(keypoints2[matches[i].trainIdx]);
  }
  std::vector<cv::KeyPoint> keypoints1_best2, keypoints2_best2;
  int j=0;
  for (int i = 0; i < mask.rows; i++) {
    if (mask.at<uchar>(i)) {
      keypoints1_best2.push_back(keypoints1_best[i]);
      keypoints2_best2.push_back(keypoints2_best[i]);
      top_matches.push_back(cv::DMatch(static_cast<int>(j), static_cast<int>(j), 0));
      j++;
    }
  }
  scaledDrawMatches(im1Proc, keypoints1_best2, im2Proc, keypoints2_best2, top_matches, imMatches);
  
  // check number of matches
  return check_matches(mask);
}

// align images with FLANN algorithm
void alignImages(Mat &im1, Mat &im2, Mat &im1Reg, Mat &im1Overlay, 
                 Mat &imMatches, Mat &h, Rcpp::List &keypoints,
                 const float GOOD_MATCH_PERCENT, const int MAX_FEATURES,
                 Rcpp::String matcher,
                 const bool invert_query, const bool invert_ref,
                 const char* flipflop_query, const char* flipflop_ref,
                 const char* rotate_query, const char* rotate_ref,
                 const bool run_Affine, const bool run_TPS)
{
  
  // parameters
  cv::setRNGSeed(0);

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

  // RUN SIFT+FLANN+Homography with retry
  bool is_faulty = FALSE;
  cv::Mat mask;
  std::vector<Point2f> points1, points2;
  if(strcmp(matcher.get_cstring(), "BRUTE-FORCE") == 0){
    
    // message
    Rcout << "MESSAGE: Running BRUTE-FORCE Alignment" << endl;
    
    // run ORB
    bool check;
    check = getORBTransformationMatrix(im1Proc, im2Proc, h, mask, imMatches,
                                       points1, points2, run_Affine, 
                                       GOOD_MATCH_PERCENT, MAX_FEATURES, is_faulty);
    
  } else {
    
    // message
    Rcout << "MESSAGE: Running SIFT+FLANN Alignment" << ((run_TPS) ? " with TPS" : "") << endl;
    
    // run SIFT
    getSIFTTransformationMatrix(im1Proc, im2Proc, h, mask, imMatches,
                                points1, points2, run_Affine, is_faulty);

  }
  
  // check result
  is_faulty = check_transformation_metrics(points1, points2, im2, h, mask);
  Rcout << "MESSAGE: Registration is " << (is_faulty ? "degenerate!" : "not degenerate!") << endl;
  
  // Use homography to warp image
  Mat im1Warp, im1NormalWarp;
  if(h.rows == 2){
    warpAffine(im1Proc, im1Warp, h, im2Proc.size());
    warpAffine(im1NormalProc, im1NormalWarp, h, im2Proc.size());   
  } else if(h.rows == 3){
    warpPerspective(im1Proc, im1Warp, h, im2Proc.size());
    warpPerspective(im1NormalProc, im1NormalWarp, h, im2Proc.size());    
  } else {
    Rcout << "WARNING: No transformation was found" << endl;
    return;
  }
  
  Rcout << "DONE: warped query image" << endl;
  log_mem_usage("warp image");
  
  ///////////////////////
  /// Find Homography ///
  ///////////////////////
  
  // continue with TPS or do FLANN only
  Mat im1Reg_Warp_nonrigid;
  Mat im1Reg_NormalWarp_nonrigid;
  Mat im1Combine;
  if(is_faulty || !run_TPS){
    
    // change color map
    cv::addWeighted(im2Proc, 0.7, im1Warp, 0.3, 0, im1Combine);

    // Reverse process
    im1Reg = reversepreprocessImage(im1NormalWarp, flipflop_ref, rotate_ref);

    // return as rgb
    cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);
    cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);

    // TPS is requested (only if FLANN succeeded)
  } else {
    
    Rcout << "MESSAGE: Running Thin-Plate-Spline Alignment" << endl;
    
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
    if (h.rows == 2){
      cv::transform(filtered_points1, filtered_points1_reg, h);
    } else {
      cv::perspectiveTransform(filtered_points1, filtered_points1_reg, h);
    }
    
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
    im1Reg_NormalWarp_nonrigid_cropped.release();
    
    cv::Mat im1Reg_Warp_nonrigid_cropped  = im1Reg_Warp_nonrigid(cv::Range(0,im2Proc.size().height), cv::Range(0,im2Proc.size().width));
    im1Reg_Warp_nonrigid = im1Reg_Warp_nonrigid_cropped.clone();
    im1Reg_Warp_nonrigid_cropped.release();
    
    // change color map
    cv::addWeighted(im2Proc, 0.7, im1Reg_Warp_nonrigid, 0.3, 0, im1Combine);
    
    // Reverse process
    im1Reg = reversepreprocessImage(im1Reg_NormalWarp_nonrigid, flipflop_ref, rotate_ref);
    im1Reg_NormalWarp_nonrigid.release();
    
    // return as rgb
    cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);
    cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);
    
  }
  
  // release
  im1Combine.release();
  im2Proc.release();
  im1Warp.release();
  im1NormalWarp.release();
  im1NormalProc.release();
  
  // resize image to visualize faster later in Shiny
  im2 = resize_image(im2, 500);
  im1Overlay = resize_image(im1Overlay, 500);
}

// [[Rcpp::export]]
Rcpp::List automated_registeration_rawvector(Rcpp::RawVector& ref_image, Rcpp::RawVector& query_image,
                                             const int width1, const int height1,
                                             const int width2, const int height2,
                                             const float GOOD_MATCH_PERCENT, const int MAX_FEATURES,
                                             const bool invert_query, const bool invert_ref,
                                             Rcpp::String flipflop_query, Rcpp::String flipflop_ref,
                                             Rcpp::String rotate_query, Rcpp::String rotate_ref,
                                             Rcpp::String matcher, Rcpp::String method)
{
  log_mem_usage("pre conversion");
  
  // Return data
  Rcpp::List out(5);
  Rcpp::List out_trans(2);
  Rcpp::List keypoints(2);
  Mat imOverlay, imReg, h, imMatches;

  // Read reference image
  cv::Mat imReference = imageToMat(ref_image, width1, height1);
  Rcout << ref_image.size()  << endl;
  Rcout << ref_image.length() << endl;
  Rcout << object_size(ref_image.size() * sizeof(ref_image[0])) << endl;
  Rcout << object_size(imReference.total() * imReference.elemSize()) << endl;
  Rcout << object_size(imReference.step[0] * imReference.rows) << endl;
  
  // Read image to be aligned
  cv::Mat im = imageToMat(query_image, width2, height2);
  Rcout << query_image.size() << endl;
  Rcout << query_image.length() << endl;
  Rcout << object_size(query_image.size()) << endl;
  Rcout << object_size(im.total() * im.elemSize()) << endl;
  Rcout << object_size(im.step[0] * im.rows) << endl;
  
  // run alignment
  const bool run_TPS = (strcmp(method.get_cstring(), "Homography + Non-Rigid") == 0 || 
                        strcmp(method.get_cstring(), "Affine + Non-Rigid") == 0);
  const bool run_Affine = (strcmp(method.get_cstring(), "Affine") == 0 || 
                           strcmp(method.get_cstring(), "Affine + Non-Rigid") == 0);
  log_mem_usage("pre alignment");
  alignImages(im, imReference, imReg, imOverlay, imMatches,
              h, keypoints,
              GOOD_MATCH_PERCENT, MAX_FEATURES,
              matcher,
              invert_query, invert_ref,
              flipflop_query.get_cstring(), flipflop_ref.get_cstring(),
              rotate_query.get_cstring(), rotate_ref.get_cstring(),
              run_Affine, run_TPS);

  // transformation matrix, can be either a matrix, set of keypoints or both
  out_trans[0] = matToNumericMatrix(h.clone());
  out_trans[1] = keypoints;
  out[0] = out_trans;
  
  // destination image, registered image, keypoint matching image
  out[1] = matToImage(imReference.clone());

  // check if transformation matrix is calculated, 
  // otherwise return NULL
  if(h.rows > 1){
    // registered image
    out[2] = matToImage(imReg.clone());
    // keypoint matching image
    out[3] = matToImage(imMatches.clone());
    // overlay image
    out[4] = matToImage(imOverlay.clone());
  } else {
    out[2] = R_NilValue;
    out[3] = R_NilValue;
    out[4] = R_NilValue;
  }
  

  // release
  im.release();
  imReference.release();
  imReg.release();
  imMatches.release();
  imOverlay.release();
  
  log_mem_usage("now");
  
  // return
  return out;
}

/////////////////
/// scratch /////
/////////////////

// align images with BRUTE FORCE algorithm
void alignImagesBRUTE(Mat &im1, Mat &im2, Mat &im1Reg, Mat &im1Overlay, Mat &imMatches, Mat &h,
                      const float GOOD_MATCH_PERCENT, const int MAX_FEATURES,
                      const bool invert_query, const bool invert_ref,
                      const char* flipflop_query, const char* flipflop_ref,
                      const char* rotate_query, const char* rotate_ref,
                      const bool run_Affine)
{
  
  // Convert images to grayscale
  Mat im1Gray, im2Gray;
  cvtColor(im1, im1Gray, cv::COLOR_BGR2GRAY);
  cvtColor(im2, im2Gray, cv::COLOR_BGR2GRAY);
  
  // Variables to store keypoints and descriptors
  std::vector<KeyPoint> keypoints1, keypoints2;
  Mat descriptors1, descriptors2;
  
  // Process images
  Mat im1Proc, im2Proc, im1NormalProc;
  im1Proc = preprocessImage(im1Gray, invert_query, flipflop_query, rotate_query);
  im1NormalProc = preprocessImage(im1, FALSE, flipflop_query, rotate_query);
  im2Proc = preprocessImage(im2Gray, invert_ref, flipflop_ref, rotate_ref);
  
  // Detect ORB features and compute descriptors.
  Ptr<Feature2D> orb = ORB::create(MAX_FEATURES);
  orb->detectAndCompute(im1Proc, Mat(), keypoints1, descriptors1);
  orb->detectAndCompute(im2Proc, Mat(), keypoints2, descriptors2);
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
  
  // check variable
  Rcout << "Calculating" << (run_Affine ? " (Affine) " : " (Homography) ") << "Transformation Matrix" << endl;
  
  // Find transformation matrix
  cv::Mat mask;
  if(run_Affine){
    std::vector<uint8_t> match_mask;
    h = estimateAffine2D(points1,
                         points2,
                         match_mask,
                         cv::RANSAC);
    mask = IntVectorToMat(match_mask);
  } else {
    h = findHomography(points1,
                       points2,
                       cv::RANSAC,
                       5,
                       mask);
  }
  
  // Draw top matches and good ones only
  std::vector<cv::DMatch> top_matches;
  std::vector<cv::KeyPoint> keypoints1_best, keypoints2_best;
  for(size_t i = 0; i < matches.size(); i++ )
  {
    keypoints1_best.push_back(keypoints1[matches[i].queryIdx]);
    keypoints2_best.push_back(keypoints2[matches[i].trainIdx]);
  }
  std::vector<cv::KeyPoint> keypoints1_best2, keypoints2_best2;
  int j=0;
  for (int i = 0; i < mask.rows; i++) {
    if (mask.at<uchar>(i)) {
      keypoints1_best2.push_back(keypoints1_best[i]);
      keypoints2_best2.push_back(keypoints2_best[i]);
      top_matches.push_back(cv::DMatch(static_cast<int>(j), static_cast<int>(j), 0));
      j++;
    }
  }
  scaledDrawMatches(im1Proc, keypoints1_best2, im2Proc, keypoints2_best2, top_matches, imMatches);
  
  // Use homography to warp image
  Mat im1Warp, im1NormalWarp;
  if(h.rows == 2){
    warpAffine(im1Proc, im1Warp, h, im2Proc.size());
    warpAffine(im1NormalProc, im1NormalWarp, h, im2Proc.size());   
  } else {
    warpPerspective(im1Proc, im1Warp, h, im2Proc.size());
    warpPerspective(im1NormalProc, im1NormalWarp, h, im2Proc.size());    
  }
  
  // Reverse process
  im1Reg = reversepreprocessImage(im1NormalWarp, flipflop_ref, rotate_ref);
  
  // return as rgb
  cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);
  
  // resize image to visualize faster later in Shiny
  im2 = resize_image(im2, 500);
  im1Overlay = resize_image(im1Reg, 500);
}

// align images with FLANN algorithm
void alignImagesFLANN(Mat &im1, Mat &im2, Mat &im1Reg, Mat &im1Overlay, 
                      Mat &imMatches, Mat &h, Rcpp::List &keypoints,
                      const bool invert_query, const bool invert_ref,
                      const char* flipflop_query, const char* flipflop_ref,
                      const char* rotate_query, const char* rotate_ref,
                      const bool run_Affine, const bool run_TPS)
{
  
  // parameters
  cv::setRNGSeed(0);
  SIFTParameters params;
  
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
  
  // RUN SIFT+FLANN+Homography with retry
  bool is_faulty = FALSE;
  cv::Mat mask;
  std::vector<Point2f> points1, points2;
  // getSIFTTransformationMatrix(im1Proc, im2Proc, im1, im2, h, mask, imMatches,
  //                             points1, points2, run_Affine, params, is_faulty);
  
  // check result
  is_faulty = check_transformation_metrics(points1, points2, im2, h, mask);
  Rcout << "MESSAGE: Registration is " << (is_faulty ? "degenerate!" : "not degenerate!") << endl;
  
  // Use homography to warp image
  Mat im1Warp, im1NormalWarp;
  if(h.rows == 2){
    warpAffine(im1Proc, im1Warp, h, im2Proc.size());
    warpAffine(im1NormalProc, im1NormalWarp, h, im2Proc.size());   
  } else {
    warpPerspective(im1Proc, im1Warp, h, im2Proc.size());
    warpPerspective(im1NormalProc, im1NormalWarp, h, im2Proc.size());    
  }
  
  Rcout << "DONE: warped query image" << endl;
  
  ///////////////////////
  /// Find Homography ///
  ///////////////////////
  
  // continue with TPS or do FLANN only
  Mat im1Reg_Warp_nonrigid;
  Mat im1Reg_NormalWarp_nonrigid;
  Mat im1Combine;
  if(is_faulty || !run_TPS){
    
    // change color map
    cv::addWeighted(im2Proc, 0.7, im1Warp, 0.3, 0, im1Combine);
    
    // Reverse process
    im1Reg = reversepreprocessImage(im1NormalWarp, flipflop_ref, rotate_ref);
    
    // return as rgb
    cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);
    cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);
    
    // TPS is requested (only if FLANN succeeded)
  } else {
    
    Rcout << "MESSAGE: Running Thin-Plate-Spline Alignment" << endl;
    
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
    if (h.rows == 2){
      cv::transform(filtered_points1, filtered_points1_reg, h);
    } else {
      cv::perspectiveTransform(filtered_points1, filtered_points1_reg, h);
    }
    
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
    cv::addWeighted(im2Proc, 0.7, im1Reg_Warp_nonrigid, 0.3, 0, im1Combine);
    
    // Reverse process
    im1Reg = reversepreprocessImage(im1Reg_NormalWarp_nonrigid, flipflop_ref, rotate_ref);
    
    // return as rgb
    cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);
    cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);
  }
  
  // resize image to visualize faster later in Shiny
  im2 = resize_image(im2, 500);
  im1Overlay = resize_image(im1Overlay, 500);
}

// align images with FLANN algorithm
void alignImagesFLANN2(Mat &im1, Mat &im2, Mat &im1Reg, Mat &im1Overlay, 
                       Mat &imMatches, Mat &h, Rcpp::List &keypoints,
                       const bool invert_query, const bool invert_ref,
                       const char* flipflop_query, const char* flipflop_ref,
                       const char* rotate_query, const char* rotate_ref,
                       const bool run_Affine, const bool run_TPS)
{
  
  // parameters
  cv::setRNGSeed(0);
  SIFTParameters params;
  
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
  
  // RUN SIFT+FLANN+Homography with retry
  bool is_faulty = FALSE;
  cv::Mat mask;
  std::vector<Point2f> points1, points2;
  
  
  // Variables to store keypoints and descriptors
  std::vector<KeyPoint> keypoints1, keypoints2;
  Mat descriptors1, descriptors2;
  
  // Detect SIFT features
  // Ptr<Feature2D> sift = cv::SIFT::create(params.sift_nfeatures);
  Ptr<Feature2D> sift = cv::SIFT::create();
  // computeSIFTTiles(im1Proc, keypoints1, descriptors1, sift, params);
  // computeSIFTTiles(im2Proc, keypoints2, descriptors2, sift, params);
  sift->detectAndCompute(im1Proc, Mat(), keypoints1, descriptors1);
  sift->detectAndCompute(im2Proc, Mat(), keypoints2, descriptors2);
  
  Rcout << "MESSAGE: Generated " << keypoints1.size() << " and " << keypoints2.size() << " keypoints"  << endl;
  Rcout << "DONE: SIFT based key-points detection and descriptors computation" << endl;
  
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
  // goodMatches = get_good_matches(matches)
  std::vector<DMatch> good_matches;
  getGoodMatches_temp(matches, good_matches);
  Rcout << "DONE: get good matches by distance thresholding" << endl;
  
  ///////////////////////
  /// Find Homography ///
  ///////////////////////
  
  // Extract location of good matches
  for( size_t i = 0; i < good_matches.size(); i++ )
  {
    points1.push_back(keypoints1[good_matches[i].queryIdx].pt);
    points2.push_back(keypoints2[good_matches[i].trainIdx].pt);
  }
  
  // check variable
  Rcout << "MESSAGE: Calculating" << (run_Affine ? " (Affine) " : " (Homography) ") << "Transformation Matrix" << endl;
  
  // Find transformation matrix
  Rcout << "MESSAGE: Matching " << points1.size() << " keypoints" << endl;
  if(run_Affine){
    std::vector<uint8_t> match_mask;
    // h = estimateAffine2D(points1,
    //                      points2,
    //                      match_mask,
    //                      cv::RANSAC,
    //                      params.ransac_pixel_threshold,
    //                      params.ransac_maxIters,
    //                      params.ransac_confidence);
    h = estimateAffine2D(points1,
                         points2,
                         match_mask,
                         cv::RANSAC);
    mask = IntVectorToMat(match_mask);
  } else {
    // h = findHomography(points1,
    //                    points2,
    //                    cv::RANSAC,
    //                    params.ransac_pixel_threshold,
    //                    mask,
    //                    params.ransac_maxIters,
    //                    params.ransac_confidence);
    h = findHomography(points1, points2, RANSAC);
  } 
  
  // Draw top matches and good ones only
  std::vector<cv::DMatch> top_matches;
  std::vector<cv::KeyPoint> keypoints1_best, keypoints2_best;
  for(size_t i = 0; i < good_matches.size(); i++ )
  {
    keypoints1_best.push_back(keypoints1[good_matches[i].queryIdx]);
    keypoints2_best.push_back(keypoints2[good_matches[i].trainIdx]);
  }
  std::vector<cv::KeyPoint> keypoints1_best2, keypoints2_best2;
  int j=0;
  for (int i = 0; i < mask.rows; i++) {
    if (mask.at<uchar>(i)) {
      keypoints1_best2.push_back(keypoints1_best[i]);
      keypoints2_best2.push_back(keypoints2_best[i]);
      top_matches.push_back(cv::DMatch(static_cast<int>(j), static_cast<int>(j), 0));
      j++;
    }
  }
  // scaledDrawMatches(im1Proc, keypoints1_best2, im2Proc, keypoints2_best2, top_matches, imMatches);
  drawMatches(im1Proc, keypoints1_best2, im2Proc, keypoints2_best2, top_matches, imMatches); 
  
  // check result
  is_faulty = check_transformation_metrics(points1, points2, im2, h, mask);
  Rcout << "MESSAGE: Registration is " << (is_faulty ? "degenerate!" : "not degenerate!") << endl;
  
  // Use homography to warp image
  Mat im1Warp, im1NormalWarp;
  if(h.rows == 2){
    warpAffine(im1Proc, im1Warp, h, im2Proc.size());
    warpAffine(im1NormalProc, im1NormalWarp, h, im2Proc.size());   
  } else {
    warpPerspective(im1Proc, im1Warp, h, im2Proc.size());
    warpPerspective(im1NormalProc, im1NormalWarp, h, im2Proc.size());    
  }
  
  Rcout << "DONE: warped query image" << endl;
  
  ///////////////////////
  /// Find Homography ///
  ///////////////////////
  
  // continue with TPS or do FLANN only
  Mat im1Reg_Warp_nonrigid;
  Mat im1Reg_NormalWarp_nonrigid;
  Mat im1Combine;
  if(is_faulty || !run_TPS){
    
    // change color map
    cv::addWeighted(im2Proc, 0.7, im1Warp, 0.3, 0, im1Combine);
    
    // Reverse process
    im1Reg = reversepreprocessImage(im1NormalWarp, flipflop_ref, rotate_ref);
    
    // return as rgb
    cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);
    cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);
    
    // TPS is requested (only if FLANN succeeded)
  } else {
    
    Rcout << "MESSAGE: Running Thin-Plate-Spline Alignment" << endl;
    
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
    if (h.rows == 2){
      cv::transform(filtered_points1, filtered_points1_reg, h);
    } else {
      cv::perspectiveTransform(filtered_points1, filtered_points1_reg, h);
    }
    
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
    cv::addWeighted(im2Proc, 0.7, im1Reg_Warp_nonrigid, 0.3, 0, im1Combine);
    
    // Reverse process
    im1Reg = reversepreprocessImage(im1Reg_NormalWarp_nonrigid, flipflop_ref, rotate_ref);
    
    // return as rgb
    cvtColor(im1Combine, im1Overlay, cv::COLOR_GRAY2BGR);
    cvtColor(im2Proc, im2, cv::COLOR_GRAY2BGR);
  }
  
  // resize image to visualize faster later in Shiny
  im2 = resize_image(im2, 500);
  im1Overlay = resize_image(im1Overlay, 500);
}