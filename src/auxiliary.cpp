#include <Rcpp.h>
// #include <sys/resource.h>
// #include <sys/sysctl.h>
// #include <unistd.h>
// #include <mach/vm_statistics.h>
// #include <mach/mach.h>

// OpenCV
#include <opencv2/opencv.hpp>
#include "opencv2/shape/shape_transformer.hpp"

// Internal functions
#include "auxiliary.h"

using namespace Rcpp;
using namespace std;
using namespace cv;

////
// replace NAs in matrices
////

// [[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix replaceNaMatrix(Rcpp::NumericMatrix mat, int replace) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if (Rcpp::NumericMatrix::is_na(mat(i, j))) {
        mat(i, j) = replace;
      }
    }
  }
  return mat;
}

////
// Conversion
////

// Function to convert a cv::Mat object to a RawVector for magick images
Rcpp::RawVector matToImage(const cv::Mat &mat) {
  
  // Create RawVector object
  Rcpp::RawVector rawvec(mat.total() * mat.elemSize());
  rawvec.attr("dim") = Rcpp::Dimension(3, mat.cols, mat.rows);

  // Copy Mat data to RawVector
  std::memcpy(rawvec.begin(), mat.data, rawvec.size());

  return rawvec;
}

// Function to convert a RawVector for magick images to a cv::Mat object
cv::Mat imageToMat(Rcpp::RawVector &image_data, int width, int height) {
  
  // Create cv::Mat object
  cv::Mat mat(height, width, CV_8UC3, image_data.begin());
  
  // Convert from RGBA to BGRA
  cv::cvtColor(mat, mat, cv::COLOR_RGBA2BGR);
  
  return mat;
}

// Function to convert a cv::Mat object to a RawVector for magick images
Rcpp::IntegerVector matToMask(const cv::Mat &mat) {

  cv::Mat intMat;
  mat.convertTo(intMat, CV_32S);
  Rcpp::IntegerVector intvec(intMat.total());
  std::memcpy(
    intvec.begin(),
    intMat.data,
    static_cast<std::size_t>(intvec.size()) * sizeof(int)
  );
  intvec.attr("dim") = Rcpp::Dimension(intMat.cols, intMat.rows);
  
  return intvec;
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
Rcpp::NumericMatrix point2fToNumericMatrix(std::vector<cv::Point2f> &points) {
  int n = points.size();
  Rcpp::NumericMatrix mat(n, 2);
  for (int i = 0; i < n; i++) {
    mat(i, 0) = points[i].x;
    mat(i, 1) = points[i].y;
  }
  return mat;
}

// Function to convert a cv::Keypoint object to a std::vector<double>
std::vector<double> KeyPointToDoubleVector(std::vector<cv::KeyPoint> &points,
                                           int axes = 0) {
  int n = points.size();
  std::vector<double> vec(n);
  for (int i = 0; i < n; i++) {
    if(axes == 0){
      vec[i] = (double) points[i].pt.x; 
    } else {
      vec[i] = (double) points[i].pt.y; 
    }
  }
  return vec;
}

// Function to convert a cv::Point2f object to a std::vector<double>
std::vector<double> Point2fToDoubleVector(std::vector<cv::Point2f> &points, 
                                          int axes = 0) {
  int n = points.size();
  std::vector<double> vec(n);
  for (int i = 0; i < n; i++) {
    if(axes == 0){
      vec[i] = (double) points[i].x; 
    } else {
      vec[i] = (double) points[i].y; 
    }
  }
  return vec;
}

// Function to convert a cv::Keypoint object to a std::vector<cv:Point2f>
std::vector<cv::Point2f> KeyPointToPoint2f(std::vector<cv::KeyPoint> &keypoints) {
  int n = keypoints.size();
  std::vector<cv::Point2f> points;
  
  for (int i = 0; i < n; i++) {
    points.push_back(keypoints[i].pt);
  }
  return points;
}

// Function to convert a cv::Point2f object to a cv::Mat
std::vector<cv::Point2f> matToPoint2f(cv::Mat &mat) {
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

cv::Mat point2fToMat(std::vector<cv::Point2f> &points) {
  cv::Mat mat(points.size(), 2, CV_32F);
  
  // Iterate over the vector of Point2f
  for (size_t i = 0; i < points.size(); ++i) {
    
    // Assign x and y coordinates to the matrix
    mat.at<float>(i, 0) = points[i].x;
    mat.at<float>(i, 1) = points[i].y;
  }
  
  return mat;
}

// Function to convert a cv::Point2f object to a cv::Mat
cv::Mat IntVectorToMat(std::vector<uint8_t> &points) {
  cv::Mat mat(points.size(), 1, CV_8U);

  // Iterate over the vector of Point2f
  for (size_t i = 0; i < points.size(); ++i) {
    // mat.at<int>(i, 0) = (int) points[i]; // causing heap error
    mat.at<uint8_t>(i, 0) = points[i];
  }

  return mat;
}

////
// stats
////

// calculate standard deviation of a vector
double cppSD(std::vector<cv::KeyPoint> &points, int axes)
{
  std::vector<double> inVec = KeyPointToDoubleVector(points, axes);
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
  std::vector<double>().swap(inVec);
  return std::sqrt( sd / (n-1) );
}

double cppSD(std::vector<cv::Point2f> &points, int axes)
{
  std::vector<double> inVec = Point2fToDoubleVector(points, axes);
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
  std::vector<double>().swap(inVec);
  return std::sqrt( sd / (n-1) );
}

double meanDistances(std::vector<cv::Point2f>& pts1,
                    std::vector<cv::Point2f>& pts2)
{
  if (pts1.size() != pts2.size() || pts1.empty())
    return 0.0;

  double sumDist = 0.0;
  for (size_t i = 0; i < pts1.size(); ++i)
  {
    const double dx = pts1[i].x - pts2[i].x;
    const double dy = pts1[i].y - pts2[i].y;
    sumDist += std::sqrt(dx * dx + dy * dy);
  }
  
  return sumDist / pts1.size();
}

double medianDistances(std::vector<cv::Point2f>& pts1,
                      std::vector<cv::Point2f>& pts2)
{
  if (pts1.size() != pts2.size() || pts1.empty())
    return 0.0;
  
  std::vector<double> distances;
  distances.reserve(pts1.size());
  
  for (size_t i = 0; i < pts1.size(); ++i)
  {
    const double dx = pts1[i].x - pts2[i].x;
    const double dy = pts1[i].y - pts2[i].y;
    distances.push_back(std::sqrt(dx * dx + dy * dy));
  }
  
  const size_t n = distances.size();
  const size_t mid = n / 2;
  
  std::nth_element(distances.begin(),
                   distances.begin() + mid,
                   distances.end());
  
  if (n % 2 == 1)
  {
    return distances[mid];
  }
  else
  {
    double upper = distances[mid];
    
    std::nth_element(distances.begin(),
                     distances.begin() + mid - 1,
                     distances.end());
    
    double lower = distances[mid - 1];
    
    return (lower + upper) / 2.0;
  }
}