#include "Rcpp.h"
#include <opencv2/opencv.hpp>
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/shape/shape_transformer.hpp"

using namespace Rcpp;
using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

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
Rcpp::RawVector matToImage(cv::Mat mat) {

  // Create RawVector object
  Rcpp::RawVector rawvec(mat.total() * mat.elemSize());
  rawvec.attr("dim") = Rcpp::Dimension(3, mat.cols, mat.rows);

  // Copy Mat data to RawVector
  std::memcpy(rawvec.begin(), mat.data, rawvec.size());

  return rawvec;
}

// Function to convert a RawVector for magick images to a cv::Mat object
cv::Mat imageToMat(Rcpp::RawVector image_data, int width, int height) {
  
  // Create cv::Mat object
  cv::Mat mat(height, width, CV_8UC3, image_data.begin());
  
  // Convert from RGBA to BGRA
  cv::cvtColor(mat, mat, cv::COLOR_RGBA2BGR);
  
  return mat;
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

// Function to convert a cv::Point2f object to a cv::Mat
// std::vector<uint8_t> matToIntVector(cv::Mat mat) {
//   std::vector<uint8_t> points;
//   
//   // Iterate over the rows of the matrix
//   for (int i = 0; i < mat.rows; ++i) {
//     
//     // Create Point2f object and add it to the vector
//     points.push_back(mat.at<int>(i,0));
//   }
//   
//   return points;
// }

cv::Mat IntVectorToMat(std::vector<uint8_t> points) {
  cv::Mat mat(points.size(), 1, CV_8U);
  // cv::Mat mat;

  // Iterate over the vector of Point2f
  for (size_t i = 0; i < points.size(); ++i) {
    mat.at<int>(i, 0) = (int) points[i];
  }

  return mat;
}

// cv::Mat IntVectortoMat(std::vector<uint8_t> vec) {
//   return cv::Mat(vec).reshape(1, vec.size()); // Create a single-column matrix
// }

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