#include <Rcpp.h>
#include <sys/resource.h>
#include <sys/sysctl.h>
#include <unistd.h>
#include <mach/vm_statistics.h>
#include <mach/mach.h>

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
// memory
////

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

void log_mem_macos(const std::string& label = "") {
  mach_task_basic_info info;
  mach_msg_type_number_t size = MACH_TASK_BASIC_INFO_COUNT;
  kern_return_t kr = task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                               (task_info_t)&info, &size);
  
  if (kr != KERN_SUCCESS) {
    Rcpp::Rcerr << "[MEM " << label << "] Failed to get memory info.\n";
    return;
  }
  
  double rss_gb = static_cast<double>(info.resident_size) / (1024.0 * 1024.0 * 1024.0);
  double virt_gb = static_cast<double>(info.virtual_size) / (1024.0 * 1024.0 * 1024.0);
  
  Rcpp::Rcout << "[MEM " << label << "] Resident (RSS): "
              << rss_gb << " GB, Virtual: " << virt_gb << " GB\n";
}

double object_size_long(long bsize) {
  
  double rss_kb = bsize / 1024.0;
  double rss_mb = rss_kb / 1024.0;
  double rss_gb = rss_mb / 1024.0;
  
  return rss_gb;
}

double object_size_double(double bsize) {
  
  double rss_kb = bsize / 1024;
  double rss_mb = rss_kb / 1024;
  double rss_gb = rss_mb / 1024;
  
  return rss_gb;
}

double get_resident_bytes() {
  mach_task_basic_info info;
  mach_msg_type_number_t size = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                (task_info_t)&info, &size) != KERN_SUCCESS) {
    return 0;
  }
  return static_cast<double>(info.resident_size);
}

double bytes_to_gb(double bytes) {
  return bytes / (1024.0 * 1024.0 * 1024.0);
}

////
// Conversion
////

// Function to convert a cv::Mat object to a RawVector for magick images
Rcpp::RawVector matToImage(const cv::Mat &mat) {
  // profiler
  // MemProfiler mp("Mat -> Image");
  
  // Create RawVector object
  Rcpp::RawVector rawvec(mat.total() * mat.elemSize());
  rawvec.attr("dim") = Rcpp::Dimension(3, mat.cols, mat.rows);

  // Copy Mat data to RawVector
  std::memcpy(rawvec.begin(), mat.data, rawvec.size());

  return rawvec;
}

// Function to convert a RawVector for magick images to a cv::Mat object
cv::Mat imageToMat(Rcpp::RawVector &image_data, int width, int height) {
  // profiler
  // MemProfiler mp("Image -> Mat");
  
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
std::vector<double> KeyPointToDoubleVector(std::vector<cv::KeyPoint> &points) {
  int n = points.size();
  std::vector<double> vec(n);
  for (int i = 0; i < n; i++) {
    vec[i] = (double) points[i].pt.x;
  }
  return vec;
}

// Function to convert a cv::Point2f object to a std::vector<double>
std::vector<double> Point2fToDoubleVector(std::vector<cv::Point2f> &points) {
  int n = points.size();
  std::vector<double> vec(n);
  for (int i = 0; i < n; i++) {
    vec[i] = (double) points[i].x;
  }
  return vec;
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
double cppSD(std::vector<cv::KeyPoint> &points)
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
  std::vector<double>().swap(inVec);
  return std::sqrt( sd / (n-1) );
}

double cppSD(std::vector<cv::Point2f> &points)
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
  std::vector<double>().swap(inVec);
  return std::sqrt( sd / (n-1) );
}
