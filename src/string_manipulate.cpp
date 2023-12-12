#include <Rcpp.h>
using namespace Rcpp;

// Function to replace a pattern in a string
void replacePattern(std::string &str, const std::string &pattern, const std::string &replacement) {
  size_t pos = 0;
  while ((pos = str.find(pattern, pos)) != std::string::npos) {
    str.replace(pos, pattern.length(), replacement);
    pos += replacement.length();
  }
}

// // Function to replace a pattern in an Rcpp vector of strings
// void replacePatternInRcppVector(Rcpp::CharacterVector &textVector, const std::string &pattern, const std::string &replacement) {
//   int n = textVector.size();
//   for (int i = 0; i < n; ++i) {
//     std::string str = Rcpp::as<std::string>(textVector[i]);
//     replacePattern(str, pattern, replacement);
//     textVector[i] = str;
//   }
// }

// // Function to replace a pattern in an Rcpp vector of strings
// void replacePatternInRcppVector(Rcpp::CharacterVector& textVector, const std::string &pattern, const std::string &replacement) {
//   R_xlen_t n = textVector.size();
//   for (R_xlen_t i = 0; i < n; ++i) {
//     std::string str = Rcpp::as<std::string>(textVector[i]);
//     replacePattern(str, pattern, replacement);
//     textVector[i] = str;
//   }
// }

// [[Rcpp::export]]
Rcpp::CharacterVector replacePatternInRcppVectorWrapper(Rcpp::CharacterVector textVector, const std::string &pattern, const std::string &replacement) {
  // replacePatternInRcppVector(textVector, pattern, replacement);

  R_xlen_t n = textVector.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    std::string str = Rcpp::as<std::string>(textVector[i]);
    replacePattern(str, pattern, replacement);
    textVector[i] = str;
  }

  return textVector;
}
