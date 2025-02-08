// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// automated_registeration_rawvector
Rcpp::List automated_registeration_rawvector(Rcpp::RawVector ref_image, Rcpp::RawVector query_image, const int width1, const int height1, const int width2, const int height2, const float GOOD_MATCH_PERCENT, const int MAX_FEATURES, const bool invert_query, const bool invert_ref, Rcpp::String flipflop_query, Rcpp::String flipflop_ref, Rcpp::String rotate_query, Rcpp::String rotate_ref, Rcpp::String matcher, Rcpp::String method);
RcppExport SEXP _VoltRon_automated_registeration_rawvector(SEXP ref_imageSEXP, SEXP query_imageSEXP, SEXP width1SEXP, SEXP height1SEXP, SEXP width2SEXP, SEXP height2SEXP, SEXP GOOD_MATCH_PERCENTSEXP, SEXP MAX_FEATURESSEXP, SEXP invert_querySEXP, SEXP invert_refSEXP, SEXP flipflop_querySEXP, SEXP flipflop_refSEXP, SEXP rotate_querySEXP, SEXP rotate_refSEXP, SEXP matcherSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type ref_image(ref_imageSEXP);
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type query_image(query_imageSEXP);
    Rcpp::traits::input_parameter< const int >::type width1(width1SEXP);
    Rcpp::traits::input_parameter< const int >::type height1(height1SEXP);
    Rcpp::traits::input_parameter< const int >::type width2(width2SEXP);
    Rcpp::traits::input_parameter< const int >::type height2(height2SEXP);
    Rcpp::traits::input_parameter< const float >::type GOOD_MATCH_PERCENT(GOOD_MATCH_PERCENTSEXP);
    Rcpp::traits::input_parameter< const int >::type MAX_FEATURES(MAX_FEATURESSEXP);
    Rcpp::traits::input_parameter< const bool >::type invert_query(invert_querySEXP);
    Rcpp::traits::input_parameter< const bool >::type invert_ref(invert_refSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type flipflop_query(flipflop_querySEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type flipflop_ref(flipflop_refSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rotate_query(rotate_querySEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type rotate_ref(rotate_refSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type matcher(matcherSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(automated_registeration_rawvector(ref_image, query_image, width1, height1, width2, height2, GOOD_MATCH_PERCENT, MAX_FEATURES, invert_query, invert_ref, flipflop_query, flipflop_ref, rotate_query, rotate_ref, matcher, method));
    return rcpp_result_gen;
END_RCPP
}
// replaceNaMatrix
Rcpp::NumericMatrix replaceNaMatrix(Rcpp::NumericMatrix mat, int replace);
RcppExport SEXP _VoltRon_replaceNaMatrix(SEXP matSEXP, SEXP replaceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type replace(replaceSEXP);
    rcpp_result_gen = Rcpp::wrap(replaceNaMatrix(mat, replace));
    return rcpp_result_gen;
END_RCPP
}
// warpImage
Rcpp::RawVector warpImage(Rcpp::RawVector ref_image, Rcpp::RawVector query_image, Rcpp::List mapping, const int width1, const int height1, const int width2, const int height2);
RcppExport SEXP _VoltRon_warpImage(SEXP ref_imageSEXP, SEXP query_imageSEXP, SEXP mappingSEXP, SEXP width1SEXP, SEXP height1SEXP, SEXP width2SEXP, SEXP height2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type ref_image(ref_imageSEXP);
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type query_image(query_imageSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mapping(mappingSEXP);
    Rcpp::traits::input_parameter< const int >::type width1(width1SEXP);
    Rcpp::traits::input_parameter< const int >::type height1(height1SEXP);
    Rcpp::traits::input_parameter< const int >::type width2(width2SEXP);
    Rcpp::traits::input_parameter< const int >::type height2(height2SEXP);
    rcpp_result_gen = Rcpp::wrap(warpImage(ref_image, query_image, mapping, width1, height1, width2, height2));
    return rcpp_result_gen;
END_RCPP
}
// warpImageAuto
Rcpp::RawVector warpImageAuto(Rcpp::RawVector ref_image, Rcpp::RawVector query_image, Rcpp::List mapping, const int width1, const int height1, const int width2, const int height2);
RcppExport SEXP _VoltRon_warpImageAuto(SEXP ref_imageSEXP, SEXP query_imageSEXP, SEXP mappingSEXP, SEXP width1SEXP, SEXP height1SEXP, SEXP width2SEXP, SEXP height2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type ref_image(ref_imageSEXP);
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type query_image(query_imageSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mapping(mappingSEXP);
    Rcpp::traits::input_parameter< const int >::type width1(width1SEXP);
    Rcpp::traits::input_parameter< const int >::type height1(height1SEXP);
    Rcpp::traits::input_parameter< const int >::type width2(width2SEXP);
    Rcpp::traits::input_parameter< const int >::type height2(height2SEXP);
    rcpp_result_gen = Rcpp::wrap(warpImageAuto(ref_image, query_image, mapping, width1, height1, width2, height2));
    return rcpp_result_gen;
END_RCPP
}
// warpImageManual
Rcpp::RawVector warpImageManual(Rcpp::RawVector ref_image, Rcpp::RawVector query_image, Rcpp::List mapping, const int width1, const int height1, const int width2, const int height2);
RcppExport SEXP _VoltRon_warpImageManual(SEXP ref_imageSEXP, SEXP query_imageSEXP, SEXP mappingSEXP, SEXP width1SEXP, SEXP height1SEXP, SEXP width2SEXP, SEXP height2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type ref_image(ref_imageSEXP);
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type query_image(query_imageSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mapping(mappingSEXP);
    Rcpp::traits::input_parameter< const int >::type width1(width1SEXP);
    Rcpp::traits::input_parameter< const int >::type height1(height1SEXP);
    Rcpp::traits::input_parameter< const int >::type width2(width2SEXP);
    Rcpp::traits::input_parameter< const int >::type height2(height2SEXP);
    rcpp_result_gen = Rcpp::wrap(warpImageManual(ref_image, query_image, mapping, width1, height1, width2, height2));
    return rcpp_result_gen;
END_RCPP
}
// manual_registeration_rawvector
Rcpp::List manual_registeration_rawvector(Rcpp::RawVector ref_image, Rcpp::RawVector query_image, Rcpp::NumericMatrix reference_landmark, Rcpp::NumericMatrix query_landmark, const int width1, const int height1, const int width2, const int height2, Rcpp::String method);
RcppExport SEXP _VoltRon_manual_registeration_rawvector(SEXP ref_imageSEXP, SEXP query_imageSEXP, SEXP reference_landmarkSEXP, SEXP query_landmarkSEXP, SEXP width1SEXP, SEXP height1SEXP, SEXP width2SEXP, SEXP height2SEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type ref_image(ref_imageSEXP);
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type query_image(query_imageSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type reference_landmark(reference_landmarkSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type query_landmark(query_landmarkSEXP);
    Rcpp::traits::input_parameter< const int >::type width1(width1SEXP);
    Rcpp::traits::input_parameter< const int >::type height1(height1SEXP);
    Rcpp::traits::input_parameter< const int >::type width2(width2SEXP);
    Rcpp::traits::input_parameter< const int >::type height2(height2SEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(manual_registeration_rawvector(ref_image, query_image, reference_landmark, query_landmark, width1, height1, width2, height2, method));
    return rcpp_result_gen;
END_RCPP
}
// applyMapping
Rcpp::NumericMatrix applyMapping(Rcpp::NumericMatrix coords, Rcpp::List mapping);
RcppExport SEXP _VoltRon_applyMapping(SEXP coordsSEXP, SEXP mappingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mapping(mappingSEXP);
    rcpp_result_gen = Rcpp::wrap(applyMapping(coords, mapping));
    return rcpp_result_gen;
END_RCPP
}
// build_snn_rank
Rcpp::List build_snn_rank(Rcpp::IntegerMatrix neighbors);
RcppExport SEXP _VoltRon_build_snn_rank(SEXP neighborsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type neighbors(neighborsSEXP);
    rcpp_result_gen = Rcpp::wrap(build_snn_rank(neighbors));
    return rcpp_result_gen;
END_RCPP
}
// build_snn_number
Rcpp::List build_snn_number(Rcpp::IntegerMatrix neighbors);
RcppExport SEXP _VoltRon_build_snn_number(SEXP neighborsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type neighbors(neighborsSEXP);
    rcpp_result_gen = Rcpp::wrap(build_snn_number(neighbors));
    return rcpp_result_gen;
END_RCPP
}
// replacePatternInRcppVectorWrapper
Rcpp::CharacterVector replacePatternInRcppVectorWrapper(Rcpp::CharacterVector textVector, const std::string& pattern, const std::string& replacement);
RcppExport SEXP _VoltRon_replacePatternInRcppVectorWrapper(SEXP textVectorSEXP, SEXP patternSEXP, SEXP replacementSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type textVector(textVectorSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type pattern(patternSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type replacement(replacementSEXP);
    rcpp_result_gen = Rcpp::wrap(replacePatternInRcppVectorWrapper(textVector, pattern, replacement));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VoltRon_automated_registeration_rawvector", (DL_FUNC) &_VoltRon_automated_registeration_rawvector, 16},
    {"_VoltRon_replaceNaMatrix", (DL_FUNC) &_VoltRon_replaceNaMatrix, 2},
    {"_VoltRon_warpImage", (DL_FUNC) &_VoltRon_warpImage, 7},
    {"_VoltRon_warpImageAuto", (DL_FUNC) &_VoltRon_warpImageAuto, 7},
    {"_VoltRon_warpImageManual", (DL_FUNC) &_VoltRon_warpImageManual, 7},
    {"_VoltRon_manual_registeration_rawvector", (DL_FUNC) &_VoltRon_manual_registeration_rawvector, 9},
    {"_VoltRon_applyMapping", (DL_FUNC) &_VoltRon_applyMapping, 2},
    {"_VoltRon_build_snn_rank", (DL_FUNC) &_VoltRon_build_snn_rank, 1},
    {"_VoltRon_build_snn_number", (DL_FUNC) &_VoltRon_build_snn_number, 1},
    {"_VoltRon_replacePatternInRcppVectorWrapper", (DL_FUNC) &_VoltRon_replacePatternInRcppVectorWrapper, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_VoltRon(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
