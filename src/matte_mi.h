#ifndef VOLTRON_MATTE_MI_H
#define VOLTRON_MATTE_MI_H

#include <opencv2/core.hpp>

cv::Mat1d MatteMIMap(
    const cv::Mat& fixed,
    const cv::Mat& moving,
    const cv::Mat& mask,
    int bins = 50);

cv::Mat1d chunkedMatteMIMap(
    const cv::Mat& fixed,
    const cv::Mat& moving,
    const cv::Mat& mask,
    int bins = 50);

double MatteMI(
    const cv::Mat& fixed,
    const cv::Mat& moving,
    const cv::Mat& mask,
    int bins = 50);

#endif