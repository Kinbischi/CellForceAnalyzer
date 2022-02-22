#pragma once

#include <string>
#include <opencv2/dnn.hpp>
#include <opencv2/core.hpp>

#include "CustomImage.h"
#include "Cell.h"

class YoloNicheDetector
{
public:
	YoloNicheDetector();
	std::vector<Cell> detectCells(CustomImage, cv::Mat&,float,float);

private:
	cv::dnn::Net m_net;
	std::vector<std::string> m_classes;
	std::string m_yoloInputDir = "YoloFiles/";
};