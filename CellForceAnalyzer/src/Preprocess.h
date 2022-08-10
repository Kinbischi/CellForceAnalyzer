#pragma once
#include <vector>
#include <string>

#include "CustomImage.h"
#include "YoloNicheDetector.h"
#include "Cell.h"

class Preprocess
{
public:

	void loadImages(std::vector<CustomImage>&, std::string,std::vector<channelType>);

	std::vector<channelType> loadNiceCellImages(std::vector<CustomImage>&, std::vector<Cell>&, std::string);

	void applyYolo(std::vector<CustomImage>&, std::vector<Cell>&, std::vector<cv::Mat>&, float);

	Cell getAverageProperties(const std::vector<Cell>);

private:
	YoloNicheDetector m_yolo;
};

