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

	void applyYolo(std::vector<CustomImage>&, std::vector<Cell>&, std::vector<cv::Mat>&, float, float);

	void thresholdImages(std::vector<CustomImage>&, std::vector<CustomImage>&, std::vector<CustomImage>&, std::vector<CustomImage>&);

	bool analyseShape(cv::Mat&);
	void analyseNucleusShape(Cell& cell);
	void yapDistribution(Cell&);

private:
	YoloNicheDetector m_yolo;
};

