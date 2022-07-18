#pragma once


#include "Cell.h"
#include "helperFunctions.h"

class Analysis
{

public:

	void analyseNucleus(Cell&);
	void analyseYap(Cell&);
	void analyseActin(Cell&);

	bool isDeadCell(Cell);
	bool analyseShape(cv::Mat&);

	bool pointCloudPCA(const std::vector<cv::Point>&, const int, std::vector<cv::Point2d>&, double = 1.5, std::shared_ptr<double> = nullptr);
	bool analyseWithPCA(cv::Mat&, std::vector<double>&, int, double, cv::Mat = cv::Mat());

	void getPCAoptThresholdedImage(cv::Mat&, int, double);
	int getOptimalThresholdingForPCA(cv::Mat, int, double);

	std::vector<cv::Point> getWhitePointsFromThresholdedImage(cv::Mat);
	std::vector<cv::Point> getPointsDependingOnIntensityFromImage(cv::Mat);

};

