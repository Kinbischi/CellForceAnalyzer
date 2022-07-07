#pragma once


#include "Cell.h"

class Analysis
{

public:

	void analyseNucleus(Cell&);
	void analyseYap(Cell&);
	void analyseActin(Cell&);

	bool isDeadCell(Cell);
	bool analyseShape(cv::Mat&);

	bool pointCloudPCA(const std::vector<cv::Point>&, const int, std::vector<cv::Point2d>&, double&, double = 1.5);
	bool analyseWithPCA(cv::Mat&, std::vector<double>&, int = 15, double = 1.5, cv::Mat = cv::Mat());
	int getOptimalThresholdingForPCA(cv::Mat, int, double);

	std::vector<cv::Point> getWhitePointsFromThresholdedImage(cv::Mat);
	std::vector<cv::Point> getPointsDependingOnIntensityFromImage(cv::Mat);

};

