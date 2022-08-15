#pragma once

#include "ParametersUI.h"

class AnalysisFiberDirection
{
public:
	AnalysisFiberDirection(ParametersUI&);
	
	bool pointCloudPCA(const std::vector<cv::Point>&, std::vector<cv::Point2d>&, std::shared_ptr<double> = nullptr);
	bool analyseWithPCA(cv::Mat&, std::vector<double>&, cv::Mat = cv::Mat());

	void getPCAoptThresholdedImage(cv::Mat&);
	double calculateFiberAlignmentConstant(std::vector<int>, double, int);
	std::vector<double> testYvecs(std::vector<std::vector<int>>, double, int);
private:
	ParametersUI& params;
};

