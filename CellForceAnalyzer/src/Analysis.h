#pragma once


#include "Cell.h"
#include "ParametersFromUI.h"
#include "helperFunctions.h"

class Analysis
{

public:
	Analysis(ParametersFromUI&);

	void analyseCell(Cell&, std::vector<channelType>);

	bool isDeadCell(Cell);
	bool analyseShape(cv::Mat&);

	void edgeDetectionCanny(cv::Mat&);
	void focalAdhesiondetection(cv::Mat&);

	bool pointCloudPCA(const std::vector<cv::Point>&, std::vector<cv::Point2d>&, std::shared_ptr<double> = nullptr);
	bool analyseWithPCA(cv::Mat&, std::vector<double>&, cv::Mat = cv::Mat());

	std::vector<int> createY(std::vector<double>, std::vector<double>);
	std::vector<double> createX(std::vector<double>, bool);

	void getPCAoptThresholdedImage(cv::Mat&);
	double calculateFiberAlignmentConstant(std::vector<int>, double,int);
	std::vector<double> Analysis::testYvecs(std::vector<std::vector<int>>, double, int);

private:

	void analyseActin(Cell&);

	ParametersFromUI& m_parameters;
};

