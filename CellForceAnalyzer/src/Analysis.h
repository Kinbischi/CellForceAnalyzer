#pragma once


#include "Cell.h"
#include "ParametersUI.h"
#include "helperFunctions.h"
#include "AnalysisFiberDirection.h"

class Analysis
{

public:
	Analysis(ParametersUI&, AnalysisFiberDirection&);
	static void blurWeak(cv::Mat&);
	static void blurStrong(cv::Mat&);
	static bool thresh(cv::Mat&, int=0, int=2, bool=true);

	bool pointCloudPCA2(const std::vector<cv::Point>&, std::vector<cv::Point2d>&, std::shared_ptr<double> = nullptr);

	void analyseCell(Cell&, std::vector<channelType>);

	bool isDeadCell(Cell);
	bool analyseShape(cv::Mat&);

	void edgeDetectionCanny(cv::Mat&);
	void focalAdhesiondetection(cv::Mat&);

private:

	void analyseActin(Cell&);

	ParametersUI& params;
	AnalysisFiberDirection& analysisFiberDir;
};

