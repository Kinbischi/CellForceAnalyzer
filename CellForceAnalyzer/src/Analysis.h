#pragma once


#include "Cell.h"
#include "ParametersUI.h"
#include "helperFunctions.h"

class Analysis
{

public:
	Analysis(ParametersUI&);

	void analyseCell(Cell&, std::vector<channelType>);

	bool isDeadCell(Cell);
	bool analyseShape(cv::Mat&);

	void edgeDetectionCanny(cv::Mat&);
	void focalAdhesiondetection(cv::Mat&);

private:

	void analyseActin(Cell&);

	ParametersUI& m_parameters;
};

