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

private:

	void analyseActin(Cell&);

	ParametersFromUI& m_parameters;
};

