#pragma once
#include "AnalysisFiberDirection.h"
#include "ParametersUI.h"
#include "Cell.h"
#include "dataContainer.h"


class Plotting
{
public:

	Plotting(ParametersUI&, dataContainer&, AnalysisFiberDirection&);

	std::string getPlotTitle(plotFeatureType);
	std::vector<std::string> getPlotLegendNames(plotFeatureType);
	std::string getPlotXLabel(plotFeatureType);
	std::string getPlotYLabel(plotFeatureType);

	static std::vector<double> createX(std::vector<double>, bool);
	static std::vector<int> createY(std::vector<double>, std::vector<double>);

	int getCellImageToPlot(cv::Mat&);
	void plotSomething(std::vector<double>, std::vector<double>, plotFeatureType, std::string);
	int plot();
	void plotData(std::vector<double>, bool, std::vector<double> = std::vector<double>());

private:

	ParametersUI& params;
	dataContainer& data;
	AnalysisFiberDirection& analysisFiberDir;

};

