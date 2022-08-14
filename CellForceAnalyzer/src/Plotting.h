#pragma once
#include "AnalysisFiberDirection.h"
#include "ParametersFromUI.h"
#include "Cell.h"


class Plotting
{
public:
	Plotting(ParametersFromUI&, AnalysisFiberDirection&, std::vector<Cell>&);

	std::string getPlotTitle(plotFeatureType);
	std::vector<std::string> getPlotLegendNames(plotFeatureType);
	std::string getPlotXLabel(plotFeatureType);
	std::string getPlotYLabel(plotFeatureType);

	static std::vector<double> createX(std::vector<double>, bool);
	static std::vector<int> createY(std::vector<double>, std::vector<double>);

	//int getImageToPlot(cv::Mat&);
	void plotSomething(std::vector<double>, std::vector<double>, plotFeatureType, std::string);
	void plot();
	void plotData(std::vector<double>, bool, std::vector<double> = std::vector<double>());


private:
	AnalysisFiberDirection& analysisFiberDir;
	ParametersFromUI& params;
	std::vector<Cell>& cellImages;

};

