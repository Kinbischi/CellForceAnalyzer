#pragma once
#include "AnalysisFiberDirection.h"
#include "ParametersFromUI.h"
#include "Cell.h"


class Plotting
{
public:
	Plotting(ParametersFromUI&, AnalysisFiberDirection&, std::vector<Cell>&);

	static std::vector<double> createX(std::vector<double>, bool);
	static std::vector<int> createY(std::vector<double>, std::vector<double>);

	void plot();
	void plotData(std::vector<double>, bool, std::vector<double> = std::vector<double>());

private:
	AnalysisFiberDirection& analysisFiberDir;
	ParametersFromUI& params;
	std::vector<Cell>& cellImages;

};

