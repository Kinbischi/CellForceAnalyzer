#pragma once
#include "helperFunctions.h"
#include "Cell.h"

class ParametersFromUI
{

public:
	int showNumber;

	thresholdingType threshType;
	displayType dispType;
	plotFeatureType plotFeatType;

	int PCAsquareLength;
	double PCAminEigValRatio;
	bool suppressLowEigValRatioSquares;

	double highThreshCanny;
	double FAminDist;
	double FAminCircleConfidence;
	double FAdp;

	bool withAnalysis;
	
	//default values for double plots
	int goodOptPCAsquareLength = 8;
	int goodIntensityPCAsquareLength = 8;
	double goodOptPCAminEigValRatio = 7;
	double goodIntensityPCAminEigValRatio = 1.3;
};

