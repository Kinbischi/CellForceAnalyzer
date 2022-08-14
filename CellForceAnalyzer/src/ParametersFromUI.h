#pragma once
#include "helperFunctions.h"
#include "Cell.h"


//TODO: struct?
class ParametersFromUI
{

public:
	// depiction
	int showNumber;
	double scaleFactor;
	bool replacingMode;
	channelType channel;

	bool cellArrays;
	bool singleCells;
	bool deletedCells;
	bool cellArraysWithBoxes;

	thresholdingType threshType;
	displayType dispType;
	plotFeatureType plotFeatType;

	// analysis conducted
	bool fiberPCA;
	bool edgeDetection;
	bool FAdetection;
	bool variousAnalysis;

	// parameters:
	// thresholding
	int manualThreshold;
	// PCA
	int PCAsquareLength;
	double PCAminEigValRatio;
	bool suppressLowEigValRatioSquares;

	// other analysis
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

