#pragma once


#include "Cell.h"

class Analysis
{

public:

	void analyseNucleusShape(Cell&);
	void analyseYapInNucleus(Cell&);
	bool isDeadCell(Cell);
	void analyseActin(Cell&);
	Cell getAverageProperties(std::vector<Cell>);
	bool analyseShape(cv::Mat&);

	bool pointCloudPCA(const std::vector<cv::Point>&, const int, std::vector<cv::Point2d>&);
	bool analyseWithPCA(cv::Mat&, std::vector<double>&);
	bool analyseWithPCA2(cv::Mat&, std::vector<double>&, int&);

	std::vector<cv::Point> getWhitePointsFromThresholdedImage(cv::Mat);

};

