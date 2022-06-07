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

	double getPCAorientation(const std::vector<cv::Point>&, std::vector<cv::Point>&);

	bool getPCAorientation2(const std::vector<cv::Point>&, std::vector<cv::Point2d>&, std::vector<double>&, cv::Point&);
	void drawAxis(cv::Mat&, cv::Point, cv::Point, cv::Scalar, const float);

	void drawAxis2(cv::Mat&, cv::Point, std::vector<cv::Point2d>, std::vector<double>, cv::Scalar, int=2, const float = 0.2);
	std::vector<cv::Point> getWhitePointsFromThresholdedImage(cv::Mat);

};

