#pragma once


#include "Cell.h"

class Analysis
{

public:

	void analyseNucleusShape(Cell&);
	void analyseYapInNucleus(Cell&);
	std::vector<int> removeBadCells(std::vector<Cell>&);
	void analyseActin(Cell&);
	double findDataPointWithMostNeighbours(std::vector<double>, double);
	bool analyseShape(cv::Mat&);

};

