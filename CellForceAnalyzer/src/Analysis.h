#pragma once


#include "Cell.h"

class Analysis
{

public:

	void analyseNucleusShape(Cell&);
	void analyseYapInNucleus(Cell&);
	void removeBadCells(std::vector<Cell>&);
	void analyseActin(Cell&);
	double findPointWithMostNeighbours(std::vector<double>, double);
	bool analyseShape(cv::Mat&);

};

