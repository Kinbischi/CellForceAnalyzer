#include "Cell.h"

// whenever an object of a derived class is created, a base class constructor is called
// constructor of cell: all members from CustomImage are initialized using the constructor from CustomImage
Cell::Cell(CustomImage& other) :CustomImage{ other } {};

double Cell::getQuantity(plotType type)
{
	double out;
    switch (type)
    {
    case plotType::nuclArea:
        out = nucleus_area;
        break;
    case plotType::nuclCircularity:
        out = nucleus_circularity;
        break;
    case plotType::nuclRoundness:
        out = nucleus_roundness;
        break;
    case plotType::yapInNucleus:
        out = yap_inNucleus;
        break;
    case plotType::actArea:
        out = actin_area;
        break;
    case plotType::actDensity:
        out = actin_density;
        break;
    case plotType::actMaxLength:
        out = actin_maxLength;
        break;
    case plotType::actMainAngle:
        out = actin_mainAngle;
        break;
    }
	return out;
}

