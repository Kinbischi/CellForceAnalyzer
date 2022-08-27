#include "Cell.h"

// whenever an object of a derived class is created, a base class constructor is called
// constructor of cell: all members from CustomImage are initialized using the constructor from CustomImage

using namespace std;

Cell::Cell(CustomImage& other) :CustomImage{ other } {};

double Cell::getQuantity(plotFeatureType type)
{
	double out;
    switch (type)
    {
    case plotFeatureType::nuclArea:
        out = nucleus_area;
        break;
    case plotFeatureType::nuclCircularity:
        out = nucleus_circularity;
        break;
    case plotFeatureType::nuclRoundness:
        out = nucleus_roundness;
        break;
    case plotFeatureType::yapInNucleus:
        out = yap_percentageInNucleus;
        break;
    case plotFeatureType::actArea:
        out = actin_area;
        break;
    case plotFeatureType::actDensity:
        out = actin_avgIntensity;
        break;
    case plotFeatureType::actMaxLength:
        out = actin_maxLength;
        break;
    case plotFeatureType::actMainAngle:
        out = actin_mainAngle;
        break;
    }
	return out;
}



