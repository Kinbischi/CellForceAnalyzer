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
        out = yap_inNucleus;
        break;
    case plotFeatureType::actArea:
        out = actin_area;
        break;
    case plotFeatureType::actDensity:
        out = actin_density;
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

string Cell::getPlotTitle(plotFeatureType CFType)
{
    string name;
    switch (CFType)
    {
    case plotFeatureType::actFibersOptThresh:
        name = "Actin fiber analysis (optimally thresholded)";
        break;
    case plotFeatureType::actFibersIntensity:
        name = "Actin fiber analysis (intensity-based)";
        break;
    case plotFeatureType::actFibersBoth:
        name = "Actin fiber analysis";
        break;
    case plotFeatureType::nuclArea:
        name = "Nucleus area";
        break;
    case plotFeatureType::nuclCircularity:
        name = "Nucleus circularity";
        break;
    case plotFeatureType::nuclRoundness:
        name = "Nucleus roundness";
        break;
    case plotFeatureType::yapInNucleus:
        name = "YAP in nucleus";
        break;
    case plotFeatureType::actArea:
        name = "Actin area";
        break;
    case plotFeatureType::actDensity:
        name = "Actin density";
        break;
    case plotFeatureType::actMaxLength:
        name = "Actin maximum stretch";
        break;
    case plotFeatureType::actMainAngle:
        name = "Actin main angle";
        break;
    }
    return name;
}

string Cell::getPlotXLabel(plotFeatureType CFType)
{
    string name;
    switch (CFType)
    {
    case plotFeatureType::actFibersOptThresh:
        name = "Angle (in degrees)";
        break;
    case plotFeatureType::actFibersIntensity:
        name = "Angle (in degrees)";
        break;
    case plotFeatureType::actFibersBoth:
        name = "Angle (in degrees)";
        break;
    default:
        name = getPlotTitle(CFType);
        break;
    }
    return name;
}

string Cell::getPlotYLabel(plotFeatureType CFType)
{
    string name;
    switch (CFType)
    {
    case plotFeatureType::actFibersOptThresh:
        name = "Number of squares";
        break;
    case plotFeatureType::actFibersIntensity:
        name = "Number of squares";
        break;
    case plotFeatureType::actFibersBoth:
        name = "Number of squares";
        break;
    default:
        name = "Number of occurrences";
        break;
    }
    return name;
}

vector<string> Cell::getPlotLegendNames(plotFeatureType CFType)
{
    vector<string> names(2);
    switch (CFType)
    {
    case plotFeatureType::actFibersBoth:
        names[0] = "Optimally thresholded image PCA";
        names[1] = "Intensity-based PCA";
        break;
   
    }
    return names;
}

