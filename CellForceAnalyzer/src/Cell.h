#pragma once
#include "CustomImage.h"
#include <string>



enum class plotFeatureType { actFibersOptThresh, actFibersIntensity, actFibersBoth, nuclArea, nuclCircularity, nuclRoundness, actArea, actDensity, actMaxLength, actMainAngle, yapInNucleus };

class Cell :
    public CustomImage
{
public:
    Cell(CustomImage&);
    Cell() = default;

    double getQuantity(plotFeatureType);

    int nucleus_area;
    double nucleus_circularity;
    double nucleus_roundness;
    
    double yap_inNucleus; // percentage of the intensity
    
    int actin_area;
    double actin_density;
    double actin_maxLength;
    std::vector<double> actin_fibreAnglesPCA; // in degree from 0 to 180 counterclockwise (starting from the right/ x axis)
    double actin_mainAngle;
    double actin_fiberAlignment;

};

