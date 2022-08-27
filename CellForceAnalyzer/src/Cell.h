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

    //Nucleus
    int nucleus_area;
    double nucleus_avgIntensity;
    double nucleus_circularity;
    double nucleus_roundness;
    
    //YAP
    double yap_percentageInNucleus;
    double yap_avgIntensityInNucleus;
    double yap_avgIntensityOutsideNucleus;

    //Actin
    int actin_area;
    double actin_avgIntensity;
    double actin_maxLength;
    std::vector<double> actin_fibreAnglesPCA; // in degree from 0 to 180 counterclockwise (starting from the right/ x axis)
    double actin_mainAngle;
    double actin_fiberAlignmentValue;

};

