#pragma once
#include "CustomImage.h"

enum class plotType { actFibers, nuclArea, nuclCircularity, nuclRoundness, yapInNucleus, actArea, actDensity, actMaxLength, actMainAngle};

class Cell :
    public CustomImage
{
public:
    Cell(CustomImage&);
    Cell() = default;

    double getQuantity(plotType);

    int nucleus_area;
    double nucleus_circularity;
    double nucleus_roundness;
    
    double yap_inNucleus; //in percent
    
    int actin_area;
    double actin_density;
    double actin_maxLength;
    std::vector<double> actin_fibreAnglesPCA; // in degree from 0 to 180 counterclockwise (starting from the right/ x axis)
    double actin_mainAngle;

};

