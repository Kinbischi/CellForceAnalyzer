#pragma once
#include "CustomImage.h"
class Cell :
    public CustomImage
{
public:
    Cell(CustomImage&);
    Cell() = default;


    double nucleus_circularity;
    double nucleus_roundness;
    int nucleus_area;

    double yap_inNucleus; //in percent

    double actin_density;
    int actin_area;
    double actin_maxLength;
    double actin_PCAangle;

};

