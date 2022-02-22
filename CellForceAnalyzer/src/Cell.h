#pragma once
#include "CustomImage.h"
class Cell :
    public CustomImage
{
public:
    Cell(CustomImage&);
    
    double nucleus_circularity;
    double nucleus_roundness;
    int nucleus_area;
    double yapInNucleus; //in percent

};

