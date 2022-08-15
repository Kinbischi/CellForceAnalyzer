#pragma once

struct dataContainer
{
    //vector with order in which image channels are loaded
    std::vector<channelType> channels;

    //data
    std::vector<CustomImage> arrayImages;
    std::vector<cv::Mat> arrayImages_withYoloBoxes;
    std::vector<Cell> cellImages;
    std::vector<Cell> deletedCellImages;

    //Cell object with all the averages stored
    Cell averageAllCells;
};