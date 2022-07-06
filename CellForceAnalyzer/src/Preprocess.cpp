#include "Preprocess.h"

#include <filesystem>
#include <string>
#include <set>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "CustomImage.h"
#include "helperFunctions.h"

using namespace std;
using namespace cv;
namespace fs = std::filesystem;

// if there is numbers in the start of the image => can happen that images are not loaded in the same order as in the windows folder
void Preprocess::loadImages(std::vector<CustomImage>& images, string inpDir, std::vector<channelType> channelsOrder)
{
    fs::path inPath = inpDir;
    vector<Mat> matImages;
    int i = 0;
    for (const auto& entry : fs::directory_iterator(inPath))
    {
        string imgPath = entry.path().string();
        Mat img = imread(imgPath, IMREAD_UNCHANGED);
        auto dep = img.depth(); //if 0 => 8bit image(uchar), if 2 => 16 bit image(ushort)

        if (img.channels()==3) //if images are loaded as rgb images and the channel that we are interested is in only one of the rgb channels
        {
            Mat bgr[3];
            split(img, bgr);

            int maxNonZeros=0;
            for (int i = 0; i < 3; i++)
            {
                Mat channel = bgr[i];
                int nonZeroPixels = countNonZero(channel);
                if (nonZeroPixels>maxNonZeros)
                {
                    maxNonZeros = nonZeroPixels;
                    img = channel;
                }
            }
        }
        matImages.push_back(img);

        if (i == channelsOrder.size() - 1)
        {
            //whenever all of the (e.g. 4) image channels have been loaded, a CustomImage is created
            string imgName = entry.path().filename().string();
            CustomImage image(matImages, channelsOrder, imgName);

            images.push_back(image);

            i = -1;
            matImages.clear();
        }
        i++;
    }
}

void Preprocess::applyYolo(std::vector<CustomImage>& arrayImages, vector<Cell>& cellImages, vector<Mat>& arrayImages_withYoloBoxes, float confThreshold)
{
    for (const auto& entry : arrayImages)
    {
        Mat image_yolo_analysed;
        vector<Cell> detectedCells = m_yolo.detectCells(entry, image_yolo_analysed, confThreshold);
        cellImages.insert(cellImages.end(), detectedCells.begin(), detectedCells.end());
        arrayImages_withYoloBoxes.push_back(image_yolo_analysed);
    }
}

//better name?
double findDataPointWithMostNeighbours(vector<double> input, double margin)
{
    int i = 0;
    vector<double> neighboursCount(input.size());
    for (auto elem : input)
    {
        for (auto neighbour : input)
        {
            if (neighbour<elem + margin && neighbour>elem - margin)
            {
                neighboursCount[i]++;
            }
        }
        i++;
    }

    vector<double> results;

    double maxNeighbours = *max_element(neighboursCount.begin(), neighboursCount.end());
    for (int j = 0; j < neighboursCount.size(); j++)
    {
        if (neighboursCount[j] == maxNeighbours)
        {
            results.push_back(input[j]);
        }
    }

    return help::median(results);
}

Cell Preprocess::getAverageProperties(const vector<Cell> cells)
{
    double summedNucleusArea = 0, summedNucleusCircularity = 0, summedNucleusRoundness = 0, summedActinArea = 0, summedActinDensity = 0,
        summedActinMaxLength = 0, summedYapInNucleus = 0;
    vector<double> actinPCAangles;
    Cell averageAllCells;
    for (auto cell : cells)
    {
        summedNucleusArea += cell.nucleus_area;
        summedNucleusCircularity += cell.nucleus_circularity;
        summedNucleusRoundness += cell.nucleus_roundness;

        summedActinArea += cell.actin_area;
        summedActinDensity += cell.actin_density;
        summedActinMaxLength += cell.actin_maxLength;
        actinPCAangles.push_back(cell.actin_mainAngle);

        summedYapInNucleus += cell.yap_inNucleus;
    }

    averageAllCells.nucleus_area = summedNucleusArea / cells.size();
    averageAllCells.nucleus_circularity = summedNucleusCircularity / cells.size();
    averageAllCells.nucleus_roundness = summedNucleusRoundness / cells.size();
    averageAllCells.actin_area = summedActinArea / cells.size();
    averageAllCells.actin_density = summedActinDensity / cells.size();
    averageAllCells.actin_maxLength = summedActinMaxLength / cells.size();
    if (!actinPCAangles.empty()) { averageAllCells.actin_mainAngle = static_cast<int>(findDataPointWithMostNeighbours(actinPCAangles, 10)); }
    averageAllCells.yap_inNucleus = summedYapInNucleus / cells.size();

    return averageAllCells;
}

