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
void Preprocess::loadImages(vector<CustomImage>& images, string inpDir, vector<channelType> channelsOrder)
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

vector<channelType> Preprocess::loadNiceCellImages(vector<CustomImage>& arrayImages,vector<Cell>& cellImages, string inpDir)
{
    vector<string> imageNames;
    vector<vector<Rect>> rectangles;
    vector<vector<channelType>> channels;

    //Oksana's cell
    vector<channelType> channelsOrder1;
    channelsOrder1.push_back(channelType::nucleus);
    channelsOrder1.push_back(channelType::actin);
    channelsOrder1.push_back(channelType::brightfield);

    imageNames.push_back("50x50_niches_RGD_3mM_2022_02_22__14_11_14_Maximum intensity projection_Filter.tif");
    vector<Rect> rVec;
    rVec.push_back(Rect(220, 400, 300, 300));
    rectangles.push_back(rVec);
    rVec.clear();
    channels.push_back(channelsOrder1);

    //Leslie's cells
    vector<channelType> channelsOrder2;
    channelsOrder2.push_back(channelType::nucleus);
    channelsOrder2.push_back(channelType::brightfield);
    channelsOrder2.push_back(channelType::actin);

    imageNames.push_back("300321_PEGvsGELMAlif.lif - Image006-1 (RGB).tif");
    rVec.push_back(Rect(1135, 1200, 270, 480));
    rVec.push_back(Rect(140, 1180, 340, 340));
    rVec.push_back(Rect(1250, 200, 750, 250));
    rVec.push_back(Rect(880, 1630, 350, 250));
    rectangles.push_back(rVec);
    rVec.clear();
    channels.push_back(channelsOrder2);

    imageNames.push_back("300321_PEGvsGELMAlif.lif - Image001-1 (RGB).tif");

    rVec.push_back(Rect(1020, 420, 300, 250));
    rVec.push_back(Rect(10, 1240, 270, 490));
    rVec.push_back(Rect(840, 1080, 270, 260));
    rVec.push_back(Rect(1640, 1130, 310, 240));
    rVec.push_back(Rect(1280, 1090, 220, 310));
    rVec.push_back(Rect(200, 1150, 330, 200));
    rectangles.push_back(rVec);
    rVec.clear();
    channels.push_back(channelsOrder2);

    imageNames.push_back("300321_PEGvsGELMAlif.lif - Image003-1 (RGB).tif");
    rVec.push_back(Rect(330, 600, 280, 350));
    rVec.push_back(Rect(910, 790, 190, 460));
    rectangles.push_back(rVec);
    rVec.clear();
    channels.push_back(channelsOrder2);

    imageNames.push_back("300321_PEGvsGELMAlif.lif - Image004 (RGB).tif");
    rVec.push_back(Rect(1300, 780, 500, 620));
    rectangles.push_back(rVec);
    rVec.clear();
    channels.push_back(channelsOrder2);

    for (int i = 0; i < imageNames.size(); i++)
    {
        Mat img = imread(inpDir + imageNames[i], IMREAD_UNCHANGED);

        Mat bgr[3];
        split(img, bgr);
        vector<Mat> matImages;
        matImages.push_back(bgr[0]);
        matImages.push_back(bgr[1]);
        matImages.push_back(bgr[2]);
        CustomImage image(matImages, channels[i], imageNames[i]);
        arrayImages.push_back(image);

        for (int j = 0; j < rectangles[i].size(); j++)
        {
            if (imageNames[i].find(".tif") != std::string::npos)
            {
                imageNames[i].erase(imageNames[i].find(".tif"));
            }
            string name_cell = imageNames[i] + "_cell" + to_string(j);
            CustomImage ce = image.cutImageOut(rectangles[i][j], name_cell);
            Cell cell(ce);
            cellImages.push_back(cell);
        }
    }
    return channelsOrder2;
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
        summedActinMaxLength = 0, summedYapInNucleus = 0, summedActinFiberAlignment = 0;
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
        summedActinFiberAlignment += cell.actin_fiberAlignment;
    }

    averageAllCells.nucleus_area = summedNucleusArea / cells.size();
    averageAllCells.nucleus_circularity = summedNucleusCircularity / cells.size();
    averageAllCells.nucleus_roundness = summedNucleusRoundness / cells.size();
    averageAllCells.actin_area = summedActinArea / cells.size();
    averageAllCells.actin_density = summedActinDensity / cells.size();
    averageAllCells.actin_maxLength = summedActinMaxLength / cells.size();
    if (!actinPCAangles.empty()) { averageAllCells.actin_mainAngle = static_cast<int>(findDataPointWithMostNeighbours(actinPCAangles, 10)); }
    averageAllCells.yap_inNucleus = summedYapInNucleus / cells.size();
    averageAllCells.actin_fiberAlignment = summedActinFiberAlignment / cells.size();

    return averageAllCells;
}

