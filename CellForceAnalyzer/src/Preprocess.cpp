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

void Preprocess::applyYolo(std::vector<CustomImage>& arrayImages, vector<Cell>& cellImages, vector<Mat>& arrayImages_withYoloBoxes, float confThreshold, float nmsThreshold)
{
    for (const auto& entry : arrayImages)
    {
        Mat image_yolo_analysed;
        vector<Cell> detectedCells = m_yolo.detectCells(entry, image_yolo_analysed, confThreshold, nmsThreshold);
        cellImages.insert(cellImages.end(), detectedCells.begin(), detectedCells.end());
        arrayImages_withYoloBoxes.push_back(image_yolo_analysed);
    }
}


// deprecated
void Preprocess::thresholdImages(std::vector<CustomImage>& arrayImages, vector<CustomImage>& cellImages, vector<CustomImage>& arrayImages_thresholded, vector<CustomImage>& cellImages_thresholded)
{
    for (auto& elem : arrayImages)
    {
        arrayImages_thresholded.push_back(elem.getThresholdedImage());
    }

    for (auto& elem : cellImages)
    {
        cellImages_thresholded.push_back(elem.getThresholdedImage());
    }
}


