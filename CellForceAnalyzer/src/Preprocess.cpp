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
        

        if (img.channels()==3) //if channels are in one of the rgb channels
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


void getRadii(vector<Point> contour, Point centroid, vector<Point>& outPoints, double& innerRadius, double& outerRadius)
{
    Point maxP;
    Point minP;
    outerRadius = 0;
    innerRadius = 10000; //random high radius that gets replaced

    for (auto point : contour)
    {
        double distance = cv::norm(centroid - point);
        if (distance > outerRadius)
        {
            outerRadius = distance;
            maxP = point;
        }
        if (distance < innerRadius)
        {
            innerRadius = distance;
            minP = point;
        }
    }
    outPoints[0] = minP;
    outPoints[1] = maxP;
}

vector<Point> getLongestContour(vector<vector<Point> > contours)
{
    vector<Point> maxContour;
    int longestSize = 0;

    for (auto elem : contours)
    {
        if (elem.size() > longestSize)
        {
            maxContour = elem;
            longestSize = elem.size();
        }
    }
    return maxContour;
}

void getContourAndCentroidOfLargestWhiteRegion(Mat img, vector<vector<Point> > contours, vector<Point>& longestContour, Point& centroidLargest)
{
    Mat labels, stats, centroids;
    int numberOfLabels = cv::connectedComponentsWithStats(img, labels, stats, centroids);

    //get label of Largest White region (white: 255 => in Binary images: every nonZero)
    Mat nonZeros;
    cv::findNonZero(img, nonZeros);

    MatIterator_<Point> it, end;
    set<int> labelList;
    for (it = nonZeros.begin<Point>(), end = nonZeros.end<Point>(); it != end; ++it)
    {
        Point p = *it;
        int label = labels.at<int>(p);
        labelList.insert(label);
    }

    int maxLabel = 0;
    int maxArea = 0;
    for (int label : labelList)
    {
        int area = stats.at<int>(label, CC_STAT_AREA);
        if (area > maxArea)
        {
            maxLabel = label;
            maxArea = area;
        }
    }

    // centroid from Largest white region
    centroidLargest = Point(centroids.at<double>(maxLabel, 0), centroids.at<double>(maxLabel, 1));
    
    // get longest contour
    int longestSize = 0;
    for (auto elem : contours)
    {
        if (elem.size() > longestSize)
        {
            longestContour = elem;
            longestSize = elem.size();
        }
    }

    //the contours are around the white objects => it is assumed that the longest contour contains the largest white area 
    //=> could lead to wrong results in roundness for specific cases! (e.g. for two objects with similar area but different contour)
    //=> contour of one object is taken and centroid of the other
    //no optimal opencv function that connects contours and enclosed areas and their centroids
    //=>countourArea and CC_STAT_AREA give different results and can not be compared!

    //=>small test to counter that (checks whether centroid is within the largest contour)
    double isInContour = cv::pointPolygonTest(longestContour, centroidLargest, false);
    assert(isInContour + 1);
}

// function used for depiction in show image function, draws centroid, inner & outer radii,...
bool Preprocess::analyseShape(Mat& img)
{
    bool successful = help::thresh(img);
    if (!successful)
    {
        return false;
    }

    vector<vector<Point> > contours;
    cv::findContours(img, contours, RETR_LIST, CHAIN_APPROX_NONE);

    Point centroidLargest;
    vector<Point> longestContour;
    getContourAndCentroidOfLargestWhiteRegion(img, contours, longestContour, centroidLargest);

    vector<Point> radiusPoints(2);
    double innerRadius, outerRadius;
    getRadii(longestContour, centroidLargest, radiusPoints, innerRadius, outerRadius);

    cv::cvtColor(img, img, COLOR_GRAY2RGB);
    cv::drawContours(img, contours, -1, Scalar(0, 255, 0));
    cv::drawMarker(img, centroidLargest, Scalar(255, 0, 0), 0, 5);
    cv::drawMarker(img, radiusPoints[0], Scalar(0, 0, 255), 0, 5);
    cv::drawMarker(img, radiusPoints[1], Scalar(0, 0, 255), 0, 5);

    return true;
}

// function used to get circularity and roundness for every cell
void Preprocess::analyseNucleusShape(Cell& cell)
{
    Mat img = cell.m_nucleusChannel.clone();
    help::thresh(img);

    vector<vector<Point> > contours;
    cv::findContours(img, contours, RETR_LIST, CHAIN_APPROX_NONE);

    Point centroidLargest;
    vector<Point> longestContour;
    getContourAndCentroidOfLargestWhiteRegion(img, contours, longestContour, centroidLargest);

    //circularity
    double perimeter = cv::arcLength(longestContour, true);
    double area = cv::contourArea(longestContour); // area inside contour
    double circularity = 4 * help::M_PI * area / (perimeter * perimeter);

    //roundness
    vector<Point> radiusPoints(2);
    double innerRadius, outerRadius;
    getRadii(longestContour, centroidLargest, radiusPoints, innerRadius, outerRadius);
    double roundness = innerRadius / outerRadius;

    cell.nucleus_circularity = circularity;
    cell.nucleus_roundness = roundness;
    cell.nucleus_area = area;
}

//TODO: normalize cells seperately, threshold cells again and store them in cell obj?? => used multiple times
//get largest connected area only??
//assert cells are good??


double analyzeYap(Mat yap, Mat nucleus)
{
    auto my = help::elementsMap<uchar>(yap);
    auto mn = help::elementsMap<uchar>(nucleus);

    //assert is same dimens
    int yapInNucleusCount = 0;
    int yapOutsideNucleusCount = 0;
    int yapArea = 0;
    for (int i = 0; i < yap.rows; i++)
    {
        for (int j = 0; j < yap.cols; j++)
        {
            uchar n = nucleus.at<uchar>(i,j);
            uchar y = yap.at<uchar>(i,j);
            
            if (y == CustomImage::max8bit)
            {
                yapArea++;
                if (n == CustomImage::max8bit)
                {
                    yapInNucleusCount++;
                }
                if (n == 0)
                {
                    yapOutsideNucleusCount++;
                }
            }
        }
    }

    double percentageOfYapInNucleus = double(yapInNucleusCount) / yapArea;
    double percentageOfYapOutsideNucleus = double(yapOutsideNucleusCount) / yapArea;
    return percentageOfYapInNucleus;
}

void Preprocess::yapDistribution(Cell& cell)
{
    Mat yap = cell.m_yapChannel.clone();
    help::thresh(yap);
    Mat nucleus = cell.m_nucleusChannel.clone();
    help::thresh(nucleus);
    double yapInNucleus = analyzeYap(yap, nucleus);
    cell.yapInNucleus = yapInNucleus;
}
