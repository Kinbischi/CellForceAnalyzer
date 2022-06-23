#include "Analysis.h"

#include <vector>
#include <numeric>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "CustomImage.h"
#include "helperFunctions.h"

using namespace std;
using namespace cv;



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

// works with 8 bit images
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
    // 
    //double isInContour = cv::pointPolygonTest(longestContour, centroidLargest, false);
    //assert(isInContour + 1);
}

// finds furthest distance between points of a contour
double findFurthestDistance(vector<Point> contour, vector<Point>& furthestPoints)
{
    vector<Point> hullPoints;
    cv::convexHull(contour, hullPoints);

    int start = 1;
    double maxDistance = 0;
    Point firstPoint, secondPoint;
    for (int i = 0; i < hullPoints.size() - 1; i++)
    {
        for (int j = start; j < hullPoints.size(); j++)
        {
            double distance = cv::norm(hullPoints[i] - hullPoints[j]);
            if (distance > maxDistance)
            {
                maxDistance = distance;
                firstPoint = hullPoints[i];
                secondPoint = hullPoints[j];
            }
        }
        start++;
    }
    furthestPoints.push_back(firstPoint);
    furthestPoints.push_back(secondPoint);
    return maxDistance;
}

// function used to get circularity and roundness for every cell

// TODO:
// fitEllipse() => approximates ellipse into shape and returns rectangle => could be better than inner/outer radii
void Analysis::analyseNucleusShape(Cell& cell)
{
    Mat nucleusThresh = cell.m_nucleusChannel;
    help::thresh(nucleusThresh);

    vector<vector<Point> > contours;
    cv::findContours(nucleusThresh, contours, RETR_LIST, CHAIN_APPROX_NONE);
    //if (contours.empty()) { return; }
    Point centroidLargest;
    vector<Point> longestContour;
    getContourAndCentroidOfLargestWhiteRegion(nucleusThresh, contours, longestContour, centroidLargest);

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

//get largest connected area only??
// TODO change yap from only thresholded to also weighted with pixel values?
//TODO: vielleicht eine gesamtzell definition welche vom actin, yap, nucleus das thresholding vereinigt?
double analyzePercentageInNucleus(Mat something, Mat nucleus)
{

    something = something.clone();
    nucleus = nucleus.clone();
    auto d = nucleus.depth();
    auto d2 = something.depth();
    Mat something2 = something.clone();

    //auto m1 = help::elementsMap<ushort>(something);
    //auto m2 = help::elementsMap<ushort>(nucleus);

    help::thresh(something);
    help::thresh(nucleus);
    help::showWindow(something);
    auto d3 = nucleus.depth();
    
    help::scaleData(something2);
    help::thresh(something2);
    //auto m5 = help::elementsMap<uchar>(something2);
    help::showWindow(something2);
    help::showWindow(something);

    //auto m3 = help::elementsMap<uchar>(something);
    //auto m4 = help::elementsMap<uchar>(nucleus);
    

    int somethingInNucleusCount = 0;
    int somethingOutsideNucleusCount = 0;
    int somethingArea = 0;
    for (int i = 0; i < something.rows; i++)
    {
        for (int j = 0; j < something.cols; j++)
        {
            uchar n = nucleus.at<uchar>(i, j);
            uchar s = something.at<uchar>(i, j);

            if (s == CustomImage::max8bit)
            {
                somethingArea++;
                if (n == CustomImage::max8bit)
                {
                    somethingInNucleusCount++;
                }
                if (n == 0)
                {
                    somethingOutsideNucleusCount++;
                }
            }
        }
    }

    double percentageOfSomethingInNucleus = double(somethingInNucleusCount) / somethingArea;
    double percentageOfSomethingOutsideNucleus = double(somethingOutsideNucleusCount) / somethingArea;
    return percentageOfSomethingInNucleus;
}

void Analysis::analyseYapInNucleus(Cell& cell)
{
    double yapInNucleus = analyzePercentageInNucleus(cell.m_yapChannel, cell.m_nucleusChannel);
    cell.yap_inNucleus = yapInNucleus;
}

bool Analysis::isDeadCell(Cell cell)
{
    // if most of the actin lies within the nucleus => cell is dead
    //TODO: generell zu klein gegenüber bild? => area zurück geben

    double actinInNucleusPercentage = analyzePercentageInNucleus(cell.m_actinChannel, cell.m_nucleusChannel);
    if (actinInNucleusPercentage > 0.5) 
    {
        return true;
    }
    else 
    {
        return false;
    }
}

vector<double> analyseAreaAndDensity(Mat channel)
{
    Mat channelThresh = channel;
    help::thresh(channelThresh);
    int pixelSum = 0;
    int numberOfPixels = 0;

    for (int i = 0; i < channel.rows; i++)
    {
        for (int j = 0; j < channel.cols; j++)
        {
            if (channelThresh.at<uchar>(i, j) != 0)
            {
                if (channel.depth() == CV_8U)
                {
                    pixelSum += channel.at<uchar>(i, j);
                }
                if (channel.depth() == CV_16U)
                {
                    pixelSum += channel.at<ushort>(i, j);
                }
                numberOfPixels++;
            }
        }
    }

    double actinDensity = double(pixelSum) / numberOfPixels;
    vector<double> result;
    result.push_back(numberOfPixels);
    result.push_back(actinDensity);
    return result;
}

bool Analysis::pointCloudPCA(const vector<Point>& pts, const int squareLength, vector<Point2d>& eigen_vecs)
{
    if (pts.size() < 3)
    {
        return false;
    }
    // the pca analysis uses another input format for the data
    Mat data_pts = Mat(pts.size(), 2, CV_64F);
    for (int i = 0; i < data_pts.rows; i++)
    {
        data_pts.at<double>(i, 0) = pts[i].x;
        data_pts.at<double>(i, 1) = pts[i].y;
    }

    //Perform PCA analysis
    PCA pca_analysis(data_pts, Mat(), PCA::DATA_AS_ROW);
    // center = Point(static_cast<int>(pca_analysis.mean.at<double>(0, 0)), static_cast<int>(pca_analysis.mean.at<double>(0, 1)));

    //Store the eigenvalues and eigenvectors
    vector<double> eigen_val(2);
    for (int i = 0; i < 2; i++)
    {
        eigen_vecs[i] = Point2d(pca_analysis.eigenvectors.at<double>(i, 0), pca_analysis.eigenvectors.at<double>(i, 1));
        eigen_val[i] = pca_analysis.eigenvalues.at<double>(i);
    }

    // for fiber analysis: analysis only worked if:
    // - at least 10 percent of image need to be covered with fibers (white points)
    // - in the analysis one principal component is clearly more dominant than the other => this can be seen via the eigenvalues
    int subImgArea = squareLength * squareLength;
    double ratio = eigen_val[0] / eigen_val[1];
    if (pts.size() < 0.1 * subImgArea || ratio < 1.5)
    {
        return false;
    }

    return true;
}

vector<Point> Analysis::getWhitePointsFromThresholdedImage(Mat img)
{
    vector<Point> pointCloud;
    for (int i = 0; i < img.rows; i++)
    {
        for (int j = 0; j < img.cols; j++)
        {
            if (img.at<uchar>(i,j)==255)
            {
                Point p = Point(j, i);
                pointCloud.push_back(p);
            }
        }
    }
    return pointCloud;
}

double getAngleFromVectors(vector<Point2d> vectors)
{
    double angle = atan2(vectors[0].y, vectors[0].x); // orientation in radians
    angle = angle * (360 / (2 * CV_PI));              // to degrees
    angle = angle * (-1);                             // invert angle so that it is counterclock wise
    if (angle < -25)                                  // angle goes from -25 deg to 155 deg => so that 0 and 90 deg are depicted well in plot
    {
        angle += 180;
    }
    return angle;
}

bool Analysis::analyseWithPCA(Mat& img, vector<double>& resultingAngles)
{
    int lengthX = 15; //length of mini-squares that are fed into the pca
    int lengthY = 15;
    
    int edgeX = (img.cols % lengthX)/2;
    int edgeY = (img.rows % lengthY)/2;
    
    Mat imgThresh = img;
    bool successful = help::thresh(imgThresh);
    if (!successful) { return false; }

    help::scaleData(img);
    if (img.depth() == CV_16U) { img.convertTo(img, CV_8U, 1 / 256.0); }
    cv::cvtColor(img, img, COLOR_GRAY2RGB);

    for (int i = 0; i < img.cols/lengthX; i++)
    {
        for (int j = 0; j < img.rows/lengthY; j++)
        {
            Rect rect(edgeX + i*lengthX, edgeY + j*lengthY, lengthX, lengthY);
            Mat subImg = imgThresh(rect);

            vector<Point> points = getWhitePointsFromThresholdedImage(subImg);
            
            vector<Point2d> eigen_vecs(2);
            bool worked = pointCloudPCA(points, lengthX, eigen_vecs);

            if (worked)
            {
                Point center = Point(edgeX + (i + 0.5) * lengthX, edgeY + (j + 0.5) * lengthY);
                rectangle(img, rect, Scalar(0, 255, 128));
                help::drawDoubleArrow(img, center, eigen_vecs, Scalar(0, 128, 255), lengthX);

                double angle = getAngleFromVectors(eigen_vecs);
                resultingAngles.push_back(angle);
            }
        }
    }
    if (resultingAngles.size() < 5) { return false; }
    return true;
}


//TODO introduced worked
void Analysis::analyseActin(Cell& cell)
{
    Mat actin = cell.m_actinChannel.clone();
    Mat actinThresh = actin;
    help::thresh(actinThresh);

    vector<double> result = analyseAreaAndDensity(actin);
    cell.actin_area = result[0];
    cell.actin_density = result[1];

    vector<vector<Point> > contours;
    cv::findContours(actinThresh, contours, RETR_LIST, CHAIN_APPROX_NONE);
    if (contours.empty()) { return; }
    Point centroidLargest;
    vector<Point> longestContour;
    getContourAndCentroidOfLargestWhiteRegion(actinThresh, contours, longestContour, centroidLargest);

    vector<Point> furthestPoints;
    double maxLength = findFurthestDistance(longestContour, furthestPoints);
    cell.actin_maxLength = maxLength;

    vector<Point> points = getWhitePointsFromThresholdedImage(actinThresh);
    vector<Point2d> eigen_vecs(2);
    bool worked = pointCloudPCA(points, actin.cols, eigen_vecs);
    double mainAngle = getAngleFromVectors(eigen_vecs);
    cell.actin_mainAngle = mainAngle;

    vector<double> angles;
    analyseWithPCA(actin, angles);
    cell.actin_fibreAnglesPCA = angles;
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

//TODO: get that out of here??
Cell Analysis::getAverageProperties(vector<Cell> cells)
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

// function used for depiction in show image function, draws centroid, inner & outer radii,...
bool Analysis::analyseShape(Mat& img)
{
    Mat imgThresh = img;
    bool successful = help::thresh(imgThresh);
    if (!successful) { return false; }

    vector<vector<Point> > contours;
    cv::findContours(imgThresh, contours, RETR_LIST, CHAIN_APPROX_NONE);
    if (contours.empty()) { return false; }
    Point centroidLargest;
    vector<Point> longestContour;
    getContourAndCentroidOfLargestWhiteRegion(imgThresh, contours, longestContour, centroidLargest);

    vector<Point> radiusPoints(2);
    double innerRadius, outerRadius;
    getRadii(longestContour, centroidLargest, radiusPoints, innerRadius, outerRadius);

    vector<Point> furthestPoints;
    findFurthestDistance(longestContour, furthestPoints);

    vector<Point> points = getWhitePointsFromThresholdedImage(imgThresh);
    vector<Point2d> eigen_vecs(2);
    bool worked = pointCloudPCA(points, img.cols, eigen_vecs);
    Point ctr = Point(img.cols / 2, img.rows / 2);

    help::scaleData(img);
    if (img.depth() == CV_16U) { img.convertTo(img, CV_8U, 1 / 256.0); }
    cv::cvtColor(img, img, COLOR_GRAY2RGB);

    help::drawDoubleArrow(img, ctr, eigen_vecs, Scalar(0, 128, 255), img.cols);
    cv::drawContours(img, contours, -1, Scalar(0, 255, 0));
    cv::drawMarker(img, centroidLargest, Scalar(255, 0, 0), 0, 5);
    cv::drawMarker(img, radiusPoints[0], Scalar(0, 0, 255), 0, 5);
    cv::drawMarker(img, radiusPoints[1], Scalar(0, 0, 255), 0, 5);
    cv::line(img, furthestPoints[0], furthestPoints[1], Scalar(255, 255, 0));

    return true;
}