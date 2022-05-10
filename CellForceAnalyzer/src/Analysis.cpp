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
    Mat img = cell.m_nucleusChannel.clone();
    help::thresh(img);

    vector<vector<Point> > contours;
    cv::findContours(img, contours, RETR_LIST, CHAIN_APPROX_NONE);
    //if (contours.empty()) { return; }
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

//get largest connected area only??

double analyzePercentageInNucleus(Mat something, Mat nucleus)
{

    something = something.clone();
    nucleus = nucleus.clone();
    help::thresh(something);
    help::thresh(nucleus);

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

vector<int> Analysis::removeBadCells(vector<Cell>& cells)
{
    vector<int> removedCells;
    // if most of the actin lies within the nucleus => cell is dead
    for (int i = 0; i < cells.size(); i++)
    {
        double actinInNucleusPercentage = analyzePercentageInNucleus(cells[i].m_actinChannel, cells[i].m_nucleusChannel);
        if (actinInNucleusPercentage > 0.5) //TODO: generell zu klein gegenüber bild? => area zurück geben
        {
            removedCells.push_back(i);
            cells.erase(cells.begin() + i);
        }
    }
    return removedCells;
}

vector<double> analyseAreaAndDensity(Mat channel)
{
    Mat channelThresh = channel.clone();
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

void drawAxis(Mat& img, Point p, Point q, Scalar colour, const float scale = 0.2)
{
    double angle = atan2((double)p.y - q.y, (double)p.x - q.x); // angle in radians
    double hypotenuse = sqrt((double)(p.y - q.y) * (p.y - q.y) + (p.x - q.x) * (p.x - q.x));
    // Here we lengthen the arrow by a factor of scale
    q.x = (int)(p.x - scale * hypotenuse * cos(angle));
    q.y = (int)(p.y - scale * hypotenuse * sin(angle));
    line(img, p, q, colour, 1);
    // create the arrow hooks
    p.x = (int)(q.x + 9 * cos(angle + CV_PI / 4));
    p.y = (int)(q.y + 9 * sin(angle + CV_PI / 4));
    line(img, p, q, colour, 1);
    p.x = (int)(q.x + 9 * cos(angle - CV_PI / 4));
    p.y = (int)(q.y + 9 * sin(angle - CV_PI / 4));
    line(img, p, q, colour, 1);
}

//TODO: better not with contours but with all the points?
double getPCAorientation(const vector<Point>& pts, vector<Point>& resVec)
{
    // the pca analysis uses another input format for the contour data
    Mat data_pts = Mat(pts.size(), 2, CV_64F);
    for (int i = 0; i < data_pts.rows; i++)
    {
        data_pts.at<double>(i, 0) = pts[i].x;
        data_pts.at<double>(i, 1) = pts[i].y;
    }

    //Perform PCA analysis
    PCA pca_analysis(data_pts, Mat(), PCA::DATA_AS_ROW);
    //Store the center of the object
    Point cntr = Point(static_cast<int>(pca_analysis.mean.at<double>(0, 0)), static_cast<int>(pca_analysis.mean.at<double>(0, 1)));

    //Store the eigenvalues and eigenvectors
    vector<Point2d> eigen_vecs(2);
    vector<double> eigen_val(2);
    for (int i = 0; i < 2; i++)
    {
        eigen_vecs[i] = Point2d(pca_analysis.eigenvectors.at<double>(i, 0), pca_analysis.eigenvectors.at<double>(i, 1));
        eigen_val[i] = pca_analysis.eigenvalues.at<double>(i);
    }

    // to Draw the principal components
    Point p1 = cntr + 0.02 * Point(static_cast<int>(eigen_vecs[0].x * eigen_val[0]), static_cast<int>(eigen_vecs[0].y * eigen_val[0]));
    Point p2 = cntr - 0.02 * Point(static_cast<int>(eigen_vecs[1].x * eigen_val[1]), static_cast<int>(eigen_vecs[1].y * eigen_val[1]));

    double angle = atan2(eigen_vecs[0].y, eigen_vecs[0].x); // orientation in radians

    angle = angle * (360 / (2 * CV_PI));

    resVec[0] = cntr;
    resVec[1] = p1;
    return angle;
}


void Analysis::analyseActin(Cell& cell)
{
    Mat actin = cell.m_actinChannel;
    Mat actinThresh = actin.clone();
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

    vector<Point> resPCAaxis(2);
    double angle = getPCAorientation(longestContour, resPCAaxis);
    cell.actin_PCAangle = angle;

}

double median(vector<double> v)
{
    int n = v.size() / 2;
    nth_element(v.begin(), v.begin() + n, v.end());
    return v[n];
}

double average(vector<double> v)
{
    return accumulate(v.begin(), v.end(), 0) / v.size();
}

//better name?
double Analysis::findDataPointWithMostNeighbours(vector<double> input, double margin)
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

    return median(results);
}

// function used for depiction in show image function, draws centroid, inner & outer radii,...
bool Analysis::analyseShape(Mat& img)
{
    bool successful = help::thresh(img);
    if (!successful)
    {
        return false;
    }

    vector<vector<Point> > contours;
    cv::findContours(img, contours, RETR_LIST, CHAIN_APPROX_NONE);
    if (contours.empty()) { return false; }
    Point centroidLargest;
    vector<Point> longestContour;
    getContourAndCentroidOfLargestWhiteRegion(img, contours, longestContour, centroidLargest);

    vector<Point> radiusPoints(2);
    double innerRadius, outerRadius;
    getRadii(longestContour, centroidLargest, radiusPoints, innerRadius, outerRadius);

    vector<Point> furthestPoints;
    findFurthestDistance(longestContour, furthestPoints);

    vector<Point> resPCAaxis(2);
    getPCAorientation(longestContour, resPCAaxis);

    cv::cvtColor(img, img, COLOR_GRAY2RGB);

    cv::drawContours(img, contours, -1, Scalar(0, 255, 0));
    cv::drawMarker(img, centroidLargest, Scalar(255, 0, 0), 0, 5);
    cv::drawMarker(img, radiusPoints[0], Scalar(0, 0, 255), 0, 5);
    cv::drawMarker(img, radiusPoints[1], Scalar(0, 0, 255), 0, 5);
    cv::line(img, furthestPoints[0], furthestPoints[1], Scalar(255, 255, 0));
    drawAxis(img, resPCAaxis[0], resPCAaxis[1], Scalar(0, 128, 255), 1);

    return true;
}