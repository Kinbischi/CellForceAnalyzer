#include "Analysis.h"

#include <vector>
#include <numeric>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "CustomImage.h"
#include "helperFunctions.h"
#include <cmath>

using namespace std;
using namespace cv;

Analysis::Analysis(ParametersFromUI& parameters) :m_parameters(parameters) {}

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

double analyzePercentageInNucleus(Mat something, Mat nucleus)
{
    Mat somethingThresh = something;
    Mat nucleusThresh = nucleus;
    help::thresh(somethingThresh);
    help::thresh(nucleusThresh);

    int somethingInNucleusCount = 0;
    int somethingOutsideNucleusCount = 0;
    int somethingArea = 0;
    double somethingInNucleusIntensity = 0;
    double somethingOutsideNucleusIntensity = 0;
    for (int i = 0; i < somethingThresh.rows; i++)
    {
        for (int j = 0; j < somethingThresh.cols; j++)
        {
            uchar nThresh = nucleusThresh.at<uchar>(i, j);
            uchar sThresh = somethingThresh.at<uchar>(i, j);

            int sValue;
            if (something.depth() == CV_8U)
            {
                sValue = something.at<uchar>(i, j);
            }
            if (something.depth() == CV_16U)
            {
                sValue = something.at<ushort>(i, j);
            }

            if (sThresh == CustomImage::max8bit)
            {
                somethingArea++;
                if (nThresh == CustomImage::max8bit)
                {
                    somethingInNucleusCount++;
                    somethingInNucleusIntensity += sValue;
                }
                if (nThresh == 0)
                {
                    somethingOutsideNucleusCount++;
                    somethingOutsideNucleusIntensity += sValue;
                }
            }
        }
    }

    //intensity based calculation
    double percentageIntensityInNucleus = somethingInNucleusIntensity / (somethingInNucleusIntensity + somethingOutsideNucleusIntensity);
    double percentageIntensityOutsideNucleus = somethingOutsideNucleusIntensity / (somethingInNucleusIntensity + somethingOutsideNucleusIntensity);

    //pixel wise calculation
    double percentageOfSomethingInNucleus = double(somethingInNucleusCount) / somethingArea;
    double percentageOfSomethingOutsideNucleus = double(somethingOutsideNucleusCount) / somethingArea;

    return percentageIntensityInNucleus;
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


bool Analysis::pointCloudPCA(const vector<Point>& pts, vector<Point2d>& eigen_vecs, shared_ptr<double> eigValRatio )
{
    int squareLength = m_parameters.PCAsquareLength;
    double minEigVecRatio = m_parameters.PCAminEigValRatio;

    if (pts.size() < help::minPercentagePointsPCA * squareLength*squareLength || pts.size() < help::minPointsPCA) // enough white pixels
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

if (eigValRatio != nullptr)
{
    *eigValRatio = ratio;
}

if (ratio < minEigVecRatio)
{
    return false;
}

return true;
}

vector<Point> getWhitePointsFromThresholdedImage(Mat img)
{
    vector<Point> pointCloud;
    for (int i = 0; i < img.rows; i++)
    {
        for (int j = 0; j < img.cols; j++)
        {
            if (img.at<uchar>(i, j) == 255)
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
    if (angle < (help::startAnglePlot - help::spacingAnglePlot / 2))// angle goes from 15 deg to 180 deg
    {
        angle += 180;
    }
    return angle;
}

vector<Point> getPointsDependingOnIntensityFromImage(Mat img)
{
    vector<Point> pointCloud;
    for (int i = 0; i < img.rows; i++)
    {
        for (int j = 0; j < img.cols; j++)
        {
            int pixelIntensity = img.at<uchar>(i, j);
            for (int k = 0; k < pixelIntensity; k++)
            {
                Point p = Point(j, i);
                pointCloud.push_back(p);
            }
        }
    }
    return pointCloud;
}

bool squareIsNotBackground(Mat img)
{
    int area = img.rows * img.cols;
    double summedPixelIntensity = 0;
    for (int i = 0; i < img.rows; i++)
    {
        for (int j = 0; j < img.cols; j++)
        {
            int pixelIntensity = img.at<uchar>(i, j);
            summedPixelIntensity += pixelIntensity;
        }
    }

    if (summedPixelIntensity / area < help::minAverageIntensityForNotBackground) { return false; }
    else { return true; }
}

void Analysis::edgeDetectionCanny(Mat& image)
{
    cv::blur(image, image, Size(3, 3));
    int kernel_size = 3;
    cv::Canny(image, image, m_parameters.highThreshCanny / 2, m_parameters.highThreshCanny, kernel_size);
}

void Analysis::focalAdhesiondetection(Mat& image)
{
    Mat blurredImage;
    GaussianBlur(image, blurredImage, Size(9, 9), 2, 2);

    vector<Vec3f> circles;
    HoughCircles(blurredImage, circles, HOUGH_GRADIENT, m_parameters.FAdp, m_parameters.FAminDist, 
        m_parameters.highThreshCanny, m_parameters.FAminCircleConfidence, 0, 10);

    cvtColor(image, image, COLOR_GRAY2RGB);

    for (size_t i = 0; i < circles.size(); i++)
    {
        Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
        int radius = cvRound(circles[i][2]);

        circle(image, center, radius, Scalar(0, 0, 255), 1);
    }
}

void Analysis::getPCAoptThresholdedImage(Mat& img)
{
    int squareLength = m_parameters.PCAsquareLength;
    int lengthX = squareLength; //length of mini-squares that are fed into the pca
    int lengthY = squareLength;

    int edgeX = (img.cols % lengthX) / 2;
    int edgeY = (img.rows % lengthY) / 2;

    img = img.clone();
    help::scaleData(img);
    if (img.depth() == CV_16U) { img.convertTo(img, CV_8U, 1 / 256.0); }

    for (int i = 0; i < img.cols / lengthX; i++)
    {
        for (int j = 0; j < img.rows / lengthY; j++)
        {
            Rect rect(edgeX + i * lengthX, edgeY + j * lengthY, lengthX, lengthY);
            Mat subImg = img(rect);
            bool notBackground = squareIsNotBackground(subImg);

            double optEigValRatio = 0;
            int optimalThresholding;
            vector<Point2d> eigen_vecs(2);
            for (int k = 10; k < 256; k += 5) //heavy function => thresholds a lot of times!!
            {
                Mat threshImg = subImg.clone();

                help::thresh(threshImg, k, false);
                vector<Point> points = getWhitePointsFromThresholdedImage(threshImg);
                
                bool enoughPoints = true; // only one of the conditions from pointCloudPCA function is needed: at least X percent of area in white
                if (points.size() < help::minPercentagePointsPCA * squareLength * squareLength || points.size() < help::minPointsPCA) 
                {
                    enoughPoints = false;
                }

                std::shared_ptr<double> eigValRatio = std::make_shared<double>(); // shared pointer => optional arguments can not be a reference ;(
                pointCloudPCA(points, eigen_vecs, eigValRatio);

                if (*eigValRatio > optEigValRatio && enoughPoints && notBackground)
                {
                    optEigValRatio = *eigValRatio;
                    optimalThresholding = k;
                }
            }
            //TODO was machsch du da isch optEigVal überhaut existent??
            if (optEigValRatio < m_parameters.PCAminEigValRatio)
            {
                if (m_parameters.suppressLowEigValRatioSquares)
                {
                    notBackground = false;
                }
            }

            if (notBackground == true)
            {
                cv::threshold(subImg, subImg, optimalThresholding, 255, THRESH_BINARY);
            }
            else
            {
                cv::threshold(subImg, subImg, 255, 255, THRESH_BINARY);
            }

        }
    }
}


bool Analysis::analyseWithPCA(Mat& img, vector<double>& resultingAngles, Mat imgThresh)
{
    int squareLength = m_parameters.PCAsquareLength;
    int lengthX = squareLength; //length of mini-squares that are fed into the pca
    int lengthY = squareLength;
    
    int edgeX = (img.cols % lengthX) / 2;
    int edgeY = (img.rows % lengthY) / 2;

    
    help::scaleData(img);
    if (img.depth() == CV_16U) { img.convertTo(img, CV_8U, 1 / 256.0); }
    Mat imgCopy = img.clone(); // copy of 8-bit scaled image
    cv::cvtColor(img, img, COLOR_GRAY2RGB);

    for (int i = 0; i < img.cols/lengthX; i++)
    {
        for (int j = 0; j < img.rows/lengthY; j++)
        {
            Rect rect(edgeX + i*lengthX, edgeY + j*lengthY, lengthX, lengthY);

            Mat subImg;
            vector<Point> points;
            bool enoughPoints = true;
            if (imgThresh.empty()) //without thresholding => intensity mode
            {
                subImg = imgCopy(rect);
                points = getPointsDependingOnIntensityFromImage(subImg);
                if (points.size() < help::minAverageIntensityForNotBackground * squareLength * squareLength)
                {
                    enoughPoints = false;
                }
            }
            else
            {
                subImg = imgThresh(rect);
                points = getWhitePointsFromThresholdedImage(subImg);
            }
            
            vector<Point2d> eigen_vecs(2);
            bool PCAworked = pointCloudPCA(points, eigen_vecs); // TODO: all conditions?

            if (PCAworked && enoughPoints)
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

double Analysis::calculateFiberAlignmentConstant(vector<int> y, double exponent, int mode)
{
    double result = -1;
    double yExpSum = 0;
    double ySum = 0;
    for (int i = 0; i < y.size(); i++)
    {
        ySum += y[i];
        yExpSum += pow(y[i], exponent);
    }

    if (mode == 0)
    {
        result = yExpSum / pow(ySum, exponent);
    }
    if (mode == 1)
    {
        result = 1 - pow(ySum, exponent)/(yExpSum*pow(y.size(),exponent-1));
    }
    if (mode == 2)
    {
        double temp = pow(y.size(), exponent - 1);
        result = yExpSum / pow(ySum, exponent) * temp / (temp - 1) - 1 / (temp - 1);
    }
    
    return result;
    //double yAverage = ySum / y.size();
    //return (y.size() * pow(yAverage, exponent)) / expSum;
}

vector<double> Analysis::testYvecs(vector<vector<int>> yVecs, double exponent, int mode)
{
    vector<double> result;
    for (int i = 0; i < yVecs.size(); i++)
    {
        auto y = yVecs[i];
        double res = calculateFiberAlignmentConstant(y, exponent, mode);
        result.push_back(res); 
    }
    return result;
}


vector<double> Analysis::createX(vector<double> data1, bool plottingAngles)
{
    double start, end, spacing;
    vector<double> x;
    if (plottingAngles)
    {
        start = help::startAnglePlot;
        end = 180 + start;
        spacing = help::spacingAnglePlot;
    }
    else
    {
        int pltDetail = 12;  
        double min = *std::min_element(data1.begin(), data1.end());
        double max = *std::max_element(data1.begin(), data1.end());
        spacing = (max - min) / pltDetail;
        start = min - spacing;
        end = max + spacing;
    }

    for (double i = start; i < end; i += spacing)
    {
        x.push_back(i);
    }
    return x;
}

vector<int> Analysis::createY(vector<double> data, vector<double> x)
{
    std::vector<int> y(x.size());
    double spacing = x[1] - x[0];

    for (int i = 0; i < data.size(); i++)
    {
        for (int j = 0; j < x.size(); j++)
        {
            if (data[i] >= x[j] - spacing / 2 && data[i] < x[j] + spacing / 2)
            {
                y[j]++;
            }
        }
    }
    return y;
}

// falls altered => cell info müsste als resultVector zurück gebeben
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

    bool worked = pointCloudPCA(points, eigen_vecs); //TOCHANGE
    double mainAngle = getAngleFromVectors(eigen_vecs);
    cell.actin_mainAngle = mainAngle;

    vector<double> fiberAngles;
    Mat optThresholdedActin = actin.clone();
    getPCAoptThresholdedImage(optThresholdedActin);
    analyseWithPCA(actin, fiberAngles, optThresholdedActin); //changes pic!!
    cell.actin_fibreAnglesPCA = fiberAngles;

    double spacing;
    vector<double> x = createX(fiberAngles,true);
    vector<int> yFiberAngles = createY(fiberAngles,x);
    cell.actin_fiberAlignment = calculateFiberAlignmentConstant(yFiberAngles,3,1);
}

// function used for depiction in show image function, draws centroid, inner & outer radii,...
// not very beautiful but does the job
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
    bool worked = pointCloudPCA(points, eigen_vecs); //TOCHANGE
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


// TODO: fitEllipse() => approximates ellipse into shape and returns rectangle => could be better than inner/outer radii
// 
// function used to get circularity and roundness for every cell
void analyseNucleus(Cell& cell)
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

void analyseYap(Cell& cell)
{
    double yapInNucleus = analyzePercentageInNucleus(cell.m_yapChannel, cell.m_nucleusChannel);
    cell.yap_inNucleus = yapInNucleus;
}


// TODO: wenn fixed (e.g.x20) => regulate nucleus size => between x and y pixels
// actin detatched from nucleus?
bool Analysis::isDeadCell(Cell cell)
{
    // 1)
    // triggers if most of the actin lies within the nucleus
    // or nucleus is huge => perhaps 2 molten nucli together
    double actinInNucleusPercentage = analyzePercentageInNucleus(cell.m_actinChannel, cell.m_nucleusChannel);
    if (actinInNucleusPercentage > 0.4)
    {
        return true;
    }
    // actin area should be bigger than nucleus area
    vector<double> result;
    result = analyseAreaAndDensity(cell.m_actinChannel);
    double actin_area = result[0];
    result = analyseAreaAndDensity(cell.m_nucleusChannel);
    double nucleus_area = result[0];
    auto hi = actin_area / nucleus_area;
    if (actin_area / nucleus_area < 2)
    {
        return true;
    }
    

    // 2)
    // if thresholding of nucleus has multiple areas => multiple cells in niche or very very blurred nucleus
    // spots are with small area are not considered
    Mat nucleusThresh = cell.m_nucleusChannel;
    help::thresh(nucleusThresh);
    Mat labels, stats, centroids;
    int nucleusParts = cv::connectedComponentsWithStats(nucleusThresh, labels, stats, centroids);

    if (nucleusParts == 1) // no nucleus, only background
    {
        return true;
    }
    if (nucleusParts!=2) // more than background and nucleus
    {
        int bigNucleusParts = 0;
        for (int i = 1; i < nucleusParts; i++) // 0 not used => background label
        {
            int area = stats.at<int>(i, CC_STAT_AREA);
            if (area > 10)
            {
                bigNucleusParts++;
            }
        }
        if (bigNucleusParts>1)
        {
            return true;
        }
    }

    return false;
}

void Analysis::analyseCell(Cell& cell, std::vector<channelType> channels)
{
    bool hasNucleus, hasActin, hasYap;
    hasNucleus = std::count(channels.begin(), channels.end(), channelType::nucleus);
    hasActin = std::count(channels.begin(), channels.end(), channelType::actin);
    hasYap = std::count(channels.begin(), channels.end(), channelType::yap);
    if (hasNucleus)
    {
        analyseNucleus(cell);
    }
    if (hasActin)
    {
        analyseActin(cell);
    }
    if (hasYap&&hasNucleus)
    {
        analyseYap(cell);
    }
    
    
}