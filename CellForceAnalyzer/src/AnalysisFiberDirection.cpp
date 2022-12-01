#include "AnalysisFiberDirection.h"
#include "Analysis.h"
using namespace std;
using namespace cv;


AnalysisFiberDirection::AnalysisFiberDirection(ParametersUI& p) :params(p) {};


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

bool AnalysisFiberDirection::pointCloudPCA(const vector<Point>& pts, vector<Point2d>& eigen_vecs, shared_ptr<double> eigValRatio)
{
    int squareLength = params.PCAsquareLength;
    double minEigVecRatio = params.PCAminEigValRatio;

    if (pts.size() < help::minPercentagePointsPCA * squareLength * squareLength || pts.size() < help::minPointsPCA) // enough white pixels
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

void AnalysisFiberDirection::getPCAoptThresholdedImage(Mat& img)
{
    int squareLength = params.PCAsquareLength;
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

                Analysis::thresh(threshImg, k, 0,false);
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
            if (optEigValRatio < params.PCAminEigValRatio)
            {
                if (params.suppressLowEigValRatioSquares)
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


bool AnalysisFiberDirection::analyseWithPCA(Mat& img, vector<double>& resultingAngles, Mat imgThresh)
{
    int squareLength = params.PCAsquareLength;
    int lengthX = squareLength; //length of mini-squares that are fed into the pca
    int lengthY = squareLength;

    int edgeX = (img.cols % lengthX) / 2;
    int edgeY = (img.rows % lengthY) / 2;


    help::scaleData(img);
    if (img.depth() == CV_16U) { img.convertTo(img, CV_8U, 1 / 256.0); }
    Mat imgCopy = img.clone(); // copy of 8-bit scaled image
    cv::cvtColor(img, img, COLOR_GRAY2RGB);

    for (int i = 0; i < img.cols / lengthX; i++)
    {
        for (int j = 0; j < img.rows / lengthY; j++)
        {
            Rect rect(edgeX + i * lengthX, edgeY + j * lengthY, lengthX, lengthY);

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
                
                help::drawDoubleArrow(img, center, eigen_vecs, Scalar(0, 128, 255), lengthX);
                rectangle(img, rect, Scalar(0, 255, 128));
                double angle = getAngleFromVectors(eigen_vecs);
                resultingAngles.push_back(angle);
            }
        }
    }
    if (resultingAngles.size() < 5) { return false; }
    return true;
}

double AnalysisFiberDirection::calculateFiberAlignmentConstant(vector<int> y, double exponent, int mode)
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
        result = 1 - pow(ySum, exponent) / (yExpSum * pow(y.size(), exponent - 1));
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

vector<double> AnalysisFiberDirection::testYvecs(vector<vector<int>> yVecs, double exponent, int mode)
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