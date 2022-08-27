#include "Display.h"

#include <QMessagebox>
#include <QTablewidget>

using namespace std;
using namespace cv;

Display::Display(ParametersUI& p, dataContainer& d, Analysis& anal, AnalysisFiberDirection& analFiber) :
    params(p), data(d),analysis(anal), analysisFiberDir(analFiber) {};


int Display::getImageToShow(Mat& outImg, string& name, double& scale)
{
    if (data.arrayImages.empty() && data.cellImages.empty() && data.removedCellImages.empty()) { return 1; }

    int imageNumber = params.showNumber;
    channelType channel = params.channel;
    CustomImage image;

    if (params.cellArrays)
    {
        if (scale == -1) { scale = 0.4; }
        image = data.arrayImages[imageNumber];
        outImg = image.getChannel(channel);
        name = image.getName(channel);
    }
    if (params.singleCells)
    {
        if (data.cellImages.empty()) { return 3; }
        if (scale == -1) { scale = 4; }
        image = data.cellImages[imageNumber];
        outImg = image.getChannel(channel);
        name = image.getName(channel);
    }
    if (params.removedCells)
    {
        if (data.removedCellImages.empty()) { return 3; }
        if (scale == -1) { scale = 4; }
        image = data.removedCellImages[imageNumber];
        outImg = image.getChannel(channel);
        name = image.getName(channel);
    }
    if (params.cellArraysWithBoxes)
    {
        if (scale == -1) { scale = 0.4; }
        image = data.arrayImages[imageNumber]; // only for name
        outImg = data.arrayImages_withYoloBoxes[imageNumber];
        name = image.getName();
    }

    if (outImg.empty()) { return 2; }

    outImg = outImg.clone();

    return 0;
}

void Display::blurImg(Mat& img, string& name)
{
    if (params.blurWeak)
    {
        analysis.blurWeak(img);
        name = name + "_blurredWeak";
    }
    else if (params.blurStrong)
    {
        analysis.blurStrong(img);
        name = name + "_blurredStrong";
    }
}

void fillHoles(Mat& img)
{
    vector<vector<Point> > contours;
    cv::findContours(img, contours, RETR_CCOMP, CHAIN_APPROX_SIMPLE);

    for (int i = 0; i < contours.size(); i++)
    {
        cv::drawContours(img, contours, i, 255 ,FILLED);
    }
}

Mat Display::thresholdImage(Mat image, string& name)
{
    switch (params.threshType)
    {
    case thresholdingType::manual:

        Analysis::thresh(image, params.manualThreshold,0);
        name = name + "_threshManually" + to_string(params.manualThreshold);
        break;

    case thresholdingType::otsu:
        Analysis::thresh(image,0,0);
        name = name + "_threshOtsu";
        if (params.suppressLowEigValRatioSquares)
        {
            fillHoles(image);
        }
        
        break;

    case thresholdingType::squarePCAoptimizedThresh:
        analysisFiberDir.getPCAoptThresholdedImage(image);
        name = name + "_threshPCAopt" + "_Length" + to_string(params.PCAsquareLength) + "_minEigRatio" + to_string(params.PCAminEigValRatio);
        break;
    }
    return image;
}

int Display::prepareAnalysisToShow()
{
    double scaleFactor = params.scaleFactor;
    bool replacingMode = params.replacingMode;
    Mat image;
    string name;

    //load image
    int successful = getImageToShow(image, name, scaleFactor);
    if (!(successful==0)) { return successful; }

    if (image.channels() == 1)//this can only be executed for 1 channel images
    {
        blurImg(image, name);
        Mat thresholdedImage = thresholdImage(image, name);

        if ((params.withAnalysis && params.dispType == displayType::thresholded) || (!params.withAnalysis && params.threshType != thresholdingType::None))
        {
            image = thresholdedImage;
        }

        //apply analysis
        bool analysisSuccessful = false;
        if (params.variousAnalysis)
        {
            analysisSuccessful = analysis.analyseShape(image);
            name = name + "_roundAnalysed";
        }
        else if (params.fiberPCA)
        {
            vector<double> random;
            if (params.threshType == thresholdingType::None) //intensity mode
            {
                analysisSuccessful = analysisFiberDir.analyseWithPCA(image, random);
            }
            else
            {
                analysisSuccessful = analysisFiberDir.analyseWithPCA(image, random, thresholdedImage);
            }

            name = name + "_pcaAnalysed";
        }
        else if (params.edgeDetection)
        {
            analysis.edgeDetectionCanny(image);
            name = name + "_edge" + "_highthresh" + to_string(params.highThreshCanny);
        }
        else if (params.FAdetection)
        {
            analysis.focalAdhesiondetection(image);
            name = name + "_circles" + "_highthresh" + to_string(params.highThreshCanny) + "_minCircleConf" +
                to_string(params.FAminCircleConfidence) + "_minDist" + to_string(params.FAminDist) + "_dp" + to_string(params.FAdp);
        }
    }

    if (!params.withAnalysis)//TODO: check what happens for unscaled 16 bit image that were analyzed => probs not scaled ;(
    {
        help::scaleData(image);
    }

    string windowName = params.replacingMode ? "Image" : name;

    help::showWindow(image, scaleFactor, windowName);

    return 0;
}
