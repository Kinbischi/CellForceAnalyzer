#include "Display.h"

#include <QMessagebox>
#include <QTablewidget>

using namespace std;
using namespace cv;

Display::Display(ParametersFromUI& p, std::vector<CustomImage>& arr, std::vector<cv::Mat>& y, std::vector<Cell>& c, std::vector<Cell>& d, Analysis& anal, AnalysisFiberDirection& analFiber) :
    m_params(p), m_arrayImages(arr), m_arrayImages_withYoloBoxes(y), m_cellImages(c), m_deletedCellImages(d), m_analysis(anal), m_analysisFiberDir(analFiber) {};


int Display::getImageToShow(Mat& outImg, string& name, double& scale)
{
    if (m_arrayImages.empty()) { return 1; }

    int imageNumber = m_params.showNumber;
    channelType channel = m_params.channel;
    CustomImage image;

    if (m_params.cellArrays)
    {
        if (scale == -1) { scale = 0.4; }
        image = m_arrayImages[imageNumber];
        outImg = image.getChannel(channel);
        name = image.getName(channel);
    }
    if (m_params.singleCells)
    {
        if (m_cellImages.empty()) { return false; }
        if (scale == -1) { scale = 4; }
        image = m_cellImages[imageNumber];
        outImg = image.getChannel(channel);
        name = image.getName(channel);
    }
    if (m_params.deletedCells)
    {
        if (m_deletedCellImages.empty()) { return false; }
        if (scale == -1) { scale = 4; }
        image = m_deletedCellImages[imageNumber];
        outImg = image.getChannel(channel);
        name = image.getName(channel);
    }
    if (m_params.cellArraysWithBoxes)
    {
        if (scale == -1) { scale = 0.4; }
        image = m_arrayImages[imageNumber]; // only for name
        outImg = m_arrayImages_withYoloBoxes[imageNumber];
        name = image.getName();
    }

    if (outImg.empty()) { return 2; }

    outImg = outImg.clone();

    return 0;
}

Mat Display::thresholdImage(Mat image, string& name)
{
    switch (m_params.threshType)
    {
    case thresholdingType::manual:

        help::thresh(image, m_params.manualThreshold);
        name = name + "_threshManually" + to_string(m_params.manualThreshold);
        break;

    case thresholdingType::otsu:
        help::thresh(image);
        name = name + "_threshOtsu";
        break;

    case thresholdingType::squarePCAoptimizedThresh:
        m_analysisFiberDir.getPCAoptThresholdedImage(image);
        name = name + "_threshPCAopt" + "_Length" + to_string(m_params.PCAsquareLength) + "_minEigRatio" + to_string(m_params.PCAminEigValRatio);
        break;
    }
    return image;
}

int Display::prepareAnalysisToShow()
{
    double scaleFactor = m_params.scaleFactor; 
    bool replacingMode = m_params.replacingMode;
    Mat image;
    string name;

    //load image
    int successful = getImageToShow(image, name, scaleFactor);
    if (!(successful==0)) { return successful; }

    if (image.channels() == 1)//this can only be executed for 1 channel images
    {
        Mat thresholdedImage = thresholdImage(image, name);

        if ((m_params.withAnalysis && m_params.dispType == displayType::thresholded) || (!m_params.withAnalysis && m_params.threshType != thresholdingType::None))
        {
            image = thresholdedImage;
        }


        //apply analysis
        bool analysisSuccessful = false;
        if (m_params.variousAnalysis)
        {
            analysisSuccessful = m_analysis.analyseShape(image);
            name = name + "_roundAnalysed";
        }
        else if (m_params.fiberPCA)
        {
            vector<double> random;
            if (m_params.threshType == thresholdingType::None) //intensity mode
            {
                analysisSuccessful = m_analysisFiberDir.analyseWithPCA(image, random);
            }
            else
            {
                analysisSuccessful = m_analysisFiberDir.analyseWithPCA(image, random, thresholdedImage);
            }

            name = name + "_pcaAnalysed";
        }
        else if (m_params.edgeDetection)
        {
            m_analysis.edgeDetectionCanny(image);
            name = name + "_edge" + "_highthresh" + to_string(m_params.highThreshCanny);
        }
        else if (m_params.FAdetection)
        {
            m_analysis.focalAdhesiondetection(image);
            name = name + "_circles" + "_highthresh" + to_string(m_params.highThreshCanny) + "_minCircleConf" +
                to_string(m_params.FAminCircleConfidence) + "_minDist" + to_string(m_params.FAminDist) + "_dp" + to_string(m_params.FAdp);
        }
    }

    if (!m_params.withAnalysis)//TODO: check what happens for unscaled 16 bit image that were analyzed => probs not scaled ;(
    {
        help::scaleData(image);
    }

    string windowName = m_params.replacingMode ? "Image" : name;

    help::showWindow(image, m_params.scaleFactor, windowName);
}
