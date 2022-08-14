#include "Plotting.h"

using namespace std;
using namespace cv;


#ifndef _DEBUG

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

#endif


Plotting::Plotting(ParametersFromUI& p, AnalysisFiberDirection& a, std::vector<Cell>& c) :params(p), analysisFiberDir(a), cellImages(c) {};

string Plotting::getPlotTitle(plotFeatureType CFType)
{
    string name;
    switch (CFType)
    {
    case plotFeatureType::actFibersOptThresh:
        name = "Actin fiber analysis (optimally thresholded)";
        break;
    case plotFeatureType::actFibersIntensity:
        name = "Actin fiber analysis (intensity-based)";
        break;
    case plotFeatureType::actFibersBoth:
        name = "Actin fiber analysis";
        break;
    case plotFeatureType::nuclArea:
        name = "Nucleus area";
        break;
    case plotFeatureType::nuclCircularity:
        name = "Nucleus circularity";
        break;
    case plotFeatureType::nuclRoundness:
        name = "Nucleus roundness";
        break;
    case plotFeatureType::yapInNucleus:
        name = "YAP in nucleus";
        break;
    case plotFeatureType::actArea:
        name = "Actin area";
        break;
    case plotFeatureType::actDensity:
        name = "Actin density";
        break;
    case plotFeatureType::actMaxLength:
        name = "Actin maximum stretch";
        break;
    case plotFeatureType::actMainAngle:
        name = "Actin main angle";
        break;
    }
    return name;
}

vector<string> Plotting::getPlotLegendNames(plotFeatureType CFType)
{
    vector<string> names(2);
    switch (CFType)
    {
    case plotFeatureType::actFibersBoth:
        names[0] = "Optimally thresholded image PCA";
        names[1] = "Intensity-based PCA";
        break;

    }
    return names;
}

string Plotting::getPlotXLabel(plotFeatureType CFType)
{
    string name;
    switch (CFType)
    {
    case plotFeatureType::actFibersOptThresh:
        name = "Angle (in degrees)";
        break;
    case plotFeatureType::actFibersIntensity:
        name = "Angle (in degrees)";
        break;
    case plotFeatureType::actFibersBoth:
        name = "Angle (in degrees)";
        break;
    default:
        name = getPlotTitle(CFType);
        break;
    }
    return name;
}

string Plotting::getPlotYLabel(plotFeatureType CFType)
{
    string name;
    switch (CFType)
    {
    case plotFeatureType::actFibersOptThresh:
        name = "Number of squares";
        break;
    case plotFeatureType::actFibersIntensity:
        name = "Number of squares";
        break;
    case plotFeatureType::actFibersBoth:
        name = "Number of squares";
        break;
    default:
        name = "Number of occurrences";
        break;
    }
    return name;
}

vector<double> Plotting::createX(vector<double> data1, bool plottingAngles)
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

vector<int> Plotting::createY(vector<double> data, vector<double> x)
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



//TODO:
// set standard square size
void Plotting::plotData(vector<double> data1, bool plottingAngles, vector<double> data2)
{
    vector<double> x = createX(data1, plottingAngles);
    vector<int> y1 = createY(data1, x);

#ifndef _DEBUG
    plt::figure();

    int ylimMax = *max_element(y1.begin(), y1.end());

    vector<string> labels = getPlotLegendNames(params.plotFeatType);

    plt::plot(x, y1, { {"marker", "o"}, {"linestyle", "--"}, {"label", labels[0]} });
    if (!data2.empty())
    {
        vector<int> y2 = createY(data2, x);
        int ylimMax2 = *max_element(y2.begin(), y2.end());
        if (ylimMax2 > ylimMax) { ylimMax = ylimMax2; }
        plt::plot(x, y2, { {"marker", "o"}, {"linestyle", "--"}, {"label", labels[1]} });
        plt::legend();
    }
    ylimMax = ylimMax + std::ceil(0.08 * ylimMax);
    plt::ylim(0, ylimMax);

    plt::xlabel(getPlotXLabel(params.plotFeatType));
    plt::ylabel(getPlotYLabel(params.plotFeatType));
    plt::title(getPlotTitle(params.plotFeatType));
    plt::show();
#endif
}


void Plotting::plotSomething(vector<double> x, vector<double> y, plotFeatureType pltType, string title)
{
#ifndef _DEBUG
    plt::figure();

    plt::plot(x, y, { {"marker", "o"}, {"linestyle", "--"} });

    plt::xlabel(getPlotXLabel(pltType));
    plt::ylabel(getPlotYLabel(pltType));
    plt::ylim(0, 10);
    plt::title(title);
    //plt::savefig(title +".pdf");
    plt::show();
#endif
}

/*
int Plotting::getImageToPlot(Mat& outImg)
{
    if (arrayImages.empty()) { return 1; }

    int imageNumber = params.showNumber;
    channelType channel = params.channel;
    CustomImage image;
    
    if (params.singleCells)
    {
        if (cellImages.empty()) { return false; }
        if (scale == -1) { scale = 4; }
        image = cellImages[imageNumber];
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
}
*/

void Plotting::plot()
{
    Mat image;
    int pltType = static_cast<int>(params.plotFeatType);
    if (pltType == 0 || pltType == 1 || pltType == 2) // single cell
    {
        int cellNumber = params.showNumber;
        if (cellNumber < 0 || cellNumber >= cellImages.size())
        {
            return;
        }
        
        double randomScale;
        string randomName;
        //getImageToShow(image, randomName, randomScale);
    }

    if (params.plotFeatType == plotFeatureType::actFibersIntensity)
    {
        vector<double> resultingAngles;
        analysisFiberDir.analyseWithPCA(image, resultingAngles);

        plotData(resultingAngles, true);
    }
    else if (params.plotFeatType == plotFeatureType::actFibersOptThresh)
    {
        vector<double> resultingAngles;

        Mat thresholdedImage = image.clone();
        analysisFiberDir.getPCAoptThresholdedImage(thresholdedImage);
        analysisFiberDir.analyseWithPCA(image, resultingAngles, thresholdedImage);

        plotData(resultingAngles, true);
    }
    else if (params.plotFeatType == plotFeatureType::actFibersBoth)
    {
        vector<double> resultingAnglesoptThresh, resultingAnglesIntensity;
        Mat imageCopy = image.clone();
        Mat thresholdedImage = image.clone();

        params.PCAsquareLength = params.goodIntensityPCAsquareLength;
        params.PCAminEigValRatio = params.goodIntensityPCAminEigValRatio;
        analysisFiberDir.analyseWithPCA(image, resultingAnglesIntensity);

        params.PCAsquareLength = params.goodOptPCAsquareLength;
        params.PCAminEigValRatio = params.goodOptPCAminEigValRatio;
        analysisFiberDir.getPCAoptThresholdedImage(thresholdedImage);
        analysisFiberDir.analyseWithPCA(imageCopy, resultingAnglesoptThresh, thresholdedImage);

        plotData(resultingAnglesoptThresh, true, resultingAnglesIntensity);
    }
    else
    {
        vector<double> plottingData;
        for (int i = 0; i < cellImages.size(); i++)
        {
            plottingData.push_back(cellImages[i].getQuantity(params.plotFeatType));
        }
        plotData(plottingData, false);
    }
}
