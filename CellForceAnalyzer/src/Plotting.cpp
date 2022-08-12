#include "Plotting.h"

using namespace std;
using namespace cv;


#ifndef _DEBUG

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

#endif


Plotting::Plotting(ParametersFromUI& p, AnalysisFiberDirection& a, std::vector<Cell>& c) :params(p), analysisFiberDir(a), cellImages(c) {};



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

    vector<string> labels = Cell::getPlotLegendNames(params.plotFeatType);

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

    plt::xlabel(Cell::getPlotXLabel(params.plotFeatType));
    plt::ylabel(Cell::getPlotYLabel(params.plotFeatType));
    plt::title(Cell::getPlotTitle(params.plotFeatType));
    plt::show();
#endif
}

#ifndef _DEBUG
void plotSomething(vector<double> x, vector<double> y, plotFeatureType pltType, string title)
{
    plt::figure();

    plt::plot(x, y, { {"marker", "o"}, {"linestyle", "--"} });

    plt::xlabel(Cell::getPlotXLabel(pltType));
    plt::ylabel(Cell::getPlotYLabel(pltType));
    plt::ylim(0, 10);
    plt::title(title);
    //plt::savefig(title +".pdf");
    plt::show();
}
#endif



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
