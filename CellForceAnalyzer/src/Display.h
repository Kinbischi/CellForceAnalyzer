#pragma once
#include "Cell.h"
#include "ParametersUI.h"
#include "Analysis.h"
#include "AnalysisFiberDirection.h"
#include "dataContainer.h"


class Display
{
public:

    Display(ParametersUI&, dataContainer&, Analysis&, AnalysisFiberDirection&);

    int getImageToShow(cv::Mat&, std::string&, double&);
    void blurImg(cv::Mat&, std::string&);
    cv::Mat thresholdImage(cv::Mat, std::string&);
    int prepareAnalysisToShow();

private:

    ParametersUI& params;
    dataContainer& data;
    Analysis& analysis;
    AnalysisFiberDirection& analysisFiberDir;

};

