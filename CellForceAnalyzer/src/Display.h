#pragma once
#include "Cell.h"
#include "ParametersFromUI.h"
#include "Analysis.h"
#include "AnalysisFiberDirection.h"


class Display
{
public:
    Display(ParametersFromUI&, std::vector<CustomImage>&, std::vector<cv::Mat>&, 
        std::vector<Cell>&, std::vector<Cell>&, Analysis&, AnalysisFiberDirection&);

    int getImageToShow(cv::Mat&, std::string&, double&);
    cv::Mat thresholdImage(cv::Mat, std::string&);
    int prepareAnalysisToShow();

private:
    ParametersFromUI& m_params;

    std::vector<CustomImage>& m_arrayImages;
    std::vector<cv::Mat>& m_arrayImages_withYoloBoxes;
    std::vector<Cell>& m_cellImages;
    std::vector<Cell>& m_deletedCellImages;

    Analysis& m_analysis;
    AnalysisFiberDirection& m_analysisFiberDir;

};

