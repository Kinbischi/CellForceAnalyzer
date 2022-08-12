#include "WindowMain.h"

#include <string>
#include <set>
#include <cassert>
#include <fstream>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include <QMessagebox>
#include <QTablewidget>

#include "WindowChannels.h"
#include "helperFunctions.h"


using namespace std;
using namespace cv;


WindowMain::WindowMain(QWidget* parent)
    : QMainWindow(parent), m_analysis(m_params), m_analysisFiberDir(m_params), m_plotting(m_params, m_analysisFiberDir, m_cellImages)
{
    ui.setupUi(this);

    connect(ui.spinBox_showImage, SIGNAL(valueChanged(int)), this, SLOT(imageNumberChanged(int)));

    connect(ui.radioButton_cellArrays, SIGNAL(toggled(bool)), this, SLOT(radioButtonArraysChanged(bool)));
    connect(ui.radioButton_cellArraysWithBoxes, SIGNAL(toggled(bool)), this, SLOT(radioButtonArraysChanged(bool)));
    connect(ui.radioButton_singleCells, SIGNAL(toggled(bool)), this, SLOT(radioButtonCellsChanged(bool)));
    connect(ui.radioButton_deletedCells, SIGNAL(toggled(bool)), this, SLOT(radioButtonDeletedCellsChanged(bool)));

    //it grabs the channels (order of input channels) as default values from the .ui file
    WindowChannels channelsWindow;
    m_channels = channelsWindow.set_ChannelTypes();
    
}

void WindowMain::updateParameters()
{
    m_params.showNumber = ui.spinBox_showImage->value();

    m_params.threshType = static_cast<thresholdingType>(ui.comboBox_thresholding->currentIndex());
    m_params.dispType = static_cast<displayType>(ui.comboBox_display->currentIndex());
    m_params.plotFeatType = static_cast<plotFeatureType>(ui.comboBox_showPlot->currentIndex());

    //PCA parameters
    m_params.PCAsquareLength = ui.spinBox_squareLengthPCA->value();
    m_params.PCAminEigValRatio = ui.doubleSpinBox_minEigValRatioPCA->value();
    m_params.suppressLowEigValRatioSquares = ui.checkBox_suppressLowEigValRatioSquares->isChecked();

    //Focal adhesion detection parameters
    m_params.highThreshCanny = ui.spinBox_highThreshCanny->value();
    m_params.FAminDist = ui.spinBox_minDist->value();
    m_params.FAminCircleConfidence = ui.doubleSpinBox_circleConfidence->value();
    m_params.FAdp = ui.doubleSpinBox_dpResolution->value();
    
    if (ui.checkBox_roundnessAnalysis->isChecked() || ui.checkBox_pcaAnalysis->isChecked() 
        || ui.checkBox_edgeDetection->isChecked() || ui.checkBox_FAdetection->isChecked())
    { m_params.withAnalysis = true; }
    else 
    { m_params.withAnalysis = false; }

}

void WindowMain::imageNumberChanged(int index)
{
    m_imageNumber_show = index;
    if (ui.radioButton_singleCells->isChecked())
    {
        if (!m_cellImages.empty())
        {
            updateCellAnalysisTable(m_cellImages[index]);
        }
    }
}

void WindowMain::radioButtonArraysChanged(bool index)
{
    if (index)
    {
        if (m_arrayImages.empty()) { ui.spinBox_showImage->setMaximum(0); }
        else { ui.spinBox_showImage->setMaximum(m_arrayImages.size() - 1); }
        ui.tableWidget_cell->clearContents();
    }
}

void WindowMain::radioButtonCellsChanged(bool index)
{
    if (index)
    {
        if (m_cellImages.empty()) { ui.spinBox_showImage->setMaximum(0); }
        else { ui.spinBox_showImage->setMaximum(m_cellImages.size() - 1);
            updateCellAnalysisTable(m_cellImages[m_imageNumber_show]); }
    }
}

void WindowMain::radioButtonDeletedCellsChanged(bool index)
{
    if (index)
    {
        if (m_deletedCellImages.empty()) { ui.spinBox_showImage->setMaximum(0); }
        else { ui.spinBox_showImage->setMaximum(m_deletedCellImages.size() - 1); }
        ui.tableWidget_cell->clearContents();
    }
}

void WindowMain::updateGeneralTable()
{
    QTableWidgetItem* newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_arrayImages.size())));
    ui.tableWidget_all->setItem(0, 0, newItem);
    QTableWidgetItem* newItem2 = new QTableWidgetItem(QString::fromStdString(to_string(m_cellImages.size())));
    ui.tableWidget_all->setItem(1, 0, newItem2);
    QTableWidgetItem* newItem3 = new QTableWidgetItem(QString::fromStdString(to_string(m_deletedCellImages.size())));
    ui.tableWidget_all->setItem(2, 0, newItem3);
}

void WindowMain::updateCellAnalysisTable(Cell cell)
{
    // no delete is needed, as the QTableWidget takes ownership of the QTableWidgetItem and automized deletion
    QTableWidgetItem* newItem;
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.nucleus_area)));
    ui.tableWidget_cell->setItem(0, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.nucleus_circularity)));
    ui.tableWidget_cell->setItem(1, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.nucleus_roundness)));
    ui.tableWidget_cell->setItem(2, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.actin_area)));
    ui.tableWidget_cell->setItem(3, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.actin_density)));
    ui.tableWidget_cell->setItem(4, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.actin_fiberAlignment)));
    ui.tableWidget_cell->setItem(5, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(100*cell.yap_inNucleus)+" %"));
    ui.tableWidget_cell->setItem(6, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.actin_mainAngle)));
    ui.tableWidget_cell->setItem(7, 0, newItem);
}

void WindowMain::updateGeneralAnalysisTable()
{
    QTableWidgetItem* newItem;
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.nucleus_area)));
    ui.tableWidget_all->setItem(4, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.nucleus_circularity)));
    ui.tableWidget_all->setItem(5, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.nucleus_roundness)));
    ui.tableWidget_all->setItem(6, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.actin_area)));
    ui.tableWidget_all->setItem(7, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.actin_density)));
    ui.tableWidget_all->setItem(8, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.actin_fiberAlignment)));
    ui.tableWidget_all->setItem(9, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.yap_inNucleus)));
    ui.tableWidget_all->setItem(10, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.actin_mainAngle)));
    ui.tableWidget_all->setItem(11, 0, newItem);
}


void WindowMain::on_pushButton_channels_clicked()
{
    WindowChannels channelsWindow;
    channelsWindow.setModal(true);
    int dialogCode=channelsWindow.exec();
    if (dialogCode == QDialog::Accepted)
    {
        m_channels = channelsWindow.set_ChannelTypes();
    }
}


void WindowMain::on_pushButton_test_clicked()
{
    Mat img = imread(m_inpDir + "300321_PEGvsGELMAlif.lif - Image004 (RGB).tif", IMREAD_UNCHANGED);
    Rect r(1300, 780, 500, 620);
    rectangle(img,r, Scalar(128, 255, 0));
    vector<vector<double>> hi{ { 1,2 }, {3,4} };
   
    //help::showWindow(img,0.5);
  
    vector<vector<int>> testVectors
    { 
        { 4,4,4,4,4,4,4,4,4,4,4,4 }, //0
        { 0,0,0,0,8,0,0,0,0,0,0,0},
        {0, 0, 2, 2, 2, 2, 4, 4, 4, 4, 0, 0}, //2
        { 0,0,4,4,4,4,8,8,8,8,0,0 }, //3
        { 0,4,4,8,8,0 }, //4
        //{ 0,0,0,0,8,8,0,0,0,0,0,0},
        //{ 0,0,8,0,0,0},
        //{0,0,1,1,0,0,30,30,2,2,0,0}, //2
        //{ 0,0,2,2,0,0,60,60,4,4,0,0 }, //3
        //{ 0,2,0,60,4,0 }, //4
        { 4,6,5,7,4,5,6,6,4,4,5,6 }, //5
        { 3,6,2,8,4,3,6,3,4,4,5,8 }, //6
        { 1,2,3,4,5,6,7,8,9,10,11,12 }, //7
        { 3,2,4,5,1,11,5,3,5,4,3,6 }, //8
        { 1,2,4,2,1,15,2,3,1,1,3,1 }, //9
        { 31,22,43,2,1,150,2,30,13,14,32,11 } //10  
    };

    vector<double> x = {0,15,30,45,60,75,90,105,120,135,150,165};
    vector<double> x2 = { 0,30,60,90,120,150 };

    auto y0 = testVectors[0];
    auto y1 = testVectors[1];
    auto y2 = testVectors[2];
    auto y3 = testVectors[3];
    auto y4 = testVectors[4];
    m_params.plotFeatType = plotFeatureType::actFibersOptThresh;

    /*
    plotSomething(x, y0, plotFeatureType::actFibersOptThresh, "Example 1");
    plotSomething(x, y1, plotFeatureType::actFibersOptThresh, "Example 2");
    plotSomething(x, y2, plotFeatureType::actFibersOptThresh, "Example 3");
    plotSomething(x, y3, plotFeatureType::actFibersOptThresh, "Example 4");
    plotSomething(x2, y4, plotFeatureType::actFibersOptThresh, "Example 5");
    */

    vector<double> testRes_2_0 = m_analysisFiberDir.testYvecs(testVectors, 2, 0);
    vector<double> testRes_2_1 = m_analysisFiberDir.testYvecs(testVectors, 2, 1);
    vector<double> testRes_2_2 = m_analysisFiberDir.testYvecs(testVectors, 2, 2);
    vector<double> testRes_3_0 = m_analysisFiberDir.testYvecs(testVectors, 3, 0);
    vector<double> testRes_3_1 = m_analysisFiberDir.testYvecs(testVectors, 3, 1);
    vector<double> testRes_3_2 = m_analysisFiberDir.testYvecs(testVectors, 4, 2);
    
    m_channels = m_preprocess.loadNiceCellImages(m_arrayImages, m_cellImages, m_inpDir);
    updateGeneralTable();
    
}

void WindowMain::on_pushButton_loadImages_clicked()
{
    updateParameters();

    bool hasNucleus, hasBrightfield, hasActin;
    hasNucleus = hasBrightfield = hasActin = false;
    for (auto channel : m_channels)
    {
        if (channel == channelType::actin) {hasActin = true;}
        //if (channel == channelType::brightfield) { hasBrightfield = true; }
        //if (channel == channelType::nucleus) { hasNucleus = true; }
    }
    if (!(hasActin )) //&& hasBrightfield && hasNucleus
    {
        QMessageBox::information(this, "Sweetheart", "You need to have an actin channel. Otherwise this analysis tool kind of makes little sense^^");
        return;
    }

    string inName = ui.lineEdit_inpDir->text().toStdString();
    if (!inName.empty())
    {
        m_inpDir = help::copiedDirectoryToNiceString(inName);
    }

    float confThreshold = ui.doubleSpinBox_confThresh->value();

    m_preprocess.loadImages(m_arrayImages, m_inpDir, m_channels);

    //TODO introduce images already are cells, with/without yolo,...
    m_preprocess.applyYolo(m_arrayImages, m_cellImages, m_arrayImages_withYoloBoxes, confThreshold);

    
    bool isDead;
    int failedAnalysisCells=0;
    vector<int> cellsToDelete;
    for (int i=0; i<m_cellImages.size(); i++)
    {
        isDead = m_analysis.isDeadCell(m_cellImages[i]);
        if (isDead) 
        {
            cellsToDelete.push_back(i);
            m_deletedCellImages.push_back(m_cellImages[i]);
            continue; 
        }

        try 
        {
            m_analysis.analyseCell(m_cellImages[i], m_channels);//TODO: what analysis did not work => include e.g. the assert(isContour),...
        }
        catch (...)
        {
            failedAnalysisCells++;

            cellsToDelete.push_back(i);
            m_deletedCellImages.push_back(m_cellImages[i]);
        }
    }

    for(int j = cellsToDelete.size()-1; j>-1; j--)
    {
        m_cellImages.erase(m_cellImages.begin()+cellsToDelete[j]);
    }

    m_averageAllCells = m_preprocess.getAverageProperties(m_cellImages);
    
    updateGeneralTable();
    
    if (m_cellImages.size()>0)
    {
        updateGeneralAnalysisTable();
        
        radioButtonCellsChanged(ui.radioButton_singleCells->isChecked());
    }

    if(failedAnalysisCells>0)
    {
        auto temp = to_string(failedAnalysisCells) + " cells ended up in deleted cells";
        QMessageBox::information(this, "The analysis failed for some cells", temp.c_str());
    }
}

bool WindowMain::getImageToShow(Mat& outImg, string& name, double& scale)
{
    if (m_arrayImages.empty())
    {
        QMessageBox::information(this, "Good Morning", "Read in Images!?");
        return false;
    }
    
    int imageNumber = m_params.showNumber;
    channelType channel = static_cast<channelType>(ui.comboBox_showImage->currentIndex());
    CustomImage image;

    if (ui.radioButton_cellArrays->isChecked())
    {
        if (scale == -1) { scale = 0.4; }
        image = m_arrayImages[imageNumber];
        outImg = image.getChannel(channel);
        name = image.getName(channel);
    }
    if (ui.radioButton_singleCells->isChecked())
    {
        if (m_cellImages.empty()) { return false; }
        if (scale == -1) { scale = 4; }
        image = m_cellImages[imageNumber];
        outImg = image.getChannel(channel);
        name = image.getName(channel);
    }
    if (ui.radioButton_deletedCells->isChecked())
    {
        if (m_deletedCellImages.empty()) { return false; }
        if (scale == -1) { scale = 4; }
        image = m_deletedCellImages[imageNumber];
        outImg = image.getChannel(channel);
        name = image.getName(channel);
    }
    if (ui.radioButton_cellArraysWithBoxes->isChecked())
    {
        if (scale == -1) { scale = 0.4; }
        image = m_arrayImages[imageNumber]; // only for name
        outImg = m_arrayImages_withYoloBoxes[imageNumber];
        name = image.getName();
    }

    if (outImg.empty())
    {
        QMessageBox::information(this, "Good Morning", "This channel was not loaded!");
        return false;
    }

    outImg = outImg.clone();
}


void WindowMain::on_pushButton_showImage_clicked()
{
    updateParameters();
    double scale = ui.doubleSpinBox_scale->value();
    bool replacing = ui.checkBox_replace->isChecked();
    Mat image;
    string name;

    //load image
    bool successful = getImageToShow(image, name, scale);
    if (!successful) { return; }
    
    if (image.channels() == 1)//this can only be executed for 1 channel images
    {
        //thresholdings
        Mat thresholdedImage = image;

        switch (m_params.threshType)
        {
        case thresholdingType::manual:

            help::thresh(thresholdedImage, ui.spinBox_thresholdManually->value());
            name = name + "_threshManually" + to_string(ui.spinBox_thresholdManually->value());
            break;

        case thresholdingType::otsu:
            help::thresh(thresholdedImage);
            name = name + "_threshOtsu";
            break;

        case thresholdingType::squarePCAoptimizedThresh:
            m_analysisFiberDir.getPCAoptThresholdedImage(thresholdedImage);
            name = name + "_threshPCAopt" + "_Length" + to_string(m_params.PCAsquareLength) + "_minEigRatio" + to_string(m_params.PCAminEigValRatio);
            break;
        }

        
        if ((m_params.withAnalysis && m_params.dispType == displayType::thresholded) || (!m_params.withAnalysis && m_params.threshType != thresholdingType::None))
        {
            image = thresholdedImage;
        }


        //apply analysis
        bool analysisSuccessful = false;
        if (ui.checkBox_roundnessAnalysis->isChecked())
        {
            analysisSuccessful = m_analysis.analyseShape(image);
            name = name + "_roundAnalysed";
        }
        else if (ui.checkBox_pcaAnalysis->isChecked())
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
        else if (ui.checkBox_edgeDetection->isChecked())
        {
            m_analysis.edgeDetectionCanny(image);
            name = name + "_edge" + "_highthresh" + to_string(m_params.highThreshCanny);
        }
        else if (ui.checkBox_FAdetection->isChecked())
        {
            m_analysis.focalAdhesiondetection(image);
            name = name + "_circles" + "_highthresh" + to_string(m_params.highThreshCanny) + "_minCircleConf" + 
                to_string(m_params.FAminCircleConfidence) +"_minDist" + to_string(m_params.FAminDist) + "_dp" + to_string(m_params.FAdp);
        }
    }
    
    if (!m_params.withAnalysis)//TODO: check what happens for unscaled 16 bit image that were analyzed => probs not scaled ;(
    {
        help::scaleData(image);
    }

    string windowName = replacing ? "Image" : name;

    help::showWindow(image, scale, windowName);
}

void WindowMain::on_pushButton_showPlot_clicked()
{
    updateParameters();

    m_plotting.plot();
}

void WindowMain::on_pushButton_conductAnalysisOnSingleCell_clicked() 
{
    updateParameters();
    m_analysis.analyseCell(m_cellImages[m_params.showNumber], m_channels);
    updateCellAnalysisTable(m_cellImages[m_params.showNumber]);
}

void WindowMain::on_pushButton_conductAnalysisOnAllCells_clicked()
{
    updateParameters();
    
    for (int i =0; i<m_cellImages.size(); i++) 
    {
        m_analysis.analyseCell(m_cellImages[i], m_channels);
    }

    m_averageAllCells = m_preprocess.getAverageProperties(m_cellImages);
    updateGeneralAnalysisTable();
    updateCellAnalysisTable(m_cellImages[m_params.showNumber]);
}

void WindowMain::on_pushButton_writeOut_clicked()
{
    if (m_arrayImages.empty())
    {
        QMessageBox::information(this, "Good Morning", "Read in Images!?");
        return;
    }

    string outName = ui.lineEdit_outpDir->text().toStdString();
    if (!outName.empty())
    {
        m_outpDir = help::copiedDirectoryToNiceString(outName);
    }
    
    writeAnalysedDataToFile();

    if (ui.checkBox_cellArrays->isChecked())
    {
        for (auto entry : m_arrayImages)
        {
            Mat out = entry.createRGBimage();
            help::scaleData(out);
            if (out.depth()==CV_8U)
            {
                imwrite(m_outpDir + entry.m_name + ".jpg", out);
            }
            if (out.depth()==CV_16U)
            {
                imwrite(m_outpDir + entry.m_name + ".png", out);
            }
        }
    }
    if (ui.checkBox_singleCells->isChecked())
    {
        for (auto entry : m_cellImages)
        {
            Mat out = entry.createRGBimage();
            help::scaleData(out);
            if (out.depth() == CV_8U)
            {
                imwrite(m_outpDir + entry.m_name + ".jpg", out);
            }
            if (out.depth() == CV_16U)
            {
                imwrite(m_outpDir + entry.m_name + ".png", out);
            }
        }
    }
    if (ui.checkBox_cellArraysWithBoxes->isChecked())
    {
        for (int i=0;i< m_arrayImages_withYoloBoxes.size();i++)
        {
            auto x = m_arrayImages_withYoloBoxes[i].depth();
            imwrite(m_outpDir + m_arrayImages[i].m_name + "_yoloAnalysed.jpg", m_arrayImages_withYoloBoxes[i]);
        }
    }
}

void WindowMain::writeAnalysedDataToFile()
{
    
    string filename = "dataAnalysis.csv";
    
    std::ofstream myfile(m_outpDir + filename);

    myfile << ",Nucleus Area,Nucleus Circularity,Nucleus Roundness,Actin Area,Actin Density, Actin maxSpreadLength,Actin PCA direction,Yap in Nucleus" << "\n";

    for (auto cell : m_cellImages)
    {
        myfile << cell.m_name << "," << cell.nucleus_area << "," << cell.nucleus_circularity<< "," 
            << cell.nucleus_roundness << "," << cell.actin_area << "," << cell.actin_density << "," 
            << cell.actin_maxLength << "," << cell.actin_mainAngle << "," << cell.yap_inNucleus << "\n";
    }
    
    myfile <<"Averages," << m_averageAllCells.nucleus_area << "," << m_averageAllCells.nucleus_circularity << "," 
        << m_averageAllCells.nucleus_roundness << "," << m_averageAllCells.actin_area << "," << m_averageAllCells.actin_density << "," 
        << m_averageAllCells.actin_maxLength << "," << m_averageAllCells.actin_mainAngle << "," << m_averageAllCells.yap_inNucleus;

}

