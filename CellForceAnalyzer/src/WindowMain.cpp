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
    : QMainWindow(parent), m_preprocess(m_params,m_data), m_analysis(m_params), m_analysisFiberDir(m_params), m_plotting(m_params, m_data, m_analysisFiberDir),
    m_display(m_params, m_data, m_analysis, m_analysisFiberDir)
{
    ui.setupUi(this);

    connect(ui.spinBox_showImage, SIGNAL(valueChanged(int)), this, SLOT(imageNumberChanged(int)));

    connect(ui.radioButton_cellArrays, SIGNAL(toggled(bool)), this, SLOT(radioButtonArraysChanged(bool)));
    connect(ui.radioButton_cellArraysWithBoxes, SIGNAL(toggled(bool)), this, SLOT(radioButtonArraysChanged(bool)));
    connect(ui.radioButton_singleCells, SIGNAL(toggled(bool)), this, SLOT(radioButtonCellsChanged(bool)));
    connect(ui.radioButton_deletedCells, SIGNAL(toggled(bool)), this, SLOT(radioButtonDeletedCellsChanged(bool)));

    //it grabs the channels (order of input channels) as default values from the .ui file
    WindowChannels channelsWindow;
    m_data.channels = channelsWindow.set_ChannelTypes();
    
}

void WindowMain::updateParameters()
{
    m_params.showNumber = ui.spinBox_showImage->value();
    m_params.scaleFactor = ui.doubleSpinBox_scale->value();
    m_params.replacingMode = ui.checkBox_replace->isChecked();

    //loading
    m_params.loadAsCells = ui.checkBox_loadAsCells->isChecked();
    m_params.loadthreeChannels = ui.checkBox_3channelsPerImage->isChecked();

    m_params.channel = static_cast<channelType>(ui.comboBox_showImage->currentIndex());

    m_params.cellArrays = ui.radioButton_cellArrays->isChecked();
    m_params.singleCells = ui.radioButton_singleCells->isChecked();
    m_params.deletedCells = ui.radioButton_deletedCells->isChecked();
    m_params.cellArraysWithBoxes = ui.radioButton_cellArraysWithBoxes->isChecked();

    m_params.threshType = static_cast<thresholdingType>(ui.comboBox_thresholding->currentIndex());
    m_params.dispType = static_cast<displayType>(ui.comboBox_display->currentIndex());
    m_params.plotFeatType = static_cast<plotFeatureType>(ui.comboBox_showPlot->currentIndex());

    //analysis conducted
    m_params.fiberPCA = ui.radioButton_fiberPCA->isChecked();
    m_params.edgeDetection = ui.radioButton_edgeDetection->isChecked();
    m_params.FAdetection = ui.radioButton_FAdetection->isChecked();
    m_params.variousAnalysis = ui.radioButton_variousAnalysis->isChecked();

    //thresholdinng parameters
    m_params.manualThreshold = ui.spinBox_thresholdManually->value();

    //PCA parameters
    m_params.PCAsquareLength = ui.spinBox_squareLengthPCA->value();
    m_params.PCAminEigValRatio = ui.doubleSpinBox_minEigValRatioPCA->value();
    m_params.suppressLowEigValRatioSquares = ui.checkBox_suppressLowEigValRatioSquares->isChecked();

    //Focal adhesion detection parameters
    m_params.highThreshCanny = ui.spinBox_highThreshCanny->value();
    m_params.FAminDist = ui.spinBox_minDist->value();
    m_params.FAminCircleConfidence = ui.doubleSpinBox_circleConfidence->value();
    m_params.FAdp = ui.doubleSpinBox_dpResolution->value();
    
    if (ui.radioButton_variousAnalysis->isChecked() || ui.radioButton_fiberPCA->isChecked()
        || ui.radioButton_edgeDetection->isChecked() || ui.radioButton_FAdetection->isChecked())
    { m_params.withAnalysis = true; }
    else 
    { m_params.withAnalysis = false; }

}

void WindowMain::imageNumberChanged(int index)
{
    m_imageNumber_show = index;
    if (ui.radioButton_singleCells->isChecked())
    {
        if (!m_data.cellImages.empty())
        {
            updateCellAnalysisTable(m_data.cellImages[index]);
        }
    }
}

void WindowMain::radioButtonArraysChanged(bool index)
{
    if (index)
    {
        if (m_data.arrayImages.empty()) { ui.spinBox_showImage->setMaximum(0); }
        else { ui.spinBox_showImage->setMaximum(m_data.arrayImages.size() - 1); }
        ui.tableWidget_cell->clearContents();
    }
}

void WindowMain::radioButtonCellsChanged(bool index)
{
    if (index)
    {
        if (m_data.cellImages.empty()) { ui.spinBox_showImage->setMaximum(0); }
        else { ui.spinBox_showImage->setMaximum(m_data.cellImages.size() - 1);
            updateCellAnalysisTable(m_data.cellImages[m_imageNumber_show]); }
    }
}

void WindowMain::radioButtonDeletedCellsChanged(bool index)
{
    if (index)
    {
        if (m_data.deletedCellImages.empty()) { ui.spinBox_showImage->setMaximum(0); }
        else { ui.spinBox_showImage->setMaximum(m_data.deletedCellImages.size() - 1); }
        ui.tableWidget_cell->clearContents();
    }
}

void WindowMain::updateGeneralTable()
{
    QTableWidgetItem* newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_data.arrayImages.size())));
    ui.tableWidget_all->setItem(0, 0, newItem);
    QTableWidgetItem* newItem2 = new QTableWidgetItem(QString::fromStdString(to_string(m_data.cellImages.size())));
    ui.tableWidget_all->setItem(1, 0, newItem2);
    QTableWidgetItem* newItem3 = new QTableWidgetItem(QString::fromStdString(to_string(m_data.deletedCellImages.size())));
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
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_data.averageAllCells.nucleus_area)));
    ui.tableWidget_all->setItem(4, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_data.averageAllCells.nucleus_circularity)));
    ui.tableWidget_all->setItem(5, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_data.averageAllCells.nucleus_roundness)));
    ui.tableWidget_all->setItem(6, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_data.averageAllCells.actin_area)));
    ui.tableWidget_all->setItem(7, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_data.averageAllCells.actin_density)));
    ui.tableWidget_all->setItem(8, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_data.averageAllCells.actin_fiberAlignment)));
    ui.tableWidget_all->setItem(9, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_data.averageAllCells.yap_inNucleus)));
    ui.tableWidget_all->setItem(10, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_data.averageAllCells.actin_mainAngle)));
    ui.tableWidget_all->setItem(11, 0, newItem);
}


void WindowMain::on_pushButton_channels_clicked()
{
    WindowChannels channelsWindow;
    channelsWindow.setModal(true);
    int dialogCode=channelsWindow.exec();
    if (dialogCode == QDialog::Accepted)
    {
        m_data.channels = channelsWindow.set_ChannelTypes();
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
    
    m_data.channels = m_preprocess.loadNiceCellImages(m_data.arrayImages, m_data.cellImages, m_inpDir);
    updateGeneralTable();
    
}

void WindowMain::on_pushButton_loadImages_clicked()
{
    updateParameters();

    bool hasNucleus, hasBrightfield, hasActin;
    hasNucleus = hasBrightfield = hasActin = false;
    for (auto channel : m_data.channels)
    {
        if (channel == channelType::actin) {hasActin = true;}
        //if (channel == channelType::brightfield) { hasBrightfield = true; }
        //if (channel == channelType::nucleus) { hasNucleus = true; }
    }
    if (!(hasActin )) //&& hasBrightfield && hasNucleus
    {
        QMessageBox::information(this, "Sweetheart", "You need to have an actin channel. Otherwise this analysis tool makes little sense^^");
        return;
    }

    string inName = ui.lineEdit_inpDir->text().toStdString();
    if (!inName.empty())
    {
        m_inpDir = help::copiedDirectoryToNiceString(inName);
    }

    float confThreshold = ui.doubleSpinBox_confThresh->value();

    m_preprocess.loadImages(m_data.arrayImages, m_inpDir, m_data.channels);

    //TODO introduce images already are cells, with/without yolo,...
    m_preprocess.applyYolo(m_data.arrayImages, m_data.cellImages, m_data.arrayImages_withYoloBoxes, confThreshold);

    
    bool isDead;
    int failedAnalysisCells=0;
    vector<int> cellsToDelete;
    for (int i=0; i<m_data.cellImages.size(); i++)
    {
        isDead = m_analysis.isDeadCell(m_data.cellImages[i]);
        if (isDead) 
        {
            cellsToDelete.push_back(i);
            m_data.deletedCellImages.push_back(m_data.cellImages[i]);
            continue; 
        }

        try 
        {
            m_analysis.analyseCell(m_data.cellImages[i], m_data.channels);//TODO: what analysis did not work => include e.g. the assert(isContour),...
        }
        catch (...)
        {
            failedAnalysisCells++;

            cellsToDelete.push_back(i);
            m_data.deletedCellImages.push_back(m_data.cellImages[i]);
        }
    }

    for(int j = cellsToDelete.size()-1; j>-1; j--)
    {
        m_data.cellImages.erase(m_data.cellImages.begin()+cellsToDelete[j]);
    }

    m_data.averageAllCells = m_preprocess.getAverageProperties(m_data.cellImages);
    
    updateGeneralTable();
    
    if (m_data.cellImages.size()>0)
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

void WindowMain::getImageWorked(int successful)
{
    switch (successful)
    {
    case 1:
        QMessageBox::information(this, "Good Morning", "Read in Images!?");
        break;
    case 2:
        QMessageBox::information(this, "Good Morning", "This channel was not loaded!");
        break;
    case 3:
        QMessageBox::information(this, "Good Morning", "No such cells available");
        break;
    }
}

void WindowMain::on_pushButton_showImage_clicked()
{
    updateParameters();

    int successful = m_display.prepareAnalysisToShow();

    getImageWorked(successful);
}

void WindowMain::on_pushButton_showPlot_clicked()
{
    updateParameters();

    int successful = m_plotting.plot();

    getImageWorked(successful);
}

void WindowMain::on_pushButton_conductAnalysisOnSingleCell_clicked() 
{
    updateParameters();
    m_analysis.analyseCell(m_data.cellImages[m_params.showNumber], m_data.channels);
    updateCellAnalysisTable(m_data.cellImages[m_params.showNumber]);
}

void WindowMain::on_pushButton_conductAnalysisOnAllCells_clicked()
{
    updateParameters();
    
    for (int i =0; i<m_data.cellImages.size(); i++) 
    {
        m_analysis.analyseCell(m_data.cellImages[i], m_data.channels);
    }

    m_data.averageAllCells = m_preprocess.getAverageProperties(m_data.cellImages);
    updateGeneralAnalysisTable();
    updateCellAnalysisTable(m_data.cellImages[m_params.showNumber]);
}

void WindowMain::on_pushButton_writeOut_clicked()
{
    if (m_data.arrayImages.empty())
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
        for (auto entry : m_data.arrayImages)
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
        for (auto entry : m_data.cellImages)
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
        for (int i=0;i< m_data.arrayImages_withYoloBoxes.size();i++)
        {
            auto x = m_data.arrayImages_withYoloBoxes[i].depth();
            imwrite(m_outpDir + m_data.arrayImages[i].m_name + "_yoloAnalysed.jpg", m_data.arrayImages_withYoloBoxes[i]);
        }
    }
}

void WindowMain::writeAnalysedDataToFile()
{
    
    string filename = "dataAnalysis.csv";
    
    std::ofstream myfile(m_outpDir + filename);

    myfile << ",Nucleus Area,Nucleus Circularity,Nucleus Roundness,Actin Area,Actin Density, Actin maxSpreadLength,Actin PCA direction,Yap in Nucleus" << "\n";

    for (auto cell : m_data.cellImages)
    {
        myfile << cell.m_name << "," << cell.nucleus_area << "," << cell.nucleus_circularity<< "," 
            << cell.nucleus_roundness << "," << cell.actin_area << "," << cell.actin_density << "," 
            << cell.actin_maxLength << "," << cell.actin_mainAngle << "," << cell.yap_inNucleus << "\n";
    }
    
    myfile <<"Averages," << m_data.averageAllCells.nucleus_area << "," << m_data.averageAllCells.nucleus_circularity << "," 
        << m_data.averageAllCells.nucleus_roundness << "," << m_data.averageAllCells.actin_area << "," << m_data.averageAllCells.actin_density << "," 
        << m_data.averageAllCells.actin_maxLength << "," << m_data.averageAllCells.actin_mainAngle << "," << m_data.averageAllCells.yap_inNucleus;

}

