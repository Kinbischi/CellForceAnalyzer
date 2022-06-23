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

#ifndef _DEBUG

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

#endif

WindowMain::WindowMain(QWidget *parent)
    : QMainWindow(parent)
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

void WindowMain::imageNumberChanged(int index)
{
    m_imageNumber_show = index;
    if (ui.radioButton_singleCells->isChecked())
    {
        if (!m_cellImages.empty())
        {
            setCellTable(m_cellImages[index]);
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
            setCellTable(m_cellImages[m_imageNumber_show]); }
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

void WindowMain::setCellTable(Cell cell)
{
    // no delete is needed, as the QTableWidget takes ownership of the QTableWidgetItem and automized deletion
    QTableWidgetItem* newItem;
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.nucleus_circularity)));
    ui.tableWidget_cell->setItem(0, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.nucleus_roundness)));
    ui.tableWidget_cell->setItem(1, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.nucleus_area)));
    ui.tableWidget_cell->setItem(2, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(100*cell.yap_inNucleus)+" %"));
    ui.tableWidget_cell->setItem(3, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.actin_density)));
    ui.tableWidget_cell->setItem(4, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.actin_area)));
    ui.tableWidget_cell->setItem(5, 0, newItem);
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.actin_mainAngle)));
    ui.tableWidget_cell->setItem(6, 0, newItem);
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

void plotOccurrenceInData(vector<double> data)
{
    double min = *min_element(data.begin(), data.end());
    double max = *max_element(data.begin(), data.end());

    int plottingNumb = 12;
    double spacing = (max-min)/plottingNumb;
    double start = min - spacing;
    double end = max + spacing;

    vector<double> x;
    for (double i = start; i <= end; i = i + spacing)
    {
        x.push_back(i);
    }

    std::vector<int> y(x.size());
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
    #ifndef _DEBUG
    plt::plot(x, y);
    plt::show();
    #endif
}

void plotAngles(vector<double> angles)
{
    vector<int> x;
    double start = -15;
    double end = 180 + start;
    double spacing = 15;
    for (int i=start; i<end; i=i+spacing)
    {
        x.push_back(i);
    }

    std::vector<int> y(x.size());
    for (int i = 0;i<angles.size();i++)
    {
        for (int j = 0; j < x.size(); j++)
        {
            if (angles[i]>=x[j]-spacing/2 && angles[i]<x[j]+spacing/2)
            {
                y[j]++;
            }
        }
    }

    #ifndef _DEBUG
    plt::plot(x, y);
    plt::show();
    #endif
}

void plotDistribution(vector<Cell> cells, int type)
{
    vector<int> y;
    for (int i = 0; i < cells.size(); i++)
    {
        y.push_back(cells[i].nucleus_area);
    }
}
    
void WindowMain::on_pushButton_test_clicked()
{
    /*
    Mat img = imread(m_inpDir + "50x50_niches_RGD_3mM_2022_02_22__14_11_14_Maximum intensity projection_Filter.tif", IMREAD_UNCHANGED);
    Rect r(220, 400, 300, 300);
    std::vector<channelType> channelsOrder;
    channelsOrder.push_back(channelType::nucleus);
    channelsOrder.push_back(channelType::actin);
    channelsOrder.push_back(channelType::brightfield);
    */
    Mat img = imread(m_inpDir + "300321_PEGvsGELMAlif.lif - Image006-1 (RGB).tif", IMREAD_UNCHANGED);
    Rect r(1000, 1180, 500, 500);
    //rectangle(img,r, Scalar(128, 255, 0));
    //help::showWindow(img,0.5);
    
    std::vector<channelType> channelsOrder;
    channelsOrder.push_back(channelType::nucleus);
    channelsOrder.push_back(channelType::brightfield);
    channelsOrder.push_back(channelType::actin);

    Mat bgr[3];
    split(img, bgr);
    vector<Mat> matImages;
    matImages.push_back(bgr[0]);
    matImages.push_back(bgr[1]);
    matImages.push_back(bgr[2]);
    string name = "test";
    string name2 = "test_cell";
    CustomImage image(matImages, channelsOrder, name);
    CustomImage ce = image.cutImageOut(r, name2);
    Cell cell(ce);
    m_arrayImages.push_back(image);
    m_cellImages.push_back(cell);

    Mat actin = cell.getChannel(channelType::actin).clone();

    int maxArrowsCount = 0;
    int optimalThresholding;
    for (int i = 0;i<255;i++)
    {
        Mat actinThresh = actin.clone();
        cv::threshold(actinThresh, actinThresh, i, 255, THRESH_BINARY);
        
        vector<double> angles;
        m_analysis.analyseWithPCA(actinThresh, angles);

        int arrowsCount = angles.size();
        if (arrowsCount > maxArrowsCount)
        {
            maxArrowsCount = arrowsCount;
            optimalThresholding = i;
        }
    }
    
    Mat actinThresh = actin.clone();
    cv::threshold(actinThresh, actinThresh, optimalThresholding, 255, THRESH_BINARY);

    vector<double> random;
    m_analysis.analyseWithPCA(actinThresh, random);
    help::showWindow(actin, 2, "actin");
    help::showWindow(actinThresh, 2, "actinT");
}

void WindowMain::on_pushButton_loadImages_clicked()
{
    bool hasNucleus, hasBrightfield, hasActin;
    hasNucleus = hasBrightfield = hasActin = false;
    for (auto channel : m_channels)
    {
        if (channel == channelType::actin) {hasActin = true;}
        if (channel == channelType::brightfield) { hasBrightfield = true; }
        if (channel == channelType::nucleus) { hasNucleus = true; }
    }
    if (!(hasActin && hasBrightfield && hasNucleus)) 
    {
        QMessageBox::information(this, "Sweetheart", "You need to have a nucleus, actin and brightfield channel. At least one of them is not here^^");
        return;
    }
        
    string inName = ui.lineEdit_inpDir->text().toStdString();
    if (!inName.empty())
    {
        m_inpDir = help::copiedDirectoryToNiceString(inName);
    }

    float confThreshold = ui.doubleSpinBox_confThresh->value();
    float nmsThreshold = ui.doubleSpinBox_nonMaxThresh->value(); //TODO: get that out

    m_preprocess.loadImages(m_arrayImages, m_inpDir, m_channels);

    m_preprocess.applyYolo(m_arrayImages, m_cellImages, m_arrayImages_withYoloBoxes,confThreshold,nmsThreshold);

    
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
            m_analysis.analyseNucleusShape(m_cellImages[i]); //TODO: what analysis did not work => include e.g. the assert(isContour),...
            m_analysis.analyseYapInNucleus(m_cellImages[i]);
            m_analysis.analyseActin(m_cellImages[i]);
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

    m_averageAllCells = m_analysis.getAverageProperties(m_cellImages);
    
    QTableWidgetItem* newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_arrayImages.size())));
    ui.tableWidget_all->setItem(0, 0, newItem);
    QTableWidgetItem* newItem2 = new QTableWidgetItem(QString::fromStdString(to_string(m_cellImages.size())));
    ui.tableWidget_all->setItem(1, 0, newItem2);
    QTableWidgetItem* newItem3 = new QTableWidgetItem(QString::fromStdString(to_string(m_deletedCellImages.size())));
    ui.tableWidget_all->setItem(2, 0, newItem3);
    
    if (m_cellImages.size()>0)
    {
        QTableWidgetItem* newItem4 = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.nucleus_circularity)));
        ui.tableWidget_all->setItem(3, 0, newItem4);
        QTableWidgetItem* newItem5 = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.nucleus_roundness)));
        ui.tableWidget_all->setItem(4, 0, newItem5);
        
        radioButtonCellsChanged(ui.radioButton_singleCells->isChecked()); //TODO: change that!! not correct anymore
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
    
    int imageNumber = ui.spinBox_showImage->value();
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
        if (scale == -1) { scale = 4; }
        image = m_cellImages[imageNumber];
        outImg = image.getChannel(channel);
        name = image.getName(channel);
    }
    if (ui.radioButton_deletedCells->isChecked())
    {
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

// gescheites abfangen von else if/ if,... so dass nur eine kreuz geht? => auch bei anderer analyse original image anbieten?
// oder beide weglassen? 
void WindowMain::on_pushButton_showImage_clicked()
{
    double scale = ui.doubleSpinBox_scale->value();
    bool replacing = ui.checkBox_replace->isChecked();
    Mat image;
    string name;

    //load image
    bool successful = getImageToShow(image, name, scale);
    if (!successful) { return; }

    if (ui.checkBox_thresholded->isChecked())
    {
        bool successful = help::thresh(image);
        if (successful)
        {
            name = name + "_thresh";
        }
    }
    //apply analysis for depiction purposes
    bool analysisSuccessful = false;
    if (ui.checkBox_analysed->isChecked())
    {
        analysisSuccessful = m_analysis.analyseShape(image);
        if (analysisSuccessful)
        {
            name = name + "_roundAnalysed";
        }
    }
    else if (ui.checkBox_pcaAnalysis->isChecked())
    {
        vector<double> random;
        analysisSuccessful = m_analysis.analyseWithPCA(image, random);
        if (analysisSuccessful)
        {
            name = name + "_pcaAnalysed";
        }
    }
    else
    {
        if (!analysisSuccessful)
        {
            help::scaleData(image);
        }
    }

    string windowName = replacing ? "Image" : name;

    help::showWindow(image, scale, windowName);
}

void WindowMain::on_pushButton_showPlot_clicked()
{
    plotType type = static_cast<plotType>(ui.comboBox_showPlot->currentIndex());

    if (type == plotType::actFibers)
    {
        int cellNumber = ui.spinBox_showImage->value();

        if (cellNumber>=0 && cellNumber<=ui.spinBox_showImage->maximum()) 
        {
            plotOccurrenceInData(m_cellImages[cellNumber].actin_fibreAnglesPCA);
        }
    }
    else
    {
        vector<double> plottingData;
        for (int i = 0; i<m_cellImages.size();i++)
        {
              plottingData.push_back(m_cellImages[i].getQuantity(type));
        }
        plotOccurrenceInData(plottingData);
    }
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

    //TODO 8 bit jpg and 16 bit png
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

