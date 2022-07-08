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

void WindowMain::loadNiceCellImages()
{
    vector<string> imageNames;
    vector<vector<Rect>> rectangles;
    vector<vector<channelType>> channels;

    //Oksana's cell
    vector<channelType> channelsOrder1;
    channelsOrder1.push_back(channelType::nucleus);
    channelsOrder1.push_back(channelType::actin);
    channelsOrder1.push_back(channelType::brightfield);

    imageNames.push_back("50x50_niches_RGD_3mM_2022_02_22__14_11_14_Maximum intensity projection_Filter.tif");
    vector<Rect> rVec;
    rVec.push_back(Rect(220, 400, 300, 300));
    rectangles.push_back(rVec);
    rVec.clear();
    channels.push_back(channelsOrder1);

    //Leslie's cells
    vector<channelType> channelsOrder2;
    channelsOrder2.push_back(channelType::nucleus);
    channelsOrder2.push_back(channelType::brightfield);
    channelsOrder2.push_back(channelType::actin);

    imageNames.push_back("300321_PEGvsGELMAlif.lif - Image006-1 (RGB).tif");
    rVec.push_back(Rect(1000, 1180, 500, 500));
    rVec.push_back(Rect(140, 1180, 340, 340));
    rVec.push_back(Rect(1250, 200, 750, 250));
    rVec.push_back(Rect(880, 1630, 350, 250));
    rectangles.push_back(rVec);
    rVec.clear();
    channels.push_back(channelsOrder2);

    imageNames.push_back("300321_PEGvsGELMAlif.lif - Image001-1 (RGB).tif");
    
    rVec.push_back(Rect(1020, 420, 300, 250));
    rVec.push_back(Rect(10, 1240, 270, 530));
    rVec.push_back(Rect(840, 1080, 270, 260));
    rVec.push_back(Rect(1640, 1130, 310, 240));
    rVec.push_back(Rect(1280, 1090, 220, 310));
    rVec.push_back(Rect(200, 1150, 330, 210));
    rectangles.push_back(rVec);
    channels.push_back(channelsOrder2);

    for (int i = 0; i<imageNames.size(); i++)
    {
        Mat img = imread(m_inpDir + imageNames[i], IMREAD_UNCHANGED);

        Mat bgr[3];
        split(img, bgr);
        vector<Mat> matImages;
        matImages.push_back(bgr[0]);
        matImages.push_back(bgr[1]);
        matImages.push_back(bgr[2]);
        CustomImage image(matImages, channels[i], imageNames[i]);
        m_arrayImages.push_back(image);

        for (int j = 0; j<rectangles[i].size(); j++)
        {
            string name_cell = imageNames[i] + "_cell" + to_string(j);
            CustomImage ce = image.cutImageOut(rectangles[i][j], name_cell);
            Cell cell(ce);
            m_cellImages.push_back(cell);
        }
    }

}

void WindowMain::on_pushButton_test_clicked()
{
    
    Mat img = imread(m_inpDir + "300321_PEGvsGELMAlif.lif - Image001-1 (RGB).tif", IMREAD_UNCHANGED);
    Rect r(1020, 420, 300, 250);
    rectangle(img,r, Scalar(128, 255, 0));
    Rect r1(10, 1240, 270, 500);
    rectangle(img, r1, Scalar(128, 255, 0));
    Rect r2(840, 1080, 270, 260);
    rectangle(img, r2, Scalar(128, 255, 0));
    Rect r3(1640, 1130, 310, 240);
    rectangle(img, r3, Scalar(128, 255, 0));
    Rect r4(1280, 1090, 220, 310);
    rectangle(img, r4, Scalar(128, 255, 0));
    Rect r5(200, 1150, 330, 210);
    rectangle(img, r5, Scalar(128, 255, 0));
    help::showWindow(img,0.5);
   
    loadNiceCellImages();
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

    m_preprocess.loadImages(m_arrayImages, m_inpDir, m_channels);

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
            m_analysis.analyseNucleus(m_cellImages[i]); //TODO: what analysis did not work => include e.g. the assert(isContour),...
            m_analysis.analyseYap(m_cellImages[i]);
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

    m_averageAllCells = m_preprocess.getAverageProperties(m_cellImages);
    
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
    double scale = ui.doubleSpinBox_scale->value();
    bool replacing = ui.checkBox_replace->isChecked();
    Mat image;
    string name;

    int pca_squareLength = ui.spinBox_squareLengthPCA->value();
    double pca_minEigvalRatio = ui.doubleSpinBox_minEigValRatioPCA->value();

    //load image
    bool successful = getImageToShow(image, name, scale);
    if (!successful) { return; }

    //thresholdings
    if (ui.checkBox_thresholdManually->isChecked())
    {
        int threshValue = ui.spinBox_thresholdManually->value();
        successful = help::thresh(image,threshValue);
        if (successful) { name = name + "_threshManually"+to_string(threshValue); }
    }
    else if (ui.checkBox_thresholdOtsu->isChecked())
    {
        successful = help::thresh(image);
        if (successful) { name = name + "_threshOtsu"; }
    }
    else if (ui.checkBox_thresholdPCAoptSingle->isChecked())
    {
        int optimalThresholding = m_analysis.getOptimalThresholdingForPCA(image, pca_squareLength, pca_minEigvalRatio);
        successful = help::thresh(image, optimalThresholding);
        if (successful) { name = name + "_threshPCA"; }
    }
    else if (ui.checkBox_thresholdPCAoptSquares->isChecked())
    {
        Mat imgThresh;
        help::pcaType type = help::pcaType::squarePCAoptimizedThresh;
        m_analysis.getThresholdedImage(image, imgThresh, type, pca_squareLength, pca_minEigvalRatio, 0);
    }

    if(ui.checkBox_edgeDetection->isChecked())
    {
        cv::blur(image, image, Size(3, 3));
        int kernel_size = 3;
        int lowThreshold = 50;
        cv::Canny(image, image, lowThreshold, lowThreshold * 3, kernel_size); //TODO: what is wrong here?
    }
    if (ui.checkBox_FAdetection->isChecked())
    {
        GaussianBlur(image, image, Size(9, 9), 2, 2);
        vector<Vec3f> circles;
        HoughCircles(image, circles, HOUGH_GRADIENT, 2, image.rows / 4, 200, 100);
        for (size_t i = 0; i < circles.size(); i++)
        {
            Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
            int radius = cvRound(circles[i][2]);
            
            cvtColor(image,image, COLOR_GRAY2RGB);
            // draw the circle center
            circle(image, center, 3, Scalar(0, 255, 0), -1, 8, 0);
            // draw the circle outline
            circle(image, center, radius, Scalar(0, 0, 255), 3, 8, 0);
        }
    }
    

    //apply analysis for depiction purposes
    bool analysisSuccessful = false;
    if (ui.checkBox_analysed->isChecked())
    {
        analysisSuccessful = m_analysis.analyseShape(image);
        if (analysisSuccessful) { name = name + "_roundAnalysed"; }
    }
    else if (ui.checkBox_pcaAnalysis->isChecked())
    {
        vector<double> random;
        analysisSuccessful = m_analysis.analyseWithPCA(image, random, pca_squareLength, pca_minEigvalRatio);
        if (analysisSuccessful) { name = name + "_pcaAnalysed"; }
    }

    else //TODO: check what happens for unscaled 16 bit image that were analyzed => probs not scaled ;(
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

