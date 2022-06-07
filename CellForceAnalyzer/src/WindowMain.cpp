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
#include "matplotlibcpp.h"

using namespace std;
using namespace cv;
namespace plt = matplotlibcpp;

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
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(cell.actin_PCAangle)));
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


void WindowMain::on_pushButton_test_clicked()
{
    std::vector<double> y = { 1, 3, 2, 4 };
    std::vector<double> x = { 1, 2, 3, 4 };
    plt::plot(x, y);
    plt::show();
    /*
    //ui.tableWidget_cell->clearContents();
    Mat img = imread(m_inpDir+"50x50_niches_RGD_3mM_2022_02_22__14_11_14_Maximum intensity projection_Filter.tif", IMREAD_UNCHANGED);
    
    Rect r(220,400,300,300);
    //rectangle(img,r, Scalar(255, 178, 50));

    std::vector<channelType> channelsOrder;
    channelsOrder.push_back(channelType::nucleus);
    channelsOrder.push_back(channelType::actin);
    channelsOrder.push_back(channelType::brightfield);
    
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

    Mat actin = m_cellImages[0].getChannel(channelType::actin);
    int squaresPerSideX = 7;
    int squaresPerSideY = 7;

    //investigate bug for 20! for arrows depiction
    int lengthX = actin.rows/ squaresPerSideX;
    int lengthY = actin.cols/ squaresPerSideY;

    vector<Mat> subImages;

    Mat actinThresh = actin.clone();
    help::thresh(actinThresh);
    Mat actinThreshPic = actinThresh.clone();

    help::showWindow(actin, 2, "actin");
    waitKey(0);
    help::showWindow(actinThresh, 2, "actin");

    cv::cvtColor(actin, actin, COLOR_GRAY2RGB);
    cv::cvtColor(actinThreshPic, actinThreshPic, COLOR_GRAY2RGB);

    vector<vector<vector<double>>> arrows;

    for (int i = 0; i < squaresPerSideX; i++)
    {
        for (int j = 0; j < squaresPerSideY; j++)
        {
            Rect rect(i*lengthX, j*lengthY, lengthX, lengthY);
            Mat subImg = actinThresh(rect).clone();

            //rectangle(actinThresh, rect, Scalar(255, 178, 50));
            //help::showWindow(actinThresh,1,"actin");

            vector<Point> points = m_analysis.getWhitePointsFromThresholdedImage(subImg);
            // TODO: only get points from largest feature??

            
            //help::showWindow(subImg);
            //for (auto point : points)
            //{
            //    subImg.at<uchar>(point) = 128;
            //}
            //help::showWindow(subImg);
            

            Point center;
            vector<Point2d> eigen_vecs(2);
            vector<double> eigen_val(2);
            bool worked = m_analysis.getPCAorientation2(points, eigen_vecs, eigen_val, center);

            if (worked)// TODO change to worked! (bool)
            {
                waitKey(0);
                cv::cvtColor(subImg, subImg, COLOR_GRAY2RGB);

                m_analysis.drawAxis2(subImg, Point(lengthX/2,lengthY/2), eigen_vecs, eigen_val, Scalar(0, 128, 255), 2, 1);
                help::showWindow(subImg, 5);
                center = Point((i + 0.5) * lengthX, (j + 0.5) * lengthY);
                //m_analysis.drawAxis(actin, Point((i+0.5) * lengthX, (j+0.5) * lengthY), Point((i + 0.5) * lengthX + arrow_x, (j + 0.5) * lengthY + arrow_y), Scalar(0, 128, 255), 1);
                m_analysis.drawAxis2(actinThreshPic, center, eigen_vecs, eigen_val, Scalar(0, 128, 255),2,1);
                
                Mat actinCopy = actin.clone();
                rectangle(actinThreshPic, rect, Scalar(0, 255, 128));
                help::showWindow(actinThreshPic, 2,"actin");
                
            }
        }
    }
    */

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
    float nmsThreshold = ui.doubleSpinBox_nonMaxThresh->value();

    m_preprocess.loadImages(m_arrayImages, m_inpDir, m_channels);

    m_preprocess.applyYolo(m_arrayImages, m_cellImages, m_arrayImages_withYoloBoxes,confThreshold,nmsThreshold);

    
    bool isDead, analysisFailed;
    vector<int> cellsToDelete;
    for (int i=0; i<m_cellImages.size(); i++)
    {
        try 
        {
            isDead = m_analysis.isDeadCell(m_cellImages[i]);
            if (isDead) 
            {
                cellsToDelete.push_back(i);
                m_deletedCellImages.push_back(m_cellImages[i]);    
                continue; 
            }
            m_analysis.analyseNucleusShape(m_cellImages[i]);
            m_analysis.analyseYapInNucleus(m_cellImages[i]);
            m_analysis.analyseActin(m_cellImages[i]);
        }
        catch (...)
        {
            analysisFailed = true;

            cellsToDelete.push_back(i);
            m_deletedCellImages.push_back(m_cellImages[i]);
            QMessageBox::information(this, "Hmmm", "Not all analysis was able to be conducted");
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

void WindowMain::on_pushButton_showImage_clicked()
{
    double scale = ui.doubleSpinBox_scale->value();
    bool replacing = ui.checkBox_replace->isChecked();
    Mat image;
    string name;

    //load image
    bool successful = getImageToShow(image, name, scale);
    if (!successful) { return; }

    //apply analysis for depiction purposes
    if (ui.checkBox_analysed->isChecked())
    {
        bool successful = m_analysis.analyseShape(image);
        if (successful)
        {
            name = name + "_roundAnalysed";
        }
    }
    else if (ui.checkBox_thresholded->isChecked())
    {
        bool successful = help::thresh(image);
        if (successful)
        {
            name = name + "_thresh";
        }
    }

    help::scaleData(image);

    string windowName = replacing ? "Image" : name;

    help::showWindow(image, scale, windowName);
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
            << cell.actin_maxLength << "," << cell.actin_PCAangle << "," << cell .yap_inNucleus << "\n";
    }
    
    myfile <<"Averages," << m_averageAllCells.nucleus_area << "," << m_averageAllCells.nucleus_circularity << "," 
        << m_averageAllCells.nucleus_roundness << "," << m_averageAllCells.actin_area << "," << m_averageAllCells.actin_density << "," 
        << m_averageAllCells.actin_maxLength << "," << m_averageAllCells.actin_PCAangle << "," << m_averageAllCells.yap_inNucleus;

}

