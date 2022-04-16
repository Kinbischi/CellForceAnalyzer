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

WindowMain::WindowMain(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);

    connect(ui.spinBox_showImage, SIGNAL(valueChanged(int)), this, SLOT(imageNumberChanged(int)));

    //it grabs the channels (order of input channels) as default values from the .ui file
    WindowChannels channelsWindow;
    m_channels = channelsWindow.set_ChannelTypes();
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

/*
void WindowMain::on_pushButton_test_clicked()
{
    CustomImage testImage = m_cellImages_thresholded[0];

    Mat outImg = m_preprocess.analyseRoundShape(testImage.m_nucleusChannel);
    help::showWindow(outImg,4);
}
*/



void WindowMain::imageNumberChanged(int index)
{
    m_imageNumber_show = index;
    setCellTable(m_cellImages[index]);
}

void WindowMain::setCellTable(Cell cell)
{
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


void WindowMain::on_pushButton_loadImages_clicked()
{
    string inName = ui.lineEdit_inpDir->text().toStdString();
    if (!inName.empty())
    {
        m_inpDir = help::copiedDirectoryToNiceString(inName);
    }

    float confThreshold = ui.doubleSpinBox_confThresh->value();
    float nmsThreshold = ui.doubleSpinBox_nonMaxThresh->value();

    m_preprocess.loadImages(m_arrayImages, m_inpDir, m_channels);

    m_preprocess.applyYolo(m_arrayImages, m_cellImages, m_arrayImages_withYoloBoxes,confThreshold,nmsThreshold);

    double summedNucleusArea=0, summedNucleusCircularity=0, summedNucleusRoundness=0, summedActinArea=0, summedActinDensity=0,
        summedActinMaxLength=0, summedYapInNucleus=0;
    vector<double> actinPCAangles;
    for (auto& cell : m_cellImages)
    {
        m_analysis.analyseNucleusShape(cell);
        m_analysis.analyseYapInNucleus(cell);
        m_analysis.analyseActin(cell);

        summedNucleusArea += cell.nucleus_area;
        summedNucleusCircularity += cell.nucleus_circularity;
        summedNucleusRoundness += cell.nucleus_roundness;

        summedActinArea += cell.actin_area;
        summedActinDensity += cell.actin_density;
        summedActinMaxLength += cell.actin_maxLength;
        actinPCAangles.push_back(cell.actin_PCAangle);

        summedYapInNucleus += cell.yap_inNucleus;
    }

    //todo: this funct before?
    //load bad cells count into gui?
    m_analysis.removeBadCells(m_cellImages);
    
    if (m_cellImages.size()>0)
    {
        m_averageAllCells.nucleus_area = summedNucleusArea / m_cellImages.size();
        m_averageAllCells.nucleus_circularity = summedNucleusCircularity / m_cellImages.size();
        m_averageAllCells.nucleus_roundness = summedNucleusRoundness / m_cellImages.size();
        m_averageAllCells.actin_area = summedActinArea / m_cellImages.size();
        m_averageAllCells.actin_density = summedActinDensity / m_cellImages.size();
        m_averageAllCells.actin_maxLength = summedActinMaxLength / m_cellImages.size();
        m_averageAllCells.actin_PCAangle = static_cast<int>(m_analysis.findPointWithMostNeighbours(actinPCAangles, 10));
        m_averageAllCells.yap_inNucleus = summedYapInNucleus / m_cellImages.size();

        // no delete is needed, as the QTableWidget takes ownership of the QTableWidgetItem and automized deletion
        QTableWidgetItem* newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_arrayImages.size())));
        ui.tableWidget_all->setItem(0, 0, newItem);
        QTableWidgetItem* newItem2 = new QTableWidgetItem(QString::fromStdString(to_string(m_cellImages.size())));
        ui.tableWidget_all->setItem(1, 0, newItem2);
        QTableWidgetItem* newItem3 = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.nucleus_circularity)));
        ui.tableWidget_all->setItem(2, 0, newItem3);
        QTableWidgetItem* newItem4 = new QTableWidgetItem(QString::fromStdString(to_string(m_averageAllCells.nucleus_roundness)));
        ui.tableWidget_all->setItem(3, 0, newItem4);
    }
}

void WindowMain::on_pushButton_showImage_clicked()
{
    if (m_arrayImages.empty())
    {
        QMessageBox::information(this, "Good Morning", "Read in Images!?");
        return;
    }

    int imageNumber = ui.spinBox_showImage->value();
    int showInd = ui.comboBox_showImage->currentIndex();
    double scale = ui.doubleSpinBox_scale->value();
    bool replacing = ui.checkBox_replace->isChecked();

    CustomImage image;
    Mat outImg;
    string name;
    channelType channel = static_cast<channelType>(showInd);

    if (outImg.empty())
    {
        QMessageBox::information(this, "Good Morning", "This channel was not loaded!");
        return;
    }

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
        return;
    }

    outImg = outImg.clone();

    if (ui.checkBox_analysed->isChecked())
    {
        bool successful = m_analysis.analyseShape(outImg);
        if (successful)
        {
            name = name + "_roundAnalysed";
        }
    }
    else if (ui.checkBox_thresholded->isChecked())
    {
        bool successful = help::thresh(outImg);
        if (successful)
        {
            name = name + "_thresh";
        }
    }

    help::scaleData(outImg);

    string windowName = replacing ? "Image" : name;

    help::showWindow(outImg, scale, windowName);
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
            imwrite(m_outpDir + entry.m_name + ".jpg", entry.createRGBimage());
        } 
    }
    if (ui.checkBox_singleCells->isChecked())
    {
        for (auto entry : m_cellImages)
        {
            imwrite(m_outpDir + entry.m_name + ".jpg", entry.createRGBimage());
        }
    }
    if (ui.checkBox_cellArraysWithBoxes->isChecked())
    {
        for (int i=0;i< m_arrayImages_withYoloBoxes.size();i++)
        {
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

