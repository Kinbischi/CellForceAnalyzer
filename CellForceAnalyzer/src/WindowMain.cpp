#include "WindowMain.h"

#include <string>
#include <set>
#include <cassert>
 
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
    newItem = new QTableWidgetItem(QString::fromStdString(to_string(100*cell.yapInNucleus)+" %"));
    ui.tableWidget_cell->setItem(3, 0, newItem);
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

    double averageCircularity=0, averageRoundness=0;
    for (auto& cell : m_cellImages)
    {
        m_preprocess.analyseNucleusShape(cell);
        m_preprocess.yapDistribution(cell);

        averageCircularity += cell.nucleus_circularity;
        averageRoundness += cell.nucleus_roundness;
    }
    
    averageCircularity = averageCircularity / m_cellImages.size();
    averageRoundness = averageRoundness / m_cellImages.size();

    // no delete is needed, as the QTableWidget takes ownership of the QTableWidgetItem and automized deletion
    QTableWidgetItem* newItem = new QTableWidgetItem(QString::fromStdString(to_string(m_arrayImages.size())));
    ui.tableWidget_all->setItem(0, 0, newItem);
    QTableWidgetItem* newItem2 = new QTableWidgetItem(QString::fromStdString(to_string(m_cellImages.size())));
    ui.tableWidget_all->setItem(1, 0, newItem2);
    QTableWidgetItem* newItem3 = new QTableWidgetItem(QString::fromStdString(to_string(averageCircularity)));
    ui.tableWidget_all->setItem(2, 0, newItem3);
    QTableWidgetItem* newItem4 = new QTableWidgetItem(QString::fromStdString(to_string(averageRoundness)));
    ui.tableWidget_all->setItem(3, 0, newItem4);

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

    if (m_arrayImages.empty())
    {
        QMessageBox::information(this, "Good Morning", "Read in Images!?");
        return;
    }

    CustomImage image;
    Mat outImg;
    string name;
    channelType channel = static_cast<channelType>(showInd);

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

    if (ui.checkBox_analysed->isChecked())
    {
        outImg = outImg.clone();
        bool successful = m_preprocess.analyseShape(outImg);
        if (successful)
        {
            name = name + "_roundAnalysed";
        }
    }
    else if (ui.checkBox_thresholded->isChecked())
    {
        outImg = outImg.clone();
        bool successful = help::thresh(outImg);
        if (successful)
        {
            name = name + "_thresh";
        }
    }
    

    string windowName = replacing ? "Image" : name;

    help::showWindow(outImg, scale, windowName);
}

void WindowMain::on_pushButton_writeOut_clicked()
{
    string outName = ui.lineEdit_outpDir->text().toStdString();
    if (!outName.empty())
    {
        m_outpDir = help::copiedDirectoryToNiceString(outName);
    }
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
