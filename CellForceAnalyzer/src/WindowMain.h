#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_WindowMain.h"

#include "CustomImage.h"
#include "Preprocess.h"
#include "Cell.h"

class WindowMain : public QMainWindow
{
    Q_OBJECT

public:
    WindowMain(QWidget *parent = Q_NULLPTR);

private slots:
    void on_pushButton_channels_clicked();
    void on_pushButton_loadImages_clicked();
    void on_pushButton_showImage_clicked();
    void on_pushButton_writeOut_clicked();
    void imageNumberChanged(int);
    //void on_pushButton_test_clicked();
    void setCellTable(Cell cell);

private:
    Ui::WindowMainClass ui;

    std::string m_inpDir = "InputImages/";
    std::string m_outpDir = "OutputImages/";

    int m_imageNumber_show=0;

    //vector with order in which image channels are loaded
    std::vector<channelType> m_channels;

    std::vector<CustomImage> m_arrayImages;
    std::vector<Cell> m_cellImages;
    std::vector<cv::Mat> m_arrayImages_withYoloBoxes;

    Preprocess m_preprocess;

    
};

// TODO:
// at the moment cells are 16 bit but not normalized => normalize them?
// (before thresholding they get normalized)
// 
// 
//TODO: generally -> put const where it belongs to!

//TODO: perhaps nice to make customimg a virtual class and create an arrayImage type with yolo Mat in it??