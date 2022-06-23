#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_WindowMain.h"

#include "CustomImage.h"
#include "Preprocess.h"
#include "Cell.h"
#include "Analysis.h"

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
    void on_pushButton_showPlot_clicked();
    void on_pushButton_test_clicked();

    void imageNumberChanged(int);
    void radioButtonArraysChanged(bool);
    void radioButtonCellsChanged(bool);
    void radioButtonDeletedCellsChanged(bool);
    
    void setCellTable(Cell cell);



private:
    bool getImageToShow(cv::Mat&, std::string&, double&);

    void writeAnalysedDataToFile();


    Ui::WindowMainClass ui;

    std::string m_inpDir = "InputImages/"; // if nothing else specified in textbox => default is taken
    std::string m_outpDir = "OutputImages/"; // if nothing else specified in textbox => default is taken

    int m_imageNumber_show=0;

    //vector with order in which image channels are loaded
    std::vector<channelType> m_channels;

    std::vector<CustomImage> m_arrayImages;
    
    std::vector<cv::Mat> m_arrayImages_withYoloBoxes;

    std::vector<Cell> m_cellImages;
    std::vector<Cell> m_deletedCellImages;

    Cell m_averageAllCells;

    Preprocess m_preprocess;
    Analysis m_analysis;

    
};

// TODO:
// at the moment cells are 16 bit but not normalized => normalize them?
// (before thresholding they get normalized)
// 
// 
//TODO: generally -> put const where it belongs to!

//TODO: perhaps nice to make customimg a virtual class and create an arrayImage type with yolo Mat in it??

//TODO: absicherungs sachen
// mehrfaches laden => immer neu laden, nicht dazuzählen (eg arrays nd cells)

//TODO:
// fett wichtig => für sache wie density isch es wichtig unscaled images z becho vom thunder => wahrschiendlich 16 bit 
//=> scaling verzerrt die wahre wert => 
// generell zum abkläre => mehr staining => mehr signal? => then we would be fucked concerning e.g. densities

//Otsu thresholding is influenced depending on whether or not the image was normalized

//TODO: 255 oder max8bit??

// TODO array image yolo überlauf abwenden!
//TODO: neu laden alles vorher löschen?

//TODO: do not show boxes where borders were passed

//TODO: create experiment class => so that every time images with different niches are loaded => an experiment block is loaded

//TODO: why crash when previous plot not closed?
// why no titles possible in plots

//TODO: when should you clone? e.g. before every thresholding?
// 
// 
// 
//TODO check deleted cells from hannah data set 
//is everything working?
//Why tf are the images not stored in order in which they are loaded??
//Cells are kept where actin is completely detached from nucleus!!?? wtf get that out

//why is hannahs good looking cell twice detected but once deleted => in yolo procedure???

//TODO absturz abfangen bei show img falls 0 cells
