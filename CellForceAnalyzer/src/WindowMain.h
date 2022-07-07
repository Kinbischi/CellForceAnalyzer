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
    void loadNiceCellImages();

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


//TODO: generally -> put const where it belongs to!

//TODO: perhaps nice to make customimg a virtual class and create an arrayImage type with yolo Mat in it??
// 
//TODO: create experiment class => so that every time images with different niches are loaded => an experiment block is loaded
//TODO: multi data set => mehrfaches laden? => immer neu laden/ dazuzählen

// generell zum abkläre => mehr staining => mehr signal? => then we would be fucked concerning e.g. densities

//TODO: 255 oder max8bit??

//TODO: why crash when previous plot not closed?
// why no titles possible in plots

//TODO absturz abfangen bei show img falls 0 cells

//TODO: yolo => was trained on good rgb images (jpg) => what happens if scaleData does not work well??
// e.g. one pixel is very high intensity in an 16 bit image and scaling fucks up => would at least be seeable in show image
// but still keep it in mind => yolo can fail if one pixel is high and scaling fails => better scaling than min max??

//TODO: try fourrier, analysis fibre counting, dot counting for focal adhesions

//TODO: bug => first actin switching to array yolo => max image set to 0...

//TODO: check if 16 bit images works? e.g. pca analysis or show image

//TODO: for fun => if one thresholding/ analysis is checked => uncheck the other

//TODO: load in squareLength and minRatio from gui also for analysis

//TODO: try out intensity wise PCA => every point given to PCA as often as its intensity
// eventuell mehrere square lengths weisen den arrow am ort zu?

// try storing the ev ratios and iterate the thresholding and take the largest ev ratio direction
// TODO: try out intensity again => did you pass a rgb image?