#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_WindowMain.h"

#include "CustomImage.h"
#include "Preprocess.h"
#include "Cell.h"
#include "Analysis.h"
#include "ParametersUI.h"
#include "Plotting.h"
#include "AnalysisFiberDirection.h"
#include "Display.h"
#include "dataContainer.h"


class WindowMain : public QMainWindow
{
    Q_OBJECT

public:
    WindowMain(QWidget *parent = Q_NULLPTR);

private slots:
    void on_pushButton_channels_clicked();
    void on_pushButton_loadImages_clicked();
    void on_pushButton_applyYolo_clicked();
    void on_pushButton_showImage_clicked();
    void on_pushButton_writeOut_clicked();
    void on_pushButton_showPlot_clicked();
    void on_pushButton_test_clicked();
    void on_pushButton_conductAnalysisOnSingleCell_clicked();
    void on_pushButton_conductAnalysisOnAllCells_clicked();

    void imageNumberChanged(int);
    void radioButtonArraysChanged(bool);
    void radioButtonCellsChanged(bool);
    void radioButtonRemovedCellsChanged(bool);


private:
    void updateCellAnalysisTable(Cell cell);
    void updateGeneralAnalysisTable();
    void updateGeneralTable();

    void writeAnalysedDataToFile();
    void updateParameters();
    void getImageWorked(int);

    Ui::WindowMainClass ui;

    std::string m_inpDir = "InputImages/"; // if nothing else specified in textbox => default is taken
    std::string m_outpDir = "OutputImages/"; // if nothing else specified in textbox => default is taken

    int m_imageNumber_show=0;

    dataContainer m_data;

    Preprocess m_preprocess;
    ParametersUI m_params; //Parameters from UI
    Analysis m_analysis;
    Plotting m_plotting;
    AnalysisFiberDirection m_analysisFiberDir;
    Display m_display;
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

//TODO absturz abfangen bei show img falls 0 cells oder 0 yoloBoxes oder 0 deleted cells => perhaps try bei get image?

//TODO: yolo => was trained on good rgb images (jpg) => what happens if scaleData does not work well??
// e.g. one pixel is very high intensity in an 16 bit image and scaling fucks up => would at least be seeable in show image
// but still keep it in mind => yolo can fail if one pixel is high and scaling fails => better scaling than min max??

//TODO: try fourrierb, analysis fibre counting, dot counting for focal adhesions

//TODO: bug => first actin switching to array yolo => max image set to 0...

//TODO: check if 16 bit images works? e.g. pca analysis or show image

//TODO separate loading and analysis application

//TODO: load in all variables from  into constants in helperfile? not always in function

//Idee summieren entlang kanten => schauen ob pixel uniformely dist or not.
//PCA is not working if two fibers are at left and right edge => improvable?

//TODO create plotting class and showing class

//Cells empty parameter,.... => no if(cells.size())

//TODO pca out of analysis => fiberAnalysis.h in analysis and windowmain

//TODO for show plot => try and catch
// separate cells dead analysis, analysis, yolo,... => always: with checkbox for excluding cells that were bad (dead, analysisfailed)

//todo => handle besser wann analysen erlaubt => nicht menr if ().size()>0 usw => bool analysisConducted in parameters zb, bool imagesLoaded
//=> zb vor plotting funktion abfragen

//TODO: remove bad cells implementation
//check if image loading works as desired
// better featback if eg channel not loaded or image type not available

//TODO: try putting islands from optimally thresholded image into pca => not squares => then analysis gets more indep of squares

//Leslie würde gern wissen: => unterschied in yap localization zwischen 5%,10% and coverslip

// try blur detection of nuclei => single line of code: cv2.Laplacian(image, cv2.CV_64F).var()

//TODO update cell number after application of YOLO