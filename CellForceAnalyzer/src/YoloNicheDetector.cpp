#include "YoloNicheDetector.h"

#include <string>
#include <fstream>
#include <vector>
#include <opencv2/dnn.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "helperFunctions.h"

using namespace std;
using namespace cv;
using namespace dnn;


// Initialize the parameters
float confThreshold;//= 0.5; // Confidence threshold => minimum confidence in order to be seen as an object
float nmsThreshold;//= 0.4;  // Non-maximum suppression threshold => area overlap for which bounding boxes from same class count as same object

int inpWidth = 416;        // Width of network's input image
int inpHeight = 416;       // Height of network's input image

// Use functions:
// Remove the bounding boxes with low confidence using non-maxima suppression
vector<Rect> postprocess(Mat& img, const vector<Mat>& out,vector<string> classes);
// Draw the predicted bounding box
void drawPred(int classId, float conf, int left, int top, int right, int bottom, Mat& img, vector<string> classes);
// Helping functions
vector<string> getOutputsNames(const Net& net);
vector<string> getClassNames(string yoloInput);


YoloNicheDetector::YoloNicheDetector()
{
    string modelConfiguration = m_yoloInputDir + "yolov3_custom.cfg"; //"yolov3.cfg";
    string modelWeights = m_yoloInputDir + "yolov3_custom_last_18_2.weights"; //"yolov3.weights";
    m_classes = getClassNames(m_yoloInputDir); //biz unnötig da nur 1 class
    // Load the network
    m_net = readNetFromDarknet(modelConfiguration, modelWeights);
    // CUDA
    m_net.setPreferableBackend(DNN_BACKEND_CUDA);
    m_net.setPreferableTarget(DNN_TARGET_CUDA);
}


vector<Cell> YoloNicheDetector::detectCells(CustomImage arrayImage, Mat& arrayImage_withYoloBoxes, float confThresh)
{
    confThreshold = confThresh;
    nmsThreshold = 0.4;
    
    Mat img = arrayImage.createRGBimage().clone();

    help::scaleData(img);

    if (img.depth()==CV_16U)
    {
        img.convertTo(img, CV_8U, 1.0 / 256.0);
    }

    Mat blob;
    Mat imgWithBoxes = img.clone();
    // Create a 4D blob from a img.

    blobFromImage(imgWithBoxes, blob, 1 / 255.0, cv::Size(inpWidth, inpHeight), Scalar(0, 0, 0), true, false);
    m_net.setInput(blob);
    vector<Mat> outs;
    m_net.forward(outs, getOutputsNames(m_net));

    // Remove the bounding boxes with low confidence
    vector<Rect> boxes = postprocess(imgWithBoxes, outs, m_classes);

    arrayImage_withYoloBoxes = imgWithBoxes;

    vector<Cell> cells;
    int i = 1;
    for (const auto& box : boxes)
    {
        string name = arrayImage.m_name + "_cell" + to_string(i);
        int rows = arrayImage.m_brightfieldChannel.rows;
        int cols = arrayImage.m_brightfieldChannel.cols;

        if (box.x > 0 && box.y > 0 && box.x + box.width < cols && box.y + box.height < rows)
        {
            CustomImage c = arrayImage.cutImageOut(box,name);
            Cell cell(c);
            cells.push_back(cell);
            i++;
        }
    }
    return cells;
}


// Remove the bounding boxes with low confidence using non-maxima suppression
vector<Rect> postprocess(Mat& frame, const vector<Mat>& outs, vector<string> classes)
{
    vector<int> classIds;
    vector<float> confidences;
    vector<Rect> boxes;

    for (size_t i = 0; i < outs.size(); ++i)
    {
        // Scan through all the bounding boxes output from the network and keep only the
        // ones with high confidence scores. Assign the box's class label as the class
        // with the highest score for the box.
        float* data = (float*)outs[i].data;
        for (int j = 0; j < outs[i].rows; ++j, data += outs[i].cols)
        {
            Mat scores = outs[i].row(j).colRange(5, outs[i].cols);
            Point classIdPoint;
            double confidence;
            // Get the value and location of the maximum score
            minMaxLoc(scores, 0, &confidence, 0, &classIdPoint);
            if (confidence > confThreshold)
            {
                int centerX = (int)(data[0] * frame.cols);
                int centerY = (int)(data[1] * frame.rows);
                int width = (int)(data[2] * frame.cols);
                int height = (int)(data[3] * frame.rows);
                int left = centerX - width / 2;
                int top = centerY - height / 2;

                classIds.push_back(classIdPoint.x);
                confidences.push_back((float)confidence);
                boxes.push_back(Rect(left, top, width, height));
            }
        }
    }

    // Perform non maximum suppression to eliminate redundant overlapping boxes with lower confidences
    vector<int> indices;
    vector<Rect> finalBoxes;
    NMSBoxes(boxes, confidences, confThreshold, nmsThreshold, indices);

    for (size_t i = 0; i < indices.size(); ++i)
    {
        int idx = indices[i];
        Rect box = boxes[idx];
        finalBoxes.push_back(box);

        drawPred(classIds[idx], confidences[idx], box.x, box.y,
            box.x + box.width, box.y + box.height, frame, classes);
    }
    return finalBoxes;
}

// Draw the predicted bounding box
void drawPred(int classId, float conf, int left, int top, int right, int bottom, Mat& frame, vector<string> classes)
{
    //Draw a rectangle displaying the bounding box
    rectangle(frame, Point(left, top), Point(right, bottom), Scalar(255, 178, 50), 3);

    //Get the label for the class name and its confidence
    string label = format("%.2f", conf);
    if (!classes.empty())
    {
        CV_Assert(classId < (int)classes.size());
        label = classes[classId] + ":" + label;
    }

    //Display the label at the top of the bounding box
    int baseLine;
    Size labelSize = getTextSize(label, FONT_HERSHEY_SIMPLEX, 1.25, 2, &baseLine);
    top = max(top, labelSize.height);

    rectangle(frame, Point(left, top - round(1.5 * labelSize.height)), Point(left + round(1.2 * labelSize.width), top + baseLine), Scalar(255, 255, 255), FILLED);
    putText(frame, label, Point(left, top), FONT_HERSHEY_SIMPLEX, 1.5, Scalar(0, 0, 0), 2);
}

// Get the names of the output layers
vector<string> getOutputsNames(const Net& net)
{
    static vector<string> names;
    if (names.empty())
    {
        //Get the indices of the output layers, i.e. the layers with unconnected outputs
        vector<int> outLayers = net.getUnconnectedOutLayers();

        //get the names of all the layers in the network
        vector<string> layersNames = net.getLayerNames();

        // Get the names of the output layers in names
        names.resize(outLayers.size());
        for (size_t i = 0; i < outLayers.size(); ++i)
            names[i] = layersNames[outLayers[i] - 1];
    }
    return names;
}

vector<string> getClassNames(string yoloInput)
{
    vector<string> classes;
    string classesFile = yoloInput + "obj.names"; // "coco.names";
    ifstream ifs(classesFile.c_str());
    string line;

    while (getline(ifs, line))
    {
        classes.push_back(line);
    }
    return classes;
}