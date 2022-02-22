#include "helperFunctions.h"

#include <string>
#include <algorithm>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <set>

using namespace std;
using namespace cv;

namespace help
{

    string& copiedDirectoryToNiceString(string& s)
    {
        replace(s.begin(), s.end(), '\\', '/');
        s.push_back('/');
        return s;
    }

    void showWindow(const Mat& image, double scale, const string windowName)
    {
        double aspect_ratio = double(image.cols) / image.rows;
        namedWindow(windowName, WINDOW_NORMAL);
        int height = int(image.rows * scale);
        int width = int(image.cols * scale);
        resizeWindow(windowName, width, height);
        imshow(windowName, image);
        pollKey();
    }

    bool thresh(Mat& img)
    {
        if (img.channels()==3)
        {
            return false;
        }
        //convert to uint8 for threshold function
        if (img.depth()==CV_16U)
        {
            img.convertTo(img, CV_8U, 1 / 256.0);
        }
        
        int ksize = 5;
        cv::GaussianBlur(img, img, Size(ksize, ksize), 0, 0);

        auto thresholdValue = cv::threshold(img, img, 0, 255, THRESH_OTSU);
        return true;
    }

}