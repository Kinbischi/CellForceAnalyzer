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

    void scaleData(Mat& image)
    {
        if (image.channels() == 3)
        {
            Mat bgr[3];
            vector<Mat> allChannels;
            split(image, bgr);
            if (image.depth() == CV_8U)
            {
                for (int i = 0; i < 3; i++)
                {
                    cv::normalize(bgr[i], bgr[i], 0, 255, NORM_MINMAX, CV_8U);
                    allChannels.push_back(bgr[i]);
                }
            }
            if (image.depth() == CV_16U)
            {
                for (int i = 0; i < 3; i++)
                {
                    cv::normalize(bgr[i], bgr[i], 0, 65535, NORM_MINMAX, CV_16U);
                    allChannels.push_back(bgr[i]);
                }
            }
            cv::merge(allChannels, image);
        }
        else
        {
            if (image.depth() == CV_8U)
            {
                cv::normalize(image, image, 0, 255, NORM_MINMAX, CV_8U);
            }
            if (image.depth() == CV_16U)
            {
                cv::normalize(image, image, 0, 65535, NORM_MINMAX, CV_16U);
            }
        }
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