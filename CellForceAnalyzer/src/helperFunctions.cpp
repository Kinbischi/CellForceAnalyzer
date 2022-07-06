#include "helperFunctions.h"

#include <string>
#include <algorithm>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <set>
#include <numeric>

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

    bool thresh(Mat& img, int thresholdValue)
    {
        if (img.channels()==3)
        {
            return false;
        }

        img = img.clone();
        //scale images before every thresholding
        scaleData(img);
        //convert to uint8 for threshold function 
        if (img.depth()==CV_16U)
        {
            img.convertTo(img, CV_8U, 1 / 256.0);
        }

        if (thresholdValue!=0)
        {
            cv::threshold(img, img, thresholdValue, 255, THRESH_BINARY); //used in fiber pca => no smoothing
        }
        else
        {
            // smoothen image
            int ksize = 5;
            cv::GaussianBlur(img, img, Size(ksize, ksize), 0, 0);
            auto thresholdValue = cv::threshold(img, img, 0, 255, THRESH_OTSU);
        }
        
        return true;
    }

    double average(vector<double> v)
    {
        return accumulate(v.begin(), v.end(), 0) / v.size();
    }

    double median(vector<double> v)
    {
        int n = v.size() / 2;
        nth_element(v.begin(), v.begin() + n, v.end());
        return v[n];
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

    void drawArrow(Mat& img, Point p, Point q, Scalar colour, const float scale)
    {
        double angle = atan2((double)p.y - q.y, (double)p.x - q.x); // angle in radians
        double hypotenuse = sqrt((double)(p.y - q.y) * (p.y - q.y) + (p.x - q.x) * (p.x - q.x));
        // Here we lengthen the arrow by a factor of scale
        q.x = (int)(p.x - scale * hypotenuse * cos(angle));
        q.y = (int)(p.y - scale * hypotenuse * sin(angle));
        line(img, p, q, colour, 1);
        // create the arrow hooks
        p.x = (int)(q.x + 9 * cos(angle + CV_PI / 4));
        p.y = (int)(q.y + 9 * sin(angle + CV_PI / 4));
        line(img, p, q, colour, 1);
        p.x = (int)(q.x + 9 * cos(angle - CV_PI / 4));
        p.y = (int)(q.y + 9 * sin(angle - CV_PI / 4));
        line(img, p, q, colour, 1);
    }

    void drawDoubleArrow(Mat& img, Point ctr, vector<Point2d> eigen_vecs, Scalar colour, double scalingLength)
    {
        Point p, q, temp;
        double arrowLength = scalingLength / 2.2;
        double hookLength = scalingLength / 3.3;
        double xDir = eigen_vecs[0].x;
        double yDir = eigen_vecs[0].y;
        double angle = atan2(yDir, xDir);// angle in radians

        // factor that is used to strech the arrow to arrowLength
        double strechFactor = sqrt(arrowLength * arrowLength / (xDir * xDir + yDir * yDir));
        p = Point(ctr.x - strechFactor * xDir, ctr.y - strechFactor * yDir);
        q = Point(ctr.x + strechFactor * xDir, ctr.y + strechFactor * yDir);

        // arrow body
        line(img, p, q, colour, 1);
        // first arrow hooks
        temp.x = (int)(q.x - hookLength * cos(angle + CV_PI / 4));
        temp.y = (int)(q.y - hookLength * sin(angle + CV_PI / 4));
        line(img, temp, q, colour, 1);
        temp.x = (int)(q.x - hookLength * cos(angle - CV_PI / 4));
        temp.y = (int)(q.y - hookLength * sin(angle - CV_PI / 4));
        line(img, temp, q, colour, 1);
        // second arrow hooks
        temp.x = (int)(p.x + hookLength * cos(angle + CV_PI / 4));
        temp.y = (int)(p.y + hookLength * sin(angle + CV_PI / 4));
        line(img, p, temp, colour, 1);
        temp.x = (int)(p.x + hookLength * cos(angle - CV_PI / 4));
        temp.y = (int)(p.y + hookLength * sin(angle - CV_PI / 4));
        line(img, p, temp, colour, 1);
    }

}