#pragma once
#include <string>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <set>
#include <map>

namespace help
{
    const double M_PI = 3.14159265358979323846;

    std::string& copiedDirectoryToNiceString(std::string&);
    bool thresh(cv::Mat&, int = 0);
    void scaleData(cv::Mat&);
    double average(std::vector<double>);
    double median(std::vector<double>);
    
    
    // functions related to images
    void showWindow(const cv::Mat&, double = 1, const std::string = "Image");
    void drawArrow(cv::Mat&, cv::Point, cv::Point, cv::Scalar, const float scale = 0.2);
    void drawDoubleArrow(cv::Mat&, cv::Point, std::vector<cv::Point2d>, cv::Scalar, double);

    // the following are template functions that go into the header

    // function to determine what pixel values are present in an image
    template <typename T>
    std::set<T> elements(cv::Mat im)
    {
        std::set<T> s;
        cv::MatIterator_<T> it, end;
        for (it = im.begin<T>(), end = im.end<T>(); it != end; ++it)
        {
            T a = *it;
            s.insert(a);
        }
        return s;
    }

    // function to determine how often a pixel value is present in an image
    // does not work with 16 bit 3D images (because of Vec3b)
    template <typename T>
    std::vector<std::map<T, int>> elementsMap(cv::Mat img)
    {
        std::map<T, int> itemsCountR, itemsCountG, itemsCountB;
        std::vector<std::map<T, int>> allChannels;
        if (img.depth() != CV_8U)
        {
            return allChannels;
        }

        if (img.channels() == 1)
        {
            cv::MatIterator_<T> it, end;
            for (it = img.begin<T>(), end = img.end<T>(); it != end; ++it)
            {
                T a = *it;
                ++itemsCountB[a];
            }
        }
        if (img.channels() == 3)
        {
            for (int i = 0; i < img.rows; i++)
                for (int j = 0; j < img.cols; j++)
                {
                    auto pix = img.at<cv::Vec3b>(i, j);
                    T a = img.at<cv::Vec3b>(i, j)[0];
                    T b = img.at<cv::Vec3b>(i, j)[1];
                    T c = img.at<cv::Vec3b>(i, j)[2];
                    ++itemsCountB[a];
                    ++itemsCountG[b];
                    ++itemsCountR[c];
                }
        }
        if (!itemsCountB.empty())
        {
            allChannels.push_back(itemsCountB);
        }
        if (!itemsCountG.empty())
        {
            allChannels.push_back(itemsCountG);
        }
        if (!itemsCountR.empty())
        {
            allChannels.push_back(itemsCountR);
        }

        return allChannels;
    }
}