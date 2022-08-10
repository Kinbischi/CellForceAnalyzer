#pragma once

#include <opencv2/core.hpp>
#include <string>
#include <map>

// important enum
// In the showing function: None is used as a placeholder for the color image output (brightfield, vinculin and nucleus into RGB)
enum class channelType{brightfield, nucleus, vinculin, yap, actin, None};

class CustomImage
{
public:
	
	CustomImage(std::vector<cv::Mat>, std::vector<channelType>, std::string&);
	CustomImage()=default;
	CustomImage(const CustomImage&);

	void setChannel(cv::Mat, channelType);
	cv::Mat getChannel(channelType);

	void setName(std::string);
	std::string getName();
	std::string getName(channelType type);

	CustomImage cutImageOut(const cv::Rect&, std::string);
	cv::Mat createRGBimage();

	static const int max8bit = 255;
	static const int max16bit = 65535;


	std::string m_name;

	cv::Mat m_brightfieldChannel;
	cv::Mat m_nucleusChannel;
	cv::Mat m_vinculinChannel;
	cv::Mat m_yapChannel;
	cv::Mat m_actinChannel;

};



