#pragma once

#include <opencv2/core.hpp>
#include <string>
#include <map>

// important enum

// In the showing function: None is used as a placeholder for the color image output (brightfield, cytoplasm and nucleus into RGB)
enum class channelType{brightfield, nucleus, cytoplasm, yap, None};

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

	void CustomImage::threshold();
	CustomImage getThresholdedImage();

	static const int max8bit = 255;
	static const int max16bit = 65535; //change so that every class has its own? (either 255 or 65535)

	std::string m_name;

	//all stored as normalized (highest pixel 65535 and lowest pixel 0) for 16 bit images (ushort)
	cv::Mat m_brightfieldChannel;
	cv::Mat m_nucleusChannel;
	cv::Mat m_cytoplasmChannel;
	cv::Mat m_yapChannel;

};



