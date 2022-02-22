#include "CustomImage.h"

#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

#include "helperFunctions.h"

using namespace cv;
using namespace std;


CustomImage::CustomImage(const CustomImage& other)
{
	m_name = other.m_name;
	m_brightfieldChannel = other.m_brightfieldChannel.clone();
	m_nucleusChannel = other.m_nucleusChannel.clone();
	m_cytoplasmChannel = other.m_cytoplasmChannel.clone();
	m_yapChannel = other.m_yapChannel.clone();
}

void CustomImage::setName(std::string str)
{
	if (str.find("_ch") != std::string::npos)
	{
		m_name = str.erase(str.find("_ch"));
	}
	if (str.find(".tif") != std::string::npos)
	{
		m_name = str.erase(str.find(".tif"));
	}
	else
	{
		m_name = str;
	}
}
string CustomImage::getName()
{
	return m_name;
}

string CustomImage::getName(channelType type)
{
	string name=m_name;
	switch (type)
	{
	case channelType::brightfield:
		name += "_bf";
		break;
	case channelType::nucleus:
		name += "_nc";
		break;
	case channelType::cytoplasm:
		name += "_cy";
		break;
	case channelType::yap:
		name += "_yap";
		break;
	case channelType::None:
		name += "_rgb";
		break;
	}
	return name;
}

void CustomImage::setChannel(Mat img, channelType type)
{
	//scale image so that highest pixel is white and lowest is black
	if (img.depth() == CV_8U)
	{
		cv::normalize(img, img, 0, max8bit, NORM_MINMAX, CV_8U);
	}
	if (img.depth() == CV_16U)
	{
		cv::normalize(img, img, 0, max16bit, NORM_MINMAX, CV_16U);
	}
	
	switch (type)
	{
	case channelType::brightfield:
		m_brightfieldChannel = img;
		break;
	case channelType::nucleus:
		m_nucleusChannel = img;
		break;
	case channelType::cytoplasm:
		m_cytoplasmChannel = img;
		break;
	case channelType::yap:
		m_yapChannel = img;
		break;
	}
}

// For None => image returns the RGB image
Mat CustomImage::getChannel(channelType type)
{
	Mat img;
	switch (type)
	{
	case channelType::brightfield:
		img=m_brightfieldChannel;
		break;
	case channelType::nucleus:
		img = m_nucleusChannel;
		break;
	case channelType::cytoplasm:
		img = m_cytoplasmChannel;
		break;
	case channelType::yap:
		img = m_yapChannel;
		break;
	case channelType::None:
		img = this->createRGBimage();
		break;
	}
	return img;
}

CustomImage::CustomImage(vector<Mat> images, vector<channelType> channels,string& imgName)
{
	setName(imgName);
	for (int i = 0; i < images.size(); i++)
	{
		setChannel(images[i], channels[i]);
	}
}

void CustomImage::threshold()
{
	m_name = m_name + "_thresholded";
	help::thresh(m_brightfieldChannel);
	help::thresh(m_nucleusChannel);
	help::thresh(m_cytoplasmChannel);
	help::thresh(m_yapChannel);
}

// Not used at the moment  => delete later??
CustomImage CustomImage::getThresholdedImage()
{
	CustomImage temp(*this);
	temp.threshold();
	return temp;
}


CustomImage CustomImage::cutImageOut(const Rect& box,string inName)
{
	CustomImage boxImage;
	boxImage.setName(inName);

	Mat bfC = m_brightfieldChannel;
	Mat nC = m_nucleusChannel;
	Mat cC = m_cytoplasmChannel;
	Mat ncC = m_yapChannel;
	bfC = bfC(box).clone();
	nC = nC(box).clone();
	cC = cC(box).clone();
	ncC = ncC(box).clone();
	boxImage.setChannel(bfC,channelType::brightfield);
	boxImage.setChannel(nC, channelType::nucleus);
	boxImage.setChannel(cC, channelType::cytoplasm);
	boxImage.setChannel(ncC, channelType::yap);
	return boxImage;
}

// writes into RGB => R=brightfield, G=Cytoplasm, B=Nucleus
Mat CustomImage::createRGBimage()
{
	Mat rgbImage, imgRed, imgGreen, imgBlue;
	vector<Mat> allChannels;
	imgRed = m_brightfieldChannel;
	imgGreen = m_cytoplasmChannel;
	imgBlue = m_nucleusChannel;

	allChannels.push_back(imgBlue);
	allChannels.push_back(imgGreen);
	allChannels.push_back(imgRed);

	merge(allChannels, rgbImage);
	//Mat channels[3] = { imgBlue, imgGreen,  imgRed};
	//merge(channels, 3, rgbImage);
	return rgbImage;
}

