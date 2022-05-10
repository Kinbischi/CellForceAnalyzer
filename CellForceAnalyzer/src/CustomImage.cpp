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
	m_actinChannel = other.m_actinChannel.clone();
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
		name += "_cp";
		break;
	case channelType::yap:
		name += "_yap";
		break;
	case channelType::actin:
		name += "_actin";
		break;
	case channelType::None:
		name += "_rgb";
		break;
	}
	return name;
}

void CustomImage::setChannel(Mat img, channelType type)
{

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
	case channelType::actin:
		m_actinChannel = img;
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
	case channelType::actin:
		img = m_actinChannel;
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
	help::thresh(m_actinChannel);
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
	Mat yapC = m_yapChannel;
	Mat acC = m_actinChannel;

	if (!bfC.empty())
	{
		bfC = bfC(box).clone();
	}
	if (!nC.empty())
	{
		nC = nC(box).clone();
	}
	if (!cC.empty())
	{
		cC = cC(box).clone();
	}
	if (!yapC.empty())
	{
		yapC = yapC(box).clone();
	}
	if (!acC.empty())
	{
		acC = acC(box).clone();
	}

	boxImage.setChannel(bfC,channelType::brightfield);
	boxImage.setChannel(nC, channelType::nucleus);
	boxImage.setChannel(cC, channelType::cytoplasm);
	boxImage.setChannel(yapC, channelType::yap);
	boxImage.setChannel(acC, channelType::actin);
	return boxImage;
}

// writes into RGB => R=brightfield, G=actin, B=nucleus
Mat CustomImage::createRGBimage()
{
	Mat rgbImage, imgRed, imgGreen, imgBlue;
	vector<Mat> allChannels;
	imgRed = m_brightfieldChannel;
	imgGreen = m_actinChannel;
	imgBlue = m_nucleusChannel;

	allChannels.push_back(imgBlue);
	allChannels.push_back(imgGreen);
	allChannels.push_back(imgRed);

	merge(allChannels, rgbImage);
	//Mat channels[3] = { imgBlue, imgGreen,  imgRed};
	//merge(channels, 3, rgbImage);
	return rgbImage;
}


