#include "CustomImage.h"

#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>


using namespace cv;
using namespace std;


CustomImage::CustomImage(const CustomImage& other)
{
	m_name = other.m_name;
	m_brightfieldChannel = other.m_brightfieldChannel.clone();
	m_nucleusChannel = other.m_nucleusChannel.clone();
	m_vinculinChannel = other.m_vinculinChannel.clone();
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
	case channelType::vinculin:
		name += "_vin";
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
	case channelType::vinculin:
		m_vinculinChannel = img;
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
	case channelType::vinculin:
		img = m_vinculinChannel;
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

CustomImage CustomImage::cutImageOut(const Rect& box,string inName)
{
	CustomImage boxImage;
	boxImage.setName(inName);

	Mat bfC = m_brightfieldChannel;
	Mat nC = m_nucleusChannel;
	Mat vinC = m_vinculinChannel;
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
	if (!vinC.empty())
	{
		vinC = vinC(box).clone();
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
	boxImage.setChannel(vinC, channelType::vinculin);
	boxImage.setChannel(yapC, channelType::yap);
	boxImage.setChannel(acC, channelType::actin);
	return boxImage;
}

// writes into RGB => R=brightfield, G=actin, B=nucleus
Mat CustomImage::createRGBimage()
{
	Mat rgbImage, imgRed, imgGreen, imgBlue;
	vector<Mat> allChannels;
	if (!m_brightfieldChannel.empty())
	{
		imgRed = m_brightfieldChannel;
	}
	else if(!m_vinculinChannel.empty())
	{
		imgRed = m_vinculinChannel;
	}
	
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


