#include "WindowChannels.h"

using namespace std;


WindowChannels::WindowChannels(QWidget* parent) :
    QDialog(parent),
    ui(new Ui::WindowChannels)
{
    ui->setupUi(this);
}

WindowChannels::~WindowChannels()
{
    delete ui;
}

vector<channelType> WindowChannels::set_ChannelTypes()
{
    vector<channelType> channels;
    channels.push_back(channelType(ui->comboBox_1->currentIndex()));
    channels.push_back(channelType(ui->comboBox_2->currentIndex()));
    channels.push_back(channelType(ui->comboBox_3->currentIndex()));
    channels.push_back(channelType(ui->comboBox_4->currentIndex()));
    channels.push_back(channelType(ui->comboBox_5->currentIndex()));

    //removes None type
    channels.erase(std::remove(channels.begin(), channels.end(), channelType::None), channels.end());

    return channels;
}
