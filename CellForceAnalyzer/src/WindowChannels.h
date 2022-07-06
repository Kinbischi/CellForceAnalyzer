#pragma once

#include <QDialog>
#include "ui_WindowChannels.h"

#include "CustomImage.h"

namespace Ui {
    class WindowChannels;
}

class WindowChannels : public QDialog
{
    Q_OBJECT

public:
    explicit WindowChannels(QWidget* parent = nullptr);
    ~WindowChannels();
    std::vector<channelType> set_ChannelTypes();

private slots:

    
private:
    Ui::WindowChannels* ui;
};

