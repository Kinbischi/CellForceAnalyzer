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
    //void on_okButton_clicked();
    

    
private:
    //TODO decide whether pointer or not???
    //TODO: give channels as signal slot connection
    Ui::WindowChannels* ui;
};

