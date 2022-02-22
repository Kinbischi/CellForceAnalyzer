#include "WindowMain.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    WindowMain w;
    w.show();
    return a.exec();
}
