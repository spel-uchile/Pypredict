#include "dplwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    DPLWindow w;
    w.show();

    return a.exec();
}
