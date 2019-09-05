#include "helpwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    HelpWindow w;
    w.show();

    return a.exec();
}
