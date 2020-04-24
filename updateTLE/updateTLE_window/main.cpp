#include "updateTLE_dialog.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    updateTLE_dialog w;
    w.show();

    return a.exec();
}
