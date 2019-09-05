#ifndef DPLWINDOW_H
#define DPLWINDOW_H

#include <QMainWindow>

namespace Ui {
class DPLWindow;
}

class DPLWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit DPLWindow(QWidget *parent = 0);
    ~DPLWindow();

private:
    Ui::DPLWindow *ui;
};

#endif // DPLWINDOW_H
