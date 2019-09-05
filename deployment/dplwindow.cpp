#include "dplwindow.h"
#include "ui_dpl_window.h"

DPLWindow::DPLWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::DPLWindow)
{
    ui->setupUi(this);
}

DPLWindow::~DPLWindow()
{
    delete ui;
}
