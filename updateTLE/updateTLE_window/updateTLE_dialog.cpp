#include "updateTLE_dialog.h"
#include "ui_updateTLE_dialog.h"

updateTLE_dialog::updateTLE_dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::updateTLE_dialog)
{
    ui->setupUi(this);
}

updateTLE_dialog::~updateTLE_dialog()
{
    delete ui;
}
