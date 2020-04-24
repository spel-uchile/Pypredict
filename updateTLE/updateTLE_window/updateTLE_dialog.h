#ifndef UPDATETLE_DIALOG_H
#define UPDATETLE_DIALOG_H

#include <QDialog>

namespace Ui {
class updateTLE_dialog;
}

class updateTLE_dialog : public QDialog
{
    Q_OBJECT

public:
    explicit updateTLE_dialog(QWidget *parent = nullptr);
    ~updateTLE_dialog();

private:
    Ui::updateTLE_dialog *ui;
};

#endif // UPDATETLE_DIALOG_H
