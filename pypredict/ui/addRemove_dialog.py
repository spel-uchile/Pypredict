# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'addRemove_window.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_addRemove(object):
    def setupUi(self, addRemove):
        addRemove.setObjectName("addRemove")
        addRemove.resize(574, 308)
        self.horizontalLayout = QtWidgets.QHBoxLayout(addRemove)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.frame = QtWidgets.QFrame(addRemove)
        self.frame.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.frame)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.frame_4 = QtWidgets.QFrame(self.frame)
        self.frame_4.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_4.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_4.setObjectName("frame_4")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.frame_4)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label = QtWidgets.QLabel(self.frame_4)
        self.label.setObjectName("label")
        self.horizontalLayout_2.addWidget(self.label)
        self.srch_box = QtWidgets.QLineEdit(self.frame_4)
        self.srch_box.setObjectName("srch_box")
        self.horizontalLayout_2.addWidget(self.srch_box)
        self.verticalLayout.addWidget(self.frame_4)
        self.avail_sats_lst = QtWidgets.QListWidget(self.frame)
        self.avail_sats_lst.setObjectName("avail_sats_lst")
        self.verticalLayout.addWidget(self.avail_sats_lst)
        self.horizontalLayout.addWidget(self.frame)
        self.frame_2 = QtWidgets.QFrame(addRemove)
        self.frame_2.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame_2)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem)
        self.add_sat_bt = QtWidgets.QPushButton(self.frame_2)
        self.add_sat_bt.setObjectName("add_sat_bt")
        self.verticalLayout_2.addWidget(self.add_sat_bt)
        self.remove_sat_bt = QtWidgets.QPushButton(self.frame_2)
        self.remove_sat_bt.setObjectName("remove_sat_bt")
        self.verticalLayout_2.addWidget(self.remove_sat_bt)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem1)
        self.horizontalLayout.addWidget(self.frame_2)
        self.frame_3 = QtWidgets.QFrame(addRemove)
        self.frame_3.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.frame_3)
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.curr_sats_tab = QtWidgets.QTableWidget(self.frame_3)
        self.curr_sats_tab.setObjectName("curr_sats_tab")
        self.curr_sats_tab.setColumnCount(1)
        self.curr_sats_tab.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignVCenter)
        self.curr_sats_tab.setHorizontalHeaderItem(0, item)
        self.curr_sats_tab.horizontalHeader().setStretchLastSection(True)
        self.curr_sats_tab.verticalHeader().setStretchLastSection(False)
        self.verticalLayout_3.addWidget(self.curr_sats_tab)
        self.horizontalLayout.addWidget(self.frame_3)

        self.retranslateUi(addRemove)
        QtCore.QMetaObject.connectSlotsByName(addRemove)

    def retranslateUi(self, addRemove):
        _translate = QtCore.QCoreApplication.translate
        addRemove.setWindowTitle(_translate("addRemove", "Dialog"))
        self.label.setText(_translate("addRemove", "Search:"))
        self.add_sat_bt.setText(_translate("addRemove", "→"))
        self.remove_sat_bt.setText(_translate("addRemove", "←"))
        item = self.curr_sats_tab.horizontalHeaderItem(0)
        item.setText(_translate("addRemove", "Current satellites"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    addRemove = QtWidgets.QDialog()
    ui = Ui_addRemove()
    ui.setupUi(addRemove)
    addRemove.show()
    sys.exit(app.exec_())
