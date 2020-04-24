# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'updateTLE_dialog.ui'
#
# Created by: PyQt5 UI code generator 5.14.2
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_updateTLE(object):
    def setupUi(self, updateTLE):
        updateTLE.setObjectName("updateTLE")
        updateTLE.resize(534, 367)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(updateTLE)
        self.verticalLayout_2.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_2.setSpacing(6)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.plainTextEdit = QtWidgets.QPlainTextEdit(updateTLE)
        self.plainTextEdit.setPlainText("")
        self.plainTextEdit.setObjectName("plainTextEdit")
        self.verticalLayout_2.addWidget(self.plainTextEdit)
        self.progressBar = QtWidgets.QProgressBar(updateTLE)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.verticalLayout_2.addWidget(self.progressBar)
        self.buttonBox = QtWidgets.QDialogButtonBox(updateTLE)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout_2.addWidget(self.buttonBox)

        self.retranslateUi(updateTLE)
        QtCore.QMetaObject.connectSlotsByName(updateTLE)

    def retranslateUi(self, updateTLE):
        _translate = QtCore.QCoreApplication.translate
        updateTLE.setWindowTitle(_translate("updateTLE", "updateTLE_dialog"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    updateTLE = QtWidgets.QDialog()
    ui = Ui_updateTLE()
    ui.setupUi(updateTLE)
    updateTLE.show()
    sys.exit(app.exec_())
