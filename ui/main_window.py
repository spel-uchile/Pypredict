# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main_window.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1263, 713)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.tab)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.frame_2 = QtWidgets.QFrame(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_2.sizePolicy().hasHeightForWidth())
        self.frame_2.setSizePolicy(sizePolicy)
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.frame_2)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.dpl_button = QtWidgets.QPushButton(self.frame_2)
        self.dpl_button.setObjectName("dpl_button")
        self.horizontalLayout_3.addWidget(self.dpl_button)
        self.loc_button = QtWidgets.QPushButton(self.frame_2)
        self.loc_button.setObjectName("loc_button")
        self.horizontalLayout_3.addWidget(self.loc_button)
        self.prev_day = QtWidgets.QToolButton(self.frame_2)
        self.prev_day.setObjectName("prev_day")
        self.horizontalLayout_3.addWidget(self.prev_day)
        self.prev_min = QtWidgets.QToolButton(self.frame_2)
        self.prev_min.setObjectName("prev_min")
        self.horizontalLayout_3.addWidget(self.prev_min)
        self.play = QtWidgets.QToolButton(self.frame_2)
        self.play.setObjectName("play")
        self.horizontalLayout_3.addWidget(self.play)
        self.next_min = QtWidgets.QToolButton(self.frame_2)
        self.next_min.setObjectName("next_min")
        self.horizontalLayout_3.addWidget(self.next_min)
        self.next_day = QtWidgets.QToolButton(self.frame_2)
        self.next_day.setObjectName("next_day")
        self.horizontalLayout_3.addWidget(self.next_day)
        self.datetime = QtWidgets.QDateTimeEdit(self.frame_2)
        self.datetime.setObjectName("datetime")
        self.horizontalLayout_3.addWidget(self.datetime)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.verticalLayout_4.addWidget(self.frame_2)
        self.frame_plot = QtWidgets.QFrame(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_plot.sizePolicy().hasHeightForWidth())
        self.frame_plot.setSizePolicy(sizePolicy)
        self.frame_plot.setMaximumSize(QtCore.QSize(1180, 16777215))
        self.frame_plot.setObjectName("frame_plot")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.frame_plot)
        self.verticalLayout.setObjectName("verticalLayout")
        self.verticalLayout_4.addWidget(self.frame_plot)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.tab_2)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.Table = QtWidgets.QTableWidget(self.tab_2)
        self.Table.setColumnCount(13)
        self.Table.setObjectName("Table")
        self.Table.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.Table.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(5, item)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(6, item)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(7, item)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(8, item)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(9, item)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(10, item)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(11, item)
        item = QtWidgets.QTableWidgetItem()
        self.Table.setHorizontalHeaderItem(12, item)
        self.Table.horizontalHeader().setDefaultSectionSize(88)
        self.Table.horizontalHeader().setMinimumSectionSize(20)
        self.Table.horizontalHeader().setSortIndicatorShown(True)
        self.Table.horizontalHeader().setStretchLastSection(True)
        self.verticalLayout_5.addWidget(self.Table)
        self.tabWidget.addTab(self.tab_2, "")
        self.verticalLayout_2.addWidget(self.tabWidget)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1263, 22))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        self.menuEdit = QtWidgets.QMenu(self.menubar)
        self.menuEdit.setObjectName("menuEdit")
        self.menuView = QtWidgets.QMenu(self.menubar)
        self.menuView.setObjectName("menuView")
        self.menuChange_planet = QtWidgets.QMenu(self.menuView)
        self.menuChange_planet.setObjectName("menuChange_planet")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionExit = QtWidgets.QAction(MainWindow)
        self.actionExit.setObjectName("actionExit")
        self.actionAbout = QtWidgets.QAction(MainWindow)
        self.actionAbout.setObjectName("actionAbout")
        self.actionOpen = QtWidgets.QAction(MainWindow)
        self.actionOpen.setObjectName("actionOpen")
        self.actionSave_as = QtWidgets.QAction(MainWindow)
        self.actionSave_as.setObjectName("actionSave_as")
        self.actionUpdate_TLE_from_net = QtWidgets.QAction(MainWindow)
        self.actionUpdate_TLE_from_net.setObjectName("actionUpdate_TLE_from_net")
        self.actionUpdate_TLE_from_file = QtWidgets.QAction(MainWindow)
        self.actionUpdate_TLE_from_file.setObjectName("actionUpdate_TLE_from_file")
        self.actionEnable_database = QtWidgets.QAction(MainWindow)
        self.actionEnable_database.setObjectName("actionEnable_database")
        self.actionAdd_south_atlantic_anomaly = QtWidgets.QAction(MainWindow)
        self.actionAdd_south_atlantic_anomaly.setObjectName("actionAdd_south_atlantic_anomaly")
        self.actionRemove_coverage = QtWidgets.QAction(MainWindow)
        self.actionRemove_coverage.setObjectName("actionRemove_coverage")
        self.actionAdd_remove_satellites = QtWidgets.QAction(MainWindow)
        self.actionAdd_remove_satellites.setObjectName("actionAdd_remove_satellites")
        self.actionFullscreen = QtWidgets.QAction(MainWindow)
        self.actionFullscreen.setObjectName("actionFullscreen")
        self.actionEarth = QtWidgets.QAction(MainWindow)
        self.actionEarth.setObjectName("actionEarth")
        self.actionMars = QtWidgets.QAction(MainWindow)
        self.actionMars.setObjectName("actionMars")
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionSave_as)
        self.menuFile.addAction(self.actionExit)
        self.menuHelp.addAction(self.actionAbout)
        self.menuEdit.addAction(self.actionUpdate_TLE_from_net)
        self.menuEdit.addAction(self.actionUpdate_TLE_from_file)
        self.menuEdit.addAction(self.actionEnable_database)
        self.menuEdit.addSeparator()
        self.menuEdit.addAction(self.actionAdd_south_atlantic_anomaly)
        self.menuEdit.addAction(self.actionRemove_coverage)
        self.menuEdit.addAction(self.actionAdd_remove_satellites)
        self.menuChange_planet.addAction(self.actionEarth)
        self.menuChange_planet.addAction(self.actionMars)
        self.menuView.addAction(self.actionFullscreen)
        self.menuView.addAction(self.menuChange_planet.menuAction())
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())
        self.menubar.addAction(self.menuView.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        self.actionExit.triggered.connect(MainWindow.close)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "PyPredict"))
        self.dpl_button.setText(_translate("MainWindow", "Simulate deployment"))
        self.loc_button.setText(_translate("MainWindow", "Simulate localization"))
        self.prev_day.setText(_translate("MainWindow", "|◀◀"))
        self.prev_min.setText(_translate("MainWindow", "|◀"))
        self.play.setText(_translate("MainWindow", "▶️"))
        self.next_min.setText(_translate("MainWindow", "▶️|"))
        self.next_day.setText(_translate("MainWindow", "▶️▶️|"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Tab 1"))
        self.Table.setSortingEnabled(True)
        item = self.Table.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Satellite"))
        item = self.Table.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "Category"))
        item = self.Table.horizontalHeaderItem(2)
        item.setText(_translate("MainWindow", "Latitude"))
        item = self.Table.horizontalHeaderItem(3)
        item.setText(_translate("MainWindow", "Longitude"))
        item = self.Table.horizontalHeaderItem(4)
        item.setText(_translate("MainWindow", "Alt. [km]"))
        item = self.Table.horizontalHeaderItem(5)
        item.setText(_translate("MainWindow", "Spd. [m/s]"))
        item = self.Table.horizontalHeaderItem(6)
        item.setText(_translate("MainWindow", "a [km]"))
        item = self.Table.horizontalHeaderItem(7)
        item.setText(_translate("MainWindow", "h [km²/s]"))
        item = self.Table.horizontalHeaderItem(8)
        item.setText(_translate("MainWindow", "e"))
        item = self.Table.horizontalHeaderItem(9)
        item.setText(_translate("MainWindow", "Ω"))
        item = self.Table.horizontalHeaderItem(10)
        item.setText(_translate("MainWindow", "i"))
        item = self.Table.horizontalHeaderItem(11)
        item.setText(_translate("MainWindow", "ω"))
        item = self.Table.horizontalHeaderItem(12)
        item.setText(_translate("MainWindow", "T. Anom."))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Tab 2"))
        self.menuFile.setTitle(_translate("MainWindow", "Fi&le"))
        self.menuHelp.setTitle(_translate("MainWindow", "Help"))
        self.menuEdit.setTitle(_translate("MainWindow", "E&dit"))
        self.menuView.setTitle(_translate("MainWindow", "&View"))
        self.menuChange_planet.setTitle(_translate("MainWindow", "&Change planet"))
        self.actionExit.setText(_translate("MainWindow", "&Exit"))
        self.actionAbout.setText(_translate("MainWindow", "&About"))
        self.actionOpen.setText(_translate("MainWindow", "&Open"))
        self.actionSave_as.setText(_translate("MainWindow", "&Save as"))
        self.actionUpdate_TLE_from_net.setText(_translate("MainWindow", "&Update TLE from net"))
        self.actionUpdate_TLE_from_file.setText(_translate("MainWindow", "Update TLE &from file"))
        self.actionEnable_database.setText(_translate("MainWindow", "&Enable database"))
        self.actionAdd_south_atlantic_anomaly.setText(_translate("MainWindow", "&Add south atlantic anomaly"))
        self.actionRemove_coverage.setText(_translate("MainWindow", "&Remove coverage"))
        self.actionAdd_remove_satellites.setText(_translate("MainWindow", "A&dd/remove satellites"))
        self.actionFullscreen.setText(_translate("MainWindow", "&Fullscreen"))
        self.actionEarth.setText(_translate("MainWindow", "&Earth"))
        self.actionMars.setText(_translate("MainWindow", "&Mars"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
