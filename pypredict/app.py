"""
                                Pypredict
    Orbit prediction software. Displays the satellites' position and
    orbital parameters in real time. Simulates satellite localization
    and deployment.
    
    Copyright (C) 2018-2020, Matías Vidal Valladares, matvidal.
    Authors: Matías Vidal Valladares <matias.vidal.v@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
__version__ = "3.3.1"

from cartopy.crs import Geodetic, PlateCarree
from cartopy.geodesic import Geodesic
from datetime import datetime, timedelta
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.patches import Polygon
from matplotlib.pyplot import imread, figure
from numpy import abs, asarray, cos, ndarray, pi
from pkg_resources import resource_filename
from PyQt5 import QtWidgets, QtGui, QtCore
from pymongo import MongoClient
from pypredict.dayNightMap import Map
from pypredict.dpl import Dpl
from pypredict.navigation import Toolbar
from pypredict.sat import Sat
from pypredict.SAA import SAA
from pypredict.ui.main_window import Ui_MainWindow
from pypredict.ui.about_dialog import Ui_About
from pypredict.ui.addRemove_dialog import Ui_addRemove
from pypredict.ui.dpl_dialog import Ui_DPL
from pypredict.ui.updateTLE_dialog import Ui_updateTLE
from urllib.request import urlopen

class ApplicationWindow(QtWidgets.QMainWindow):

    __slots__ = ["Sats", "img", "mainSat", "mainSat_lats",
                 "mainSat_lngs", "ax_saa", "fig", "ax", "ax_tray",
                 "ax_sat", "ax_cov", "sat_txt", "saa_alpha",
                 "cov_alpha", "popup", "match", "avail_sats", "argos",
                 "beidou", "cubesat", "dmc", "education", "geodetic",
                 "goes", "intelsat", "iridium", "iridium_next",
                 "military", "molniya", "noaa", "oneweb", "planet",
                 "radar", "resource", "sarsat", "spire", "tdrss",
                 "tle_new", "weather", "x_comm", "active", "tle_files",
                 "map", "dpl_img", "tdoa_img", "world_map", "dpl",
                 "dmin", "canvas", "saa", "date", "db", "en_db",
                 "time_timer", "sats_timer", "canvas_timer",
                 "bg_timer", "Dialog", "table_timer", "sats_lngs",
                 "sats_lats", "usr_sats", "usr_tle_file", "toolbar",
                 "img_path", "tle_path", "forward", "backward", "pause"]

    def __init__(self, Sats):
        self.Sats = Sats
        self.sortSats()
        super(ApplicationWindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.img_path = resource_filename("pypredict","img/")
        self.setWindowIcon(QtGui.QIcon("{}favicon.png".format(self.img_path)))
        self.setWindowTitle('Pypredict')
        self.dmin = 0
        self.forward = False
        self.backward = False
        self.pause = False
        self.updateTime()
        self.setButtons()
        self.world_map = Map("{}earth_nasa_day.png".format(self.img_path),
                             "{}earth_nasa_night.png".format(self.img_path))
        self.plotData()
        self.setCanvas()
        self.setCustomStatusBar()
        self.setMenu()
        client = MongoClient("localhost", 27017)
        self.db = client["SatConstellation"]
        self.en_db = False
        self.data_gen()
        self.updateTableContent()
        self.setTableConnections()
        self.showMaximized()
        self.fig.canvas.mpl_connect('scroll_event',self.zoom)
        self.run()

    def __call__(self):
        return self

    def keyPressEvent(self, e):
        """
        Enables the fullscreen mode by pressing the F11 key.
        It exits the fullscreen mode by pressing Esc.

        Parameters
        ----------
        e : event
            Key pressed event.
        """
        if e.key() == QtCore.Qt.Key_Escape:
            self.exitFullscreen()
        if e.key() == QtCore.Qt.Key_F11:
            if self.isFullScreen():
                self.showMaximized()
            else:
                self.showFullScreen()

    def zoom(self, event):
        """
        This method handles the zoom in and zoom out of the world map
        figure. It first centers the view in the middle of the map, and
        then scales it and moves it to the mouse position.

        Parameters
        ----------
        event : event
                Mouse wheel event.
        """
        lng = event.xdata
        lat = event.ydata
        if (event.button == "up"):
            scale = 0.8
        else:
            scale = 1.6
        lng_min2 = -(abs(self.ax.get_xlim()[1]) + abs(self.ax.get_xlim()[0]))/2
        lng_max2 = (abs(self.ax.get_xlim()[1]) + abs(self.ax.get_xlim()[0]))/2
        lat_min2 = -(abs(self.ax.get_ylim()[1]) + abs(self.ax.get_ylim()[0]))/2
        lat_max2 = (abs(self.ax.get_ylim()[1]) + abs(self.ax.get_ylim()[0]))/2
        lng_min = lng_min2*scale + lng
        lng_max = lng_max2*scale + lng
        lat_min = lat_min2*scale + lat
        lat_max = lat_max2*scale + lat
        if (lng_min < -180 and lng_max > 180 or lat_min < -90 and lat_max > 90):
            self.ax.set_extent([-180, 180, -90, 90])
        else:
            if (lng_min < -180):
                lng_max += -(180 + lng_min)
                lng_min = -180
            elif (lng_max > 180):
                lng_min += -(lng_max - 180)
                lng_max = 90
            if (lat_min < -90):
                lat_max += -(90 + lat_min)
                lat_min = -90
            elif (lat_max > 90):
                lat_min += -(lat_max - 90)
                lat_max = 90
            self.ax.set_extent([lng_min, lng_max, lat_min, lat_max])

    def fillSAA(self):
        """
        Plot a fill on the south atlantic anomaly (3 sigma).
        """
        self.saa = SAA()
        self.saa_alpha = 0
        self.ax_saa, = self.ax.fill(self.saa.lng, self.saa.lat,
                                    color='g', alpha=self.saa_alpha)

    def plotData(self):
        """
        Plots the figure used by this program to display the world map.
        It also creates the axes used to display the satellites. The
        world map day and night is plotted according to the program's
        date.
        """
        self.fig = figure(figsize=(16, 8))
        self.ax = self.fig.add_axes([0, 0, 1, 1], projection=PlateCarree(),
                                    frameon=False)
        img_extent = (-180, 180, -90, 90)
        self.map = self.ax.imshow(self.world_map.fillDarkSideFromPicture(self.date),
                                  origin='upper', extent=img_extent)
        self.gridAndFormat("gray", 0.5, "white", 9)
        self.ax.spines['left'].set_visible(True)
        self.ax.spines['bottom'].set_visible(True)
        self.ax.spines['right'].set_visible(True)
        self.ax.spines['top'].set_visible(True)

    def gridAndFormat(self, gcolor, galpha, tcolor, tsize):
        """
        This method creates the grid and the labels of the angles. It
        gives has the option to change some properties of the grid and
        labels.

        Parameters
        ----------
        gcolor : str
                 Color of the grid lines.
        galpha : float
                 Alpha value of the grid lines.
        tcolor : str
                 Color of the text labels.
        tsize  : int
                 Size of the text labels.
        """
        self.ax.plot([-150, -150], [-90, 90], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([-120, -120], [-90, 90], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([-90, -90], [-90, 90], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([-60, -60], [-90, 90], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([-30, -30], [-90, 90], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([0, 0], [-90, 90], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([30, 30], [-90, 90], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([60, 60], [-90, 90], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([90, 90], [-90, 90], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([120, 120], [-90, 90], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([150, 150], [-90, 90], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([-180, 180], [60, 60], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([-180, 180], [30, 30], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([-180, 180], [0, 0], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([-180, 180], [-30, -30], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.plot([-180, 180], [-60, -60], color=gcolor,
                     alpha=galpha, linewidth=1, linestyle='-')
        self.ax.text(-150, 86.5, "150°W", color=tcolor, size=tsize, ha="center")
        self.ax.text(-120, 86.5, "120°W", color=tcolor, size=tsize, ha="center")
        self.ax.text(-90, 86.5, "90°W", color=tcolor, size=tsize, ha="center")
        self.ax.text(-60, 86.5, "60°W", color=tcolor, size=tsize, ha="center")
        self.ax.text(-30, 86.5, "30°W", color=tcolor, size=tsize, ha="center")
        self.ax.text(0, 86.5, "0°", color=tcolor, size=tsize, ha="center")
        self.ax.text(30, 86.5, "30°E", color=tcolor, size=tsize, ha="center")
        self.ax.text(60, 86.5, "60°E", color=tcolor, size=tsize, ha="center")
        self.ax.text(90, 86.5, "90°E", color=tcolor, size=tsize, ha="center")
        self.ax.text(120, 86.5, "120°E", color=tcolor, size=tsize, ha="center")
        self.ax.text(150, 86.5, "150°E", color=tcolor, size=tsize, ha="center")
        self.ax.text(-179.5, 60, "60°N", color=tcolor, size=tsize, va="center")
        self.ax.text(-179.5, 30, "30°N", color=tcolor, size=tsize, va="center")
        self.ax.text(-178, 0, "0°", color=tcolor, size=tsize, va="center")
        self.ax.text(-179.5, -30, "30°S", color=tcolor, size=tsize, va="center")
        self.ax.text(-179.5, -60, "60°S", color=tcolor, size=tsize, va="center")

    def data_gen(self):
        """
        Generates the objects that will handle the plots of the
        satellites, their coverage, their names and their trajectories.
        It also sets the south atlantic anomaly (SAA) data.
        """
        self.ax_tray, = self.ax.plot([], [], color='green',
                                     linewidth=1.4, linestyle='-',
                                     transform=Geodetic())
        self.ax_sat, = self.ax.plot([], [], 'yo', ms=6)
        self.sats_lngs = []
        self.sats_lats = []
        self.readAllSats()
        self.updateSatellites()
        self.mainSat = self.Sats[0]
        tf = int(self.mainSat.getPeriod()*3)
        lats, lngs = self.mainSat.getTrayectory(tf, 50, self.date)
        self.mainSat_lats, self.mainSat_lngs = lats, lngs
        self.ax_tray.set_data(self.mainSat_lngs, self.mainSat_lats)
        self.canvas.draw_idle()
        self.resetCov()
        self.fillSAA()

    def resetCov(self):
        """
        Resets the axis data of the satellites names and coverage.
        """
        self.ax_cov = []
        self.sat_txt = []
        self.cov_alpha = 0.2
        xy = ndarray(shape=(200,2), dtype=float)
        for i, Sat in enumerate(self.Sats):
            self.ax_cov.append(self.ax.add_patch(Polygon(xy, color='w',
                                                         transform=Geodetic(),
                                                         alpha=self.cov_alpha)))
            self.sat_txt.append(self.ax.text(Sat.getLng(date=self.date),
                                             Sat.getLat(), Sat.name, c="yellow",
                                             size=7, ha="center", va="center"))
            self.sats_lngs.append(Sat.getLng(date=self.date))
            self.sats_lats.append(Sat.getLat())

    def updateCanvas(self):
        """
        Updates the plots displayed in the canvas, such as the
        satellites names, position and coverage. This method
        is executed only if the World Map tab is selected.
        """
        tab_index = self.ui.tabWidget.currentIndex()
        worldTab = self.ui.tabWidget.tabText(tab_index) == "World Map"
        if (worldTab and self.pause is False):
            for i, Sat in  enumerate(self.Sats):
                self.sats_lngs[i] = Sat.getLng(date=self.date)
                lng = self.sats_lngs[i] - 6*(self.sats_lngs[i] > 173)
                lng += 6*(self.sats_lngs[i] < -173)
                self.sats_lats[i] = Sat.getLat()
                lat = self.sats_lats[i] - 3.5 + 6*(self.sats_lats[i] < -85)
                self.sat_txt[i].set_position((lng, lat))
                if (self.cov_alpha > 0):
                    coverage = Sat.getCoverage()*0.017453292519943295 # to rads
                    radius = coverage*Sat.getPlanetRadius()
                    xy_data = asarray(Geodesic().circle(lon=self.sats_lngs[i],
                                                        lat=self.sats_lats[i],
                                                        radius=radius,
                                                        n_samples=200,
                                                        endpoint=True))
                    self.ax_cov[i].set_xy(xy_data)
            self.ax_sat.set_data(self.sats_lngs, self.sats_lats)
            self.canvas.draw_idle()

    def changeMainSat(self, row):
        """
        It changes the program's main satellite. This satellite is the
        one which trayectory is being displayed.

        Parameters
        ----------
        row : int
              The number of the table's row.
        """
        selectedSat = self.ui.Table.item(row, 0).text()
        for sat in self.Sats:
            if (sat.name == selectedSat):
                self.mainSat = sat
        tf = int(self.mainSat.getPeriod()*3)
        lats, lngs = self.mainSat.getTrayectory(tf, 50, self.date)
        self.mainSat_lats, self.mainSat_lngs = lats, lngs
        self.ax_tray.set_data(self.mainSat_lngs, self.mainSat_lats)
        self.canvas.draw_idle()

    def setButtons(self):
        """
        Conects the program's buttons to the associated methods. Adds
        tool tips and icons.
        """
        self.ui.download_button.clicked.connect(self.updateTLEDialog)
        self.ui.cov_button.clicked.connect(self.removeCoverage)
        self.ui.sat_button.clicked.connect(self.addRemoveSat)
        self.ui.dpl_button.clicked.connect(self.deployPopup)
        self.ui.loc_button.clicked.connect(self.notAvailable)
        self.ui.backward.clicked.connect(self.speed_up_backward)
        self.ui.play.clicked.connect(self.pause_time)
        self.ui.stop.clicked.connect(self.stop_time)
        self.ui.forward.clicked.connect(self.speed_up_forward)
        self.ui.datetime.dateTimeChanged.connect(self.newDate)
        self.ui.download_button.setToolTip("Update TLE from network")
        self.ui.cov_button.setToolTip("Add/remove sat's coverage")
        self.ui.sat_button.setToolTip("Add/remove satellites")
        self.ui.dpl_button.setToolTip("Simulates the deployment\nof a satellite from another")
        self.ui.loc_button.setToolTip("Simulate localization (not implemented yet)")
        self.ui.backward.setToolTip("Fast backward")
        self.ui.play.setToolTip("Resume or pause the time")
        self.ui.stop.setToolTip("Back to current time")
        self.ui.forward.setToolTip("Fast forward")
        self.ui.loc_button.setIcon(QtGui.QIcon("{}locicon.png".format(self.img_path)))

    def deployPopup(self):
        """
        Instantiates the deployment dialog. Connects the dialog's
        buttons to their respective methods.
        """
        self.Dialog = QtWidgets.QDialog()
        self.popup = Ui_DPL()
        self.popup.setupUi(self.Dialog)
        self.Dialog.setWindowTitle("Deployment settings")
        self.showCurrentSats()
        self.popup.curr_sats_tab.itemClicked.connect(self.selectDeployer)
        self.popup.deploy_now_bt.clicked.connect(self.deploySat)
        self.dpl = Dpl()
        self.Dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.Dialog.exec_()

    def selectDeployer(self, item):
        """
        Obtains the deployer's name that the user wrote in the table.
        Sets this name in the deployer's label.

        Parameters
        ----------
        item : QTableWidgetItem
               The table's item selected by the user. In
               this case, the deployer's name.
        """
        deployer_name = item.text()
        self.popup.deployer_name_lbl.setText(deployer_name)

    def deploySat(self):
        """
        Gets all the data the user input to simulate the deployment of
        a satellite from another one. The new satellite is created and
        added to Sats list. A new axes fill and text is created for the
        new satellite's coverage and name.
        """
        deployer_name = self.popup.deployer_name_lbl.text()
        dplyr_mass = float(self.popup.mass_box1.text())
        dplyd_mass = float(self.popup.mass_box2.text())
        spdx = int(self.popup.spdx_box.text())
        spdy = int(self.popup.spdy_box.text())
        spdz = int(self.popup.spdz_box.text())
        name = self.popup.name_box.text()
        cat = self.popup.cat_box.text()
        for sat in self.Sats:
            if (sat.name == deployer_name):
                deployer = sat
                break
        newSat = self.dpl.deploy(cat, deployer, dplyr_mass,
                dplyd_mass, name, [spdx, spdy, spdz], date=self.date)
        xy = ndarray(shape=(200,2), dtype=float)
        self.ax_cov.append(self.ax.add_patch(Polygon(xy, color='w',
                                                     transform=Geodetic(),
                                                     alpha=self.cov_alpha)))
        self.sat_txt.append(self.ax.text([], [], "", color='yellow', size=7,
                            ha="center", va="center"))
        self.Sats.append(newSat)
        self.sortSats()
        self.sats_lngs.append(0.0)
        self.sats_lats.append(0.0)
        for i, Sat in enumerate(self.Sats):
            self.sat_txt[i].set_text(Sat.name)
            self.sats_lngs[i] = Sat.getLng(date=self.date)
            self.sats_lats[i] = Sat.getLat()
        self.Dialog.close()

    def speed_up_backward(self):
        """
        Enables the program to speed up the time into the past.
        """
        self.forward = False
        self.backward = True
        self.pause = False
        self.format_dt()

    def resume_time(self):
        """
        Sets the delta minutes to zero, so that the program's date is
        the current date.
        """
        self.forward = False
        self.backward = False
        self.pause = False
        self.ui.play.setText("| |")
        self.ui.play.clicked.connect(self.pause_time)
        self.format_dt()

    def pause_time(self):
        "It pauses the program's time."
        self.forward = False
        self.backward = False
        self.pause = True
        self.ui.play.setText("▶️")
        self.ui.play.clicked.connect(self.resume_time)
        self.format_dt()

    def stop_time(self):
        "Returns the time to the present."
        self.forward = False
        self.backward = False
        self.pause = False
        self.dmin = 0
        self.ui.play.setText("| |")
        self.ui.play.clicked.connect(self.pause_time)
        self.format_dt()

    def speed_up_forward(self):
        """
        Enables the program to speed up the time into the future.
        """
        self.forward = True
        self.backward = False
        self.pause = False
        self.format_dt()

    def format_dt(self):
        """
        Updates the program's date and sets the new date value into the
        QDateTimeEdit widget. It also refreshes the background image of
        the map.
        """
        self.updateTime()
        oldState = self.ui.datetime.blockSignals(True)
        self.ui.datetime.setDateTime(self.date)
        self.ui.datetime.blockSignals(oldState)
        self.refreshBackgroundImg()

    def updateTime(self, date=None):
        """
        Updates the program's date by adding the delta minutes, dmin,
        that may be changed by the user. The new date is displayed in
        the QDateTimeEdit widget.

        Parameters
        ----------
        date = datetime
               Datetime object to set the new date.
        """
        if (self.pause is False):
            if (self.forward is True):
                self.dmin += 2
                self.refreshBackgroundImg()
            if (self.backward is True):
                self.dmin -= 2
                self.refreshBackgroundImg()
            if (date is None):
                self.date = datetime.utcnow() + timedelta(minutes=self.dmin)
            else:
                self.date = date
            seconds = self.date.second# + self.date.microsecond*0.000001
            qdt = QtCore.QDateTime(self.date.year, self.date.month,
                                   self.date.day,  self.date.hour,
                                   self.date.minute, seconds)
            oldState = self.ui.datetime.blockSignals(True)
            self.ui.datetime.setDateTime(qdt)
            self.ui.datetime.blockSignals(oldState)

    def newDate(self):
        """
        This method is called when the user changes the value of the
        QDateTimeEdit widget. When this happens, the new date is
        obtained to update the program's date. If the time difference is
        greater than one minute, the background image is refreshed.
        """
        date = self.ui.datetime.dateTime().toPyDateTime()
        time_diff = (date - self.date).total_seconds()/60.0
        self.dmin += time_diff
        self.date = date
        if (abs(time_diff) >= 1):
            self.refreshBackgroundImg()

    def updateSatellites(self):
        """
        Updates the orbital parameters of all the satellites that are on
        the list Sats. If the user enabled the database, part of the
        data is dumped into the database.
        """
        if (self.pause is False):
            for Sat in self.Sats:
                Sat.updateOrbitalParameters(self.date)
                if (self.en_db):
                    col = self.db[Sat.name]
                    col.insert_one(self.formatDump(Sat))
 
    def setCanvas(self):
        """
        Creates a canvas object to display the plot of the background
        image and all the axis artists. This canvas is added as a widget
        into the QT user interface.
        """
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setParent(self)
        self.canvas.draw_idle()
        self.ui.frame_plot.layout().addWidget(self.canvas)

    def setCustomStatusBar(self):
        """
        Moves some buttons into the status bar. Adds a navigation
        toolbar into the status bar. Formats the displayed coordinates
        of the cursor.
        """
        self.toolbar = Toolbar(self.canvas, self)
        self.ui.statusbar.layout().addWidget(self.ui.download_button)
        self.ui.statusbar.layout().addWidget(self.ui.cov_button)
        self.ui.statusbar.layout().addWidget(self.ui.sat_button)
        self.ui.statusbar.layout().addWidget(self.ui.dpl_button)
        self.ui.statusbar.layout().addWidget(self.ui.loc_button)
        self.ui.statusbar.layout().addWidget(self.ui.backward)
        self.ui.statusbar.layout().addWidget(self.ui.play)
        self.ui.statusbar.layout().addWidget(self.ui.stop)
        self.ui.statusbar.layout().addWidget(self.ui.forward)
        self.ui.statusbar.layout().addWidget(self.ui.datetime)
        self.ui.statusbar.layout().addWidget(self.toolbar)
        self.ax.format_coord = self.formatCoordinates
        self.map.format_cursor_data = self.formatCursorData

    def formatCoordinates(self, lng, lat):
        """
        Overrrides the default format to show only latitude and
        longitude.
        """
        return "Lon: {:6.1f}, Lat: {:5.1f}".format(lng, lat)

    def formatCursorData(self, data):
        """
        Removes the RGBA data.
        """
        return ""

    def setTableConnections(self):
        """
        Connects the cell clicks of the user with the method that
        changes the main satellites and displays its future trayectory.
        """
        self.ui.Table.cellClicked.connect(self.changeMainSat)

    def updateTableContent(self):
        """
        If the user selects the Table tab, then all the current
        satellites are displayed in a table with their respective
        orbital parameters, names and categories.
        """
        tab_index = self.ui.tabWidget.currentIndex()
        if (self.ui.tabWidget.tabText(tab_index) == "Table"):
            rad2deg = 180/pi
            self.ui.Table.setRowCount(len(self.Sats))
            for i, Sat in enumerate(self.Sats):
                self.ui.Table.setItem(i, 0, QtWidgets.QTableWidgetItem(Sat.name))
                self.ui.Table.setItem(i, 1, QtWidgets.QTableWidgetItem(Sat.getCategory()))
                Lat = QtWidgets.QTableWidgetItem("{:0.4f}°".format(Sat.getLat()))
                Lat.setTextAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
                self.ui.Table.setItem(i, 2, Lat)
                Lng = QtWidgets.QTableWidgetItem("{:0.4f}°".format(Sat.getLng(date=self.date)))
                Lng.setTextAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
                self.ui.Table.setItem(i, 3, Lng)
                Alt = QtWidgets.QTableWidgetItem("{:0.1f}".format((Sat.getAlt()*0.001)))
                Alt.setTextAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
                self.ui.Table.setItem(i, 4, Alt)
                Spd = QtWidgets.QTableWidgetItem("{}".format(int(Sat.getSpeed())))
                Spd.setTextAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
                self.ui.Table.setItem(i, 5, Spd)
                a = QtWidgets.QTableWidgetItem("{:0.1f}".format((Sat.getSemiMajorAxis()*0.001)))
                a.setTextAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
                self.ui.Table.setItem(i, 6, a)
                H = QtWidgets.QTableWidgetItem("{:0.1f}".format(Sat.getSpecAngMomentum()*0.000001))
                H.setTextAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
                self.ui.Table.setItem(i, 7, H)
                e = QtWidgets.QTableWidgetItem("{:0.4f}".format(Sat.getEccentricity()))
                e.setTextAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
                self.ui.Table.setItem(i, 8, e)
                RAAN = QtWidgets.QTableWidgetItem("{:0.2f}°".format((Sat.getRAAN()*rad2deg)))
                RAAN.setTextAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
                self.ui.Table.setItem(i, 9, RAAN)
                incl = QtWidgets.QTableWidgetItem("{:0.2f}°".format((Sat.getInclination()*rad2deg)))
                incl.setTextAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
                self.ui.Table.setItem(i, 10, incl)
                w = QtWidgets.QTableWidgetItem("{:0.2f}°".format((Sat.getArgPerigee()*rad2deg)))
                w.setTextAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
                self.ui.Table.setItem(i, 11, w)
                theta = QtWidgets.QTableWidgetItem("{:0.2f}°".format((Sat.getAnomaly()*rad2deg)))
                theta.setTextAlignment(QtCore.Qt.AlignVCenter | QtCore.Qt.AlignRight)
                self.ui.Table.setItem(i, 12, theta)

    def formatDump(self, Sat):
        """
        Formats some basic data of a satellite for a database.

        Parameters
        ----------
        Sat : pypredict.sat.Sat object
              Object of the class Sat.
        """
        dump = {
                "id": Sat.satnumber,
                "name": Sat.name,
                "obc": {
                    "datetime": self.date.strftime("%F %H:%M:%S.%f")
                        },
                "orbit": {
                    "lat": Sat.getLat(),
                    "lon": Sat.getLng(date=self.date),
                    "alt": Sat.getAlt()
                        }
                }
        return dump

    def savePicture(self):
        """
        Creates a dialog for the user to save a picture of the world map
        with all the current satellites.
        """
        file_name = QtWidgets.QFileDialog.getSaveFileName(self, "Save picture", "",
                                                          "Images (*.png *.xpm *.jpg)")
        if file_name[0] is None:
            return
        self.fig.savefig(file_name[0], bbox_inches="tight", pad_inches=0)

    def saveTLEs(self):
        """
        Creates a dialog for the user to save all the TLE data of the
        satellites that are being displayed. This data is saved as a txt
        file in the folder that the user decides.
        """
        dlog_name = "Save TLE data into file"
        file_name = QtWidgets.QFileDialog.getSaveFileName(self, dlog_name, "",
                                                          "Text files (*.txt)")
        if file_name[0] == '':
            return
        with open(file_name[0], "w") as f:
            for Sat in self.Sats:
                f.write(Sat.getTLE())
                f.write("\n")

    def notAvailable(self):
        print("This command is not available yet")

    def enableDB(self):
        """
        Enables some basic data of the satellites to be saved in a
        database.
        """
        self.en_db = True
        self.ui.actionEnable_database.setText("Disable database")
        self.ui.actionEnable_database.triggered.connect(self.disableDB)

    def disableDB(self):
        """
        Disables the option to save some data of the satellites in a
        database.
        """
        self.en_db = False
        self.ui.actionEnable_database.setText("Enable database")
        self.ui.actionEnable_database.triggered.connect(self.enableDB)

    def removeSAA(self):
        """
        Removes the South Atlantic Anomaly (SAA) from the map by setting
        its alpha value to zero.
        """
        self.saa_alpha = 0
        self.ax_saa.set_alpha(self.saa_alpha)
        self.canvas.draw_idle()
        self.ui.actionAdd_south_atlantic_anomaly.setText("Add south atlantic anomaly")
        self.ui.actionAdd_south_atlantic_anomaly.triggered.connect(self.addSAA)

    def addSAA(self):
        """
        Displays the South Atlantic Anomaly (SAA) on the map by setting
        its alpha value to a number more than zero.
        """
        self.saa_alpha = 0.2
        self.ax_saa.set_alpha(self.saa_alpha)
        self.canvas.draw_idle()
        self.ui.actionAdd_south_atlantic_anomaly.setText("Remove south atlantic anomaly")
        self.ui.actionAdd_south_atlantic_anomaly.triggered.connect(self.removeSAA)

    def removeCoverage(self):
        """
        Removes the satellites coverage from the plot by setting their
        alpha value to zero.
        """
        self.cov_alpha = 0
        for i, Sat in enumerate(self.Sats):
            self.ax_cov[i].set_alpha(self.cov_alpha)
        self.canvas.draw_idle()
        self.ui.actionRemove_coverage.setText("Add coverage")
        self.ui.actionRemove_coverage.triggered.connect(self.addCoverage)
        self.ui.cov_button.clicked.connect(self.addCoverage)

    def addCoverage(self):
        """
        Displays the satellites coverage on the plot by setting their
        alpha value to a number greater than zero.
        """
        self.cov_alpha = 0.2
        for i, Sat in enumerate(self.Sats):
            self.ax_cov[i].set_alpha(self.cov_alpha)
        self.canvas.draw_idle()
        self.ui.actionRemove_coverage.setText("Remove coverage")
        self.ui.actionRemove_coverage.triggered.connect(self.removeCoverage)
        self.ui.cov_button.clicked.connect(self.removeCoverage)

    def searchSat(self, text):
        """
        Uses the text from the search box to search for a satellite's
        name that matches that text. All the matches are displayed in a
        list.

        Parameters
        ----------
        text : str
               String input of the user used to
               search for a satellite's name.
        """
        srch = text.upper()
        self.match = [s for s in self.avail_sats if srch in s]
        self.match.sort(reverse=False)
        self.popup.avail_sats_lst.clear()
        for sat in self.match:
            self.popup.avail_sats_lst.addItem(sat)

    def readSatsFromFile(self, file, lst):
        """
        Gets all the satellites' names from a file and adds them in a
        list.

        Parameters
        ----------
        file : str
               Path to the file that contains the
               TLE data.
        lst  : list
               List used to append the satellite's
               names from the file.
        """
        with open(file, 'r') as f:
            for count, line in enumerate(f):
                if (count % 3 == 0):
                    lst.append(line.strip())

    def readAllSats(self):
        """
        Gets all the satellites' names from some files that are on the
        TLE folder. Creates lists for the satellites of every file and
        also a list of the satellites of all the files.
        """
        self.active = []
        self.argos = []
        self.beidou = []
        self.cubesat = []
        self.dmc = []
        self.education = []
        self.geodetic = []
        self.goes = []
        self.intelsat = []
        self.iridium = []
        self.iridium_next = []
        self.military = []
        self.molniya = []
        self.noaa = []
        self.oneweb = []
        self.planet = []
        self.radar = []
        self.resource = []
        self.sarsat = []
        self.spire = []
        self.starlink = []
        self.tdrss = []
        self.tle_new = []
        self.visual = []
        self.weather = []
        self.x_comm = []
        self.tle_path = resource_filename("pypredict","data/")
        self.readSatsFromFile("{}active.txt".format(self.tle_path), self.active)
        self.readSatsFromFile("{}argos.txt".format(self.tle_path), self.argos)
        self.readSatsFromFile("{}beidou.txt".format(self.tle_path), self.beidou)
        self.readSatsFromFile("{}cubesat.txt".format(self.tle_path), self.cubesat)
        self.readSatsFromFile("{}dmc.txt".format(self.tle_path), self.dmc)
        self.readSatsFromFile("{}education.txt".format(self.tle_path), self.education)
        self.readSatsFromFile("{}geodetic.txt".format(self.tle_path), self.geodetic)
        self.readSatsFromFile("{}goes.txt".format(self.tle_path), self.goes)
        self.readSatsFromFile("{}intelsat.txt".format(self.tle_path), self.intelsat)
        self.readSatsFromFile("{}iridium.txt".format(self.tle_path), self.iridium)
        self.readSatsFromFile("{}iridium-NEXT.txt".format(self.tle_path), self.iridium_next)
        self.readSatsFromFile("{}military.txt".format(self.tle_path), self.military)
        self.readSatsFromFile("{}molniya.txt".format(self.tle_path), self.molniya)
        self.readSatsFromFile("{}noaa.txt".format(self.tle_path), self.noaa)
        self.readSatsFromFile("{}oneweb.txt".format(self.tle_path), self.oneweb)
        self.readSatsFromFile("{}planet.txt".format(self.tle_path), self.planet)
        self.readSatsFromFile("{}radar.txt".format(self.tle_path), self.radar)
        self.readSatsFromFile("{}resource.txt".format(self.tle_path), self.resource)
        self.readSatsFromFile("{}sarsat.txt".format(self.tle_path), self.sarsat)
        self.readSatsFromFile("{}spire.txt".format(self.tle_path), self.spire)
        self.readSatsFromFile("{}starlink.txt".format(self.tle_path), self.starlink)
        self.readSatsFromFile("{}tdrss.txt".format(self.tle_path), self.tdrss)
        self.readSatsFromFile("{}tle-new.txt".format(self.tle_path), self.tle_new)
        self.readSatsFromFile("{}visual.txt".format(self.tle_path), self.visual)
        self.readSatsFromFile("{}weather.txt".format(self.tle_path), self.weather)
        self.readSatsFromFile("{}x-comm.txt".format(self.tle_path), self.x_comm)
        self.avail_sats = self.active + self.argos + self.beidou + self.cubesat
        self.avail_sats += self.dmc + self.education + self.geodetic + self.goes
        self.avail_sats += self.intelsat + self.iridium + self.iridium_next
        self.avail_sats += self.military + self.molniya + self.noaa
        self.avail_sats += self.oneweb + self.planet + self.radar + self.resource
        self.avail_sats += self.sarsat + self.spire + self.starlink + self.tdrss
        self.avail_sats += self.tle_new + self.visual + self.weather + self.x_comm
        self.avail_sats = [*set(self.avail_sats)]
        self.avail_sats.sort()

    def showAvailSats(self):
        """
        Displays all the available satellites on a list in a dialog.
        """
        for sat in self.avail_sats:
            self.popup.avail_sats_lst.addItem(sat)

    def showCurrentSats(self):
        """
        Displays all the current satellites in a list in a dialog.
        """
        self.popup.curr_sats_tab.setRowCount(len(self.Sats))
        for i, sat in enumerate(self.Sats):
            sat_name = QtWidgets.QTableWidgetItem(sat.name)
            self.popup.curr_sats_tab.setItem(i, 0, sat_name)

    def addRemoveButtons(self):
        """
        Connects the buttons of the add remove satellites dialog with
        the associated methods.
        """
        self.popup.add_sat_bt.clicked.connect(self.addSat)
        self.popup.remove_sat_bt.clicked.connect(self.removeSat)

    def sortSats(self):
        """
        Sorts by name all the satellites that are currently displayed.
        """
        self.Sats.sort(key = lambda s: s.name)

    def addSat(self):
        """
        Adds the selected satellite into the list of the current
        satellites. First the TLE is found by name and category, second
        the satellite is generated, third the satellite is added, and
        finally the satellites are sorted by name.
        """
        add_sat = self.popup.avail_sats_lst.currentItem().text()
        xy = ndarray(shape=(200,2), dtype=float)
        self.ax_cov.append(self.ax.add_patch(Polygon(xy, color='w',
                                                     transform=Geodetic(),
                                                     alpha=self.cov_alpha)))
        self.sat_txt.append(self.ax.text([], [], "", color='yellow', size=7,
                            ha="center", va="center"))
        if (add_sat in self.argos):
            self.createSatFromFile(add_sat, "{}argos.txt".format(self.tle_path),
                                   "Argos Data Collection System")
        elif (add_sat in self.beidou):
            self.createSatFromFile(add_sat, "{}beidou.txt".format(self.tle_path),
                                   "Beidou")
        elif (add_sat in self.cubesat):
            self.createSatFromFile(add_sat, "{}cubesat.txt".format(self.tle_path),
                                   "CubeSat")
        elif (add_sat in self.dmc):
            self.createSatFromFile(add_sat, "{}dmc.txt".format(self.tle_path),
                                   "Disaster Monitoring")
        elif (add_sat in self.education):
            self.createSatFromFile(add_sat, "{}education.txt".format(self.tle_path),
                                   "Education")
        elif (add_sat in self.geodetic):
            self.createSatFromFile(add_sat, "{}geodetic.txt".format(self.tle_path),
                                   "Geodetic")
        elif (add_sat in self.goes):
            self.createSatFromFile(add_sat, "{}goes.txt".format(self.tle_path),
                                   "GOES")
        elif (add_sat in self.intelsat):
            self.createSatFromFile(add_sat, "{}intelsat.txt".format(self.tle_path),
                                   "Intelsat")
        elif (add_sat in self.iridium):
            self.createSatFromFile(add_sat, "{}iridium.txt".format(self.tle_path),
                                   "Iridium")
        elif (add_sat in self.iridium_next):
            self.createSatFromFile(add_sat, "{}iridium-NEXT.txt".format(self.tle_path),
                                   "Iridium Next")
        elif (add_sat in self.military):
            self.createSatFromFile(add_sat, "{}military.txt".format(self.tle_path),
                                   "Miscellaneous Military")
        elif (add_sat in self.molniya):
            self.createSatFromFile(add_sat, "{}molniya.txt".format(self.tle_path),
                                   "Molniya")
        elif (add_sat in self.noaa):
            self.createSatFromFile(add_sat, "{}noaa.txt".format(self.tle_path),
                                   "NOAA")
        elif (add_sat in self.oneweb):
            self.createSatFromFile(add_sat, "{}oneweb.txt".format(self.tle_path),
                                   "OneWeb")
        elif (add_sat in self.planet):
            self.createSatFromFile(add_sat, "{}planet.txt".format(self.tle_path),
                                   "Planet")
        elif (add_sat in self.radar):
            self.createSatFromFile(add_sat, "{}radar.txt".format(self.tle_path),
                                   "Radar Calibration")
        elif (add_sat in self.resource):
            self.createSatFromFile(add_sat, "{}resource.txt".format(self.tle_path),
                                   "Earth Resources")
        elif (add_sat in self.sarsat):
            self.createSatFromFile(add_sat, "{}sarsat.txt".format(self.tle_path),
                                   "Search & Rescue")
        elif (add_sat in self.spire):
            self.createSatFromFile(add_sat, "{}spire.txt".format(self.tle_path),
                                   "Spire")
        elif (add_sat in self.starlink):
            self.createSatFromFile(add_sat, "{}starlink.txt".format(self.tle_path),
                                   "Starlink")
        elif (add_sat in self.tdrss):
            self.createSatFromFile(add_sat, "{}tdrss.txt".format(self.tle_path),
                                   "Tracking and Data Relay")
        elif (add_sat in self.tle_new):
            self.createSatFromFile(add_sat, "{}tle-new.txt".format(self.tle_path),
                                   "Last 30 Days' Launches")
        elif (add_sat in self.visual):
            self.createSatFromFile(add_sat, "{}visual.txt".format(self.tle_path),
                                   "100 (or so) Brightest")
        elif (add_sat in self.weather):
            self.createSatFromFile(add_sat, "{}weather.txt".format(self.tle_path),
                                   "Weather")
        elif (add_sat in self.x_comm):
            self.createSatFromFile(add_sat, "{}x-comm.txt".format(self.tle_path),
                                   "Experimental")
        elif (add_sat in self.active):
            self.createSatFromFile(add_sat, "{}active.txt".format(self.tle_path),
                                   "Active Satellites")
        else:
            self.createSatFromFile(add_sat, self.usr_tle_file, "Added by user")
        self.sortSats()
        self.popup.curr_sats_tab.clearContents()
        self.popup.curr_sats_tab.setRowCount(len(self.Sats))
        self.popup.srch_box.setFocus()
        self.sats_lngs.append(0.0)
        self.sats_lats.append(0.0)
        for i, Sat in enumerate(self.Sats):
            self.popup.curr_sats_tab.setItem(i, 0, QtWidgets.QTableWidgetItem(Sat.name))
            self.sat_txt[i].set_text(Sat.name)
            self.sats_lngs[i] = Sat.getLng(date=self.date)
            self.sats_lats[i] = Sat.getLat()

    def createSatFromFile(self, sat_name, file_name, category):
        """
        Creates a new satellite given its name, TLE file and category.

        Parameters
        ----------
        sat_name  : str
                    The name of the new satellite.
        file_name : str
                    The name of the TLE file.
        category  : str
                    The category of the new satellite.
        """
        newSat = Sat(sat_name, tlepath=file_name, cat=category)
        newSat.updateOrbitalParameters(self.date)
        self.Sats.append(newSat)

    def removeSat(self):
        """
        Removes the selected satellite from the list of the current
        satellites (Sats). All the associated data and plots are removed
        as well.
        """
        del_sat = self.popup.curr_sats_tab.currentItem().text()
        for i, sat in enumerate(self.Sats):
            if (del_sat == sat.name):
                self.Sats.remove(sat)
                self.popup.curr_sats_tab.removeRow(i)
                self.ax_cov[i].remove()
                self.ax_cov.remove(self.ax_cov[i])
                self.sat_txt[i].remove()
                self.sat_txt.remove(self.sat_txt[i])
                self.sats_lngs.remove(self.sats_lngs[i])
                self.sats_lats.remove(self.sats_lats[i])
        self.ax_sat.remove()
        self.ax_sat, = self.ax.plot([], [], 'yo', ms=6)
        
    def addRemoveSat(self):
        """
        Creates a dialog for the user to add or remove a satellite.
        """
        Dialog = QtWidgets.QDialog()
        self.popup = Ui_addRemove()
        self.popup.setupUi(Dialog)
        Dialog.setWindowTitle("Add/remove satellites")
        Dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.addRemoveButtons()
        self.showCurrentSats()
        self.showAvailSats()
        self.popup.srch_box.textChanged.connect(self.searchSat)
        Dialog.exec_()

    def fullscreen(self, event=None):
        """
        This method is called when the user uses the F11 key. It calls
        a method to set the fullscreen mode.

        Parameters
        ----------
        event : event
                Key pressed event.
        """
        self.showFullScreen()
        self.ui.actionFullscreen.setText("Exit fullscreen")
        self.ui.actionFullscreen.triggered.connect(self.exitFullscreen)

    def exitFullscreen(self, event=None):
        """
        This method is called when the user presses the ESC key. It
        calls a method to exit the fullscreen mode.

        Parameters
        ----------
        event : event
                Key pressed event.
        """
        self.showMaximized()
        self.ui.actionFullscreen.setText("Fullscreen")
        self.ui.actionFullscreen.triggered.connect(self.fullscreen)

    def updateTLEDialog(self):
        """
        Creates a dialog to give the user the option to update all the
        TLE files from celestrack. The user can see the link of each
        downloaded file in a text box and can also watch the progress
        bar to know when it is finished.
        """
        self.Dialog = QtWidgets.QDialog()
        self.popup = Ui_updateTLE()
        self.popup.setupUi(self.Dialog)
        self.popup.buttonBox.accepted.connect(self.updateTLEfromNet)
        self.popup.buttonBox.rejected.connect(self.Dialog.reject)
        self.Dialog.setWindowTitle("Update TLE from network")
        self.Dialog.exec_()

    def updateTLEfromNet(self):
        """
        Updates the TLE files from the celestrack web page. It also
        handles the display of the current downloaded file and the
        progress bar of the QDialog.
        """
        self.tle_files = ["active", "argos", "beidou", "cubesat", "dmc",
                          "education", "geodetic", "goes", "intelsat",
                          "iridium", "iridium-NEXT", "military",
                          "molniya", "noaa", "oneweb", "planet", "radar",
                          "resource", "sarsat", "spire", "starlink",
                          "tdrss", "tle-new", "visual", "weather", "x-comm"]
        for i, file_name in enumerate(self.tle_files):
            link = "https://celestrak.com/NORAD/elements/{}.txt".format(file_name)
            self.popup.plainTextEdit.insertPlainText("Downloading: {}".format(link))
            with open("{}{}.txt".format(self.tle_path, file_name), "w") as dest:
                response = urlopen(link)
                dest.write(response.read().decode("utf-8"))
            self.popup.plainTextEdit.insertPlainText("\nDone!\n")
            self.popup.progressBar.setValue(int((i + 1)*100/len(self.tle_files)))
        self.popup.plainTextEdit.insertPlainText("Finished!")

    def updateTLEfromFile(self):
        """
        Creates a dialog for the user to add a txt file with TLE data
        of his choice. All the names of the new satellites are added
        into this program and they can be added as any other satellite.
        This data only remains while the program is not closed.
        """
        dlog_name = "Open file"
        file_name = QtWidgets.QFileDialog.getOpenFileName(self, dlog_name, "",
                                                          "Text files (*.txt)")
        self.usr_tle_file = file_name[0]
        self.usr_sats = []
        self.readSatsFromFile(self.usr_tle_file, self.usr_sats)
        self.avail_sats += self.usr_sats
        self.avail_sats.sort()

    def refreshBackgroundImg(self, img=None):
        """
        Updates the background image. By default it uses a picture of
        the Earth in day and night given the position of the sun in the
        sky. It is only done if the World Map tab is selected.

        Parameters
        ----------
        img : str
              Path to the background image.
        """
        tab_index = self.ui.tabWidget.currentIndex()
        if (self.ui.tabWidget.tabText(tab_index) == "World Map"):
            if (img is None):
                self.map.set_data(self.world_map.fillDarkSideFromPicture(self.date))
            else:
                self.map.set_data(imread(img))
            self.canvas.draw_idle()

    def earth(self):
        """
        Sets the world map of the Earth as a background image, and sets
        the Earth parameters on all the current satellites.
        """
        self.refreshBackgroundImg()
        for sat in self.Sats:
            sat.changePlanet()

    def mars(self):
        """
        Sets the world map of Mars as a background image, and sets the
        Mars parameters on all the current satellites.
        """
        self.refreshBackgroundImg("{}mars_nasa_day.png".format(self.img_path))
        for sat in self.Sats:
            sat.changePlanet(M=0.64171*10**24, P_r=3389500, Eq_r=3396200, 
                    Po_r=3376200, J2=0.00196045, P_w=7.08821812*10**(-5))

    def about(self):
        """
        Creates a dialog with information about this program and its
        license.
        """
        Dialog = QtWidgets.QDialog()
        ui = Ui_About()
        ui.setupUi(Dialog)
        Dialog.setWindowTitle('About Pypredict')
        pixmap = QtGui.QPixmap('{}favicon.png'.format(self.img_path))
        ui.icon.setPixmap(pixmap)
        ui.version.setText("Pypredict {}".format(__version__))
        Dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        Dialog.exec_()

    def setMenu(self):
        """
        Connects the menu options with their respective methods.
        """
        self.ui.actionSave_picture.triggered.connect(self.savePicture)
        self.ui.actionSave_TLEs.triggered.connect(self.saveTLEs)
        self.ui.actionUpdate_TLE_from_net.triggered.connect(self.updateTLEDialog)
        self.ui.actionUpdate_TLE_from_file.triggered.connect(self.updateTLEfromFile)
        self.ui.actionEnable_database.triggered.connect(self.enableDB)
        self.ui.actionAdd_south_atlantic_anomaly.triggered.connect(self.addSAA)
        self.ui.actionRemove_coverage.triggered.connect(self.removeCoverage)
        self.ui.actionAdd_remove_satellites.triggered.connect(self.addRemoveSat)
        self.ui.actionFullscreen.triggered.connect(self.fullscreen)
        self.ui.actionEarth.triggered.connect(self.earth)
        self.ui.actionMars.triggered.connect(self.mars)
        self.ui.actionAbout.triggered.connect(self.about)

    def run(self):
        """
        Executes the core methods of this program every certain amount
        of time.
        """
        self.time_timer = QtCore.QTimer()
        self.time_timer.timeout.connect(self.updateTime)
        self.time_timer.start(500)

        self.sats_timer = QtCore.QTimer()
        self.sats_timer.timeout.connect(self.updateSatellites)
        self.sats_timer.start(700)

        self.table_timer = QtCore.QTimer()
        self.table_timer.timeout.connect(self.updateTableContent)
        self.table_timer.start(700)

        self.canvas_timer = QtCore.QTimer()
        self.canvas_timer.timeout.connect(self.updateCanvas)
        self.canvas_timer.start(900)

        self.bg_timer = QtCore.QTimer()
        self.bg_timer.timeout.connect(self.refreshBackgroundImg)
        self.bg_timer.start(60000)
