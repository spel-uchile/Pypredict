"""
                                Pypredict
    Orbit prediction software. Displays the satellites' position and
    orbital parameters in real time. Simulates satellite localization
    and deployment.
    
    Copyright (C) 2018-2019, Matías Vidal Valladares, matvidal.
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
__version__ = "2.3.0"

from PyQt5 import QtWidgets, QtGui, QtCore
from time import sleep
from ui.main_window import Ui_MainWindow
from about.about_window import Ui_About
from deployment.dpl_window import Ui_DPL
from addRemove.addRemove_window import Ui_addRemove
from cartopy.crs import Geodetic, PlateCarree, RotatedPole
from dayNightMap import Map
from dpl import Dpl
from matplotlib.ticker import FixedLocator
from numpy import abs, arange, array, cos, empty, log, pi, sin, tan
from pyorbital import tlefile
from matplotlib.pyplot import imread, subplots, tight_layout
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import datetime, timedelta
from tkinter.filedialog import asksaveasfilename
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from sat import Sat
from SAA import SAA
from warnings import filterwarnings
from pymongo import MongoClient
import ssl
import json

filterwarnings("ignore", category=RuntimeWarning)

class ApplicationWindow(QtWidgets.QMainWindow):
    #@profile
    __slots__ = ["Sats", "img", "mainSat", "mainSat_lats",
                 "mainSat_lngs", "ax_saa", "fig", "ax", "ax_tray",
                 "ax_sat", "ax_cov", "sat_txt", "bg", "fg", "active_bg",
                 "saa_alpha", "cov_alpha", "top_index", "bottom_index",
                 "name_bt", "cat_lbl", "lat_lbl", "lng_lbl", "alt_lbl",
                 "a_lbl", "e_lbl", "raan_lbl", "i_lbl", "w_lbl",
                 "theta_lbl", "h_lbl", "up_bt", "down_bt", "popup",
                 "ent", "match", "avail_sats", "argos", "cubesat", "dmc",
                 "goes", "intelsat", "iridium", "iridium_next", "noaa",
                 "planet", "resource", "sarsat", "spire", "tdrss",
                 "tle_new", "weather", "molniya", "map", "dpl_img",
                 "tdoa_img", "world_map", "dpl", "cov_lat",
                 "cov_lng", "dmin", "canvas", "saa", "date", "db",
                 "en_db", "time_timer", "sats_timer", "canvas_timer",
                 "bg_timer", "Dialog", "table_timer"]

    def __init__(self, Sats):
        self.Sats = Sats
        self.sortSats()
        super(ApplicationWindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        client = MongoClient("localhost", 27017)
        self.db = client["SatConstellation"]
        self.en_db = False
        self.setWindowIcon( QtGui.QIcon("img/favicon.png") )
        self.setWindowTitle('Pypredict')
        self.saa_alpha = 0
        self.cov_alpha = 0.2
        self.saa = SAA()
        self.dmin = 0
        self.ui.datetime.setDateTime(datetime.utcnow())
        self.setButtons()
        self.world_map = Map("img/earth_nasa_day.png",
                             "img/earth_nasa_night.png")
        self.plotData()
        self.setCanvas()
        self.setMenu()
        self.data_gen()
        self.changeMainSat(0, 0)
        self.cov_lng = empty(180)
        self.cov_lat = empty(180)
        #self.ui.Table.setRowCount(len(self.Sats))
        #self.ui.Table.setItem(0,0, QtWidgets.QTableWidgetItem("Item (1,1)"))
        self.updateTableContent()
        self.setTableConnections()
        #self.tableRefresher()
        self.showMaximized()
        self.fig.canvas.mpl_connect('scroll_event',self.zoom)
        self.run()

    def __call__(self):
        return self

    def keyPressEvent(self, e):
        if e.key() == QtCore.Qt.Key_Escape:
            self.exitFullscreen()
        if e.key() == QtCore.Qt.Key_F11:
            if self.isFullScreen():
                self.showMaximized()
            else:
                self.showFullScreen()

    def zoom(self, event):
        lng = event.xdata
        lat = event.ydata
        if (event.button == "up"):
            scale = 0.6
        else:
            scale = 1.4
        lng_min = self.ax.get_xlim()[0]*scale + lng
        lng_max = self.ax.get_xlim()[1]*scale + lng
        lat_min = self.ax.get_ylim()[0]*scale + lat
        lat_max = self.ax.get_ylim()[1]*scale + lat
        if (lng_min < -180 and lng_max > 180 or lat_min < -90 and lat_max > 90):
            self.ax.set_extent([-180, 180, -90, 90], crs=PlateCarree())
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
            self.ax.set_extent([lng_min, lng_max, lat_min, lat_max],
                        crs=PlateCarree())

    def fillSAA(self):
        """
        Plot a fill on the south atlantic anomaly (3 sigma).
        """
        self.ax_saa.set_xy(self.saa.vertices)
        self.ax_saa.set_alpha(self.saa_alpha)

    def plotData(self):
        self.fig, self.ax = subplots(subplot_kw={'projection': PlateCarree()})
        img_extent = (-180, 180, -90, 90)
        self.date = datetime.utcnow()
        self.map = self.ax.imshow(self.world_map.fillDarkSideFromPicture(self.date),
                origin='upper', extent=img_extent, transform=PlateCarree())
        self.gridAndFormat()
        tight_layout(pad=-0.26)

    def gridAndFormat(self):
        gl = self.ax.gridlines(crs=PlateCarree(), draw_labels=True,
                               linewidth=1, color='gray', alpha=0.5, linestyle='-')
        gl.xlabels_top = True
        gl.xlabels_bottom = False
        gl.ylabels_right = False
        gl.ylabels_left = True
        gl.xlines = True
        gl.xlocator = FixedLocator([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
        gl.ylocator = FixedLocator([-90, -60, -30, 0, 30, 60, 90])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 9, 'color': 'white'}
        gl.ylabel_style = {'size': 9, 'color': 'white'}
        gl.xpadding = -10
        gl.ypadding = -25
        return gl

    def data_gen(self):
        self.ax_saa, = self.ax.fill([0,0], [0,0], color='green',
                                          alpha=self.saa_alpha)
        self.ax_tray, = self.ax.plot([], [], color='green',
                                     linewidth=1.5, linestyle='--',
                                     transform=Geodetic())
        self.ax_sat, = self.ax.plot([], [], 'yo', ms=6)
        self.ax_cov = []
        self.sat_txt = []
        self.resetCov()

    def resetCov(self):
        for i in range(0, len(self.sat_txt)):
            self.ax_cov[i].remove()
            self.sat_txt[i].remove()
        self.ax_cov = []
        self.sat_txt = []
        for i in range(0, len(self.Sats)):
            self.ax_cov.append(self.ax.fill([0,0], [0,0], transform=Geodetic(),
                               color='white', alpha=self.cov_alpha)[0])
            self.sat_txt.append(self.ax.text([], [], "", color='yellow', size=8,
                                transform=Geodetic(), ha="center"))

    def updateCanvas(self):
        self.fillSAA() 
        sats_lngs = [] 
        sats_lats = [] 
        for i, Sat in  enumerate(self.Sats):
            sats_lngs. append(Sat.getLng(date=self.date))
            lng = sats_lngs[i] - 5*(sats_lngs[i] > 174) + 5*(sats_lngs[i] < -174)
            sats_lats.append(Sat.getLat())
            lat = sats_lats[i] - (1 - 2*(sats_lats[i] < -85))*4
            self.sat_txt[i].set_position(array((lng, lat)))
            self.sat_txt[i].set_text(Sat.name)
            self.plotCoverage(Sat.getCoverage(), sats_lats[i], sats_lngs[i], i)
        self.ax_sat.set_data(sats_lngs, sats_lats)
        self.canvas.draw_idle()
        self.ui.datetime.setDateTime(self.date)
        #self.update()

    def plotCoverage(self, ang, sat_lat, sat_lng, n, deg2rad=pi/180, rad2deg=180/pi):
        for i in range(0, 180):
            theta = (2*i + 1)*deg2rad
            dlat = ang*deg2rad * cos(theta)
            self.cov_lat[179-i] = sat_lat*deg2rad + dlat
            dpsi = log(tan(self.cov_lat[179-i]*0.5 + pi*0.25)/tan(sat_lat*deg2rad*0.5 + pi*0.25))
            if (abs(dpsi) > 10e-12):
                q = dlat / dpsi
            else:
                q = cos(sat_lat*deg2rad)
            dlng = ang*deg2rad*sin(theta)/q
            self.cov_lng[179-i] = sat_lng*deg2rad + dlng
            if (abs(self.cov_lat[179-i]) > (6*deg2rad + pi - ang*deg2rad - abs(sat_lat*deg2rad))):
                self.cov_lat[179-i] = (2*(sat_lat > 0) - 1)*(6*deg2rad + pi - ang*deg2rad - abs(sat_lat*deg2rad))
                self.cov_lng[179-i] = sat_lng*deg2rad - ((sat_lng*deg2rad - pi) > 0)*pi + ((sat_lng*deg2rad - pi) <= 0)*pi
        self.cov_lat = self.cov_lat*rad2deg
        self.cov_lng = self.cov_lng*rad2deg
        self.ax_cov[n].set_alpha(self.cov_alpha) 
        self.ax_cov[n].set_xy(array((self.cov_lng, self.cov_lat)).transpose())

    def setTheme(self, bg, fg, active_bg):
        self.bg = bg
        self.fg = fg
        self.active_bg = active_bg
        self.root.configure(background=bg)
        style = ttk.Style()
        style.configure("BW.TLabel", foreground=self.fg, background=self.bg)
        style.map("BW.TLabel", background=[("active", self.active_bg)])
        style.configure('play.TLabel', font=('TkDefaultFont', 16, 'bold'),
                    foreground=self.fg, background=self.bg, anchor='center')
        style.map("play.TLabel", background=[("active", self.active_bg)])
        style.configure('nextPrev.TLabel', font=('TkDefaultFont', 10, 'bold'),
                    foreground=self.fg, background=self.bg, anchor='center')
        style.map("nextPrev.TLabel", background=[("active", self.active_bg)])

    def changeMainSat(self, row, col):
        self.mainSat = self.Sats[row]
        tf = int(self.mainSat.getPeriod()*3)
        self.mainSat_lats, self.mainSat_lngs = self.mainSat.getTrayectory(tf, 50, self.date)
        self.ax_tray.set_data(self.mainSat_lngs, self.mainSat_lats)
   
    def up(self):
        if (self.top_index > 0):
            self.top_index -= 1
            self.bottom_index -= 1

    def down(self):
        if (self.bottom_index < len(self.Sats)):
            self.top_index += 1
            self.bottom_index += 1

    def setButtons(self):
        self.ui.dpl_button.clicked.connect(self.deployPopup)
        self.ui.loc_button.clicked.connect(self.notAvailable)
        self.ui.prev_day.clicked.connect(self.previousDay)
        self.ui.prev_min.clicked.connect(self.previousMinute)
        self.ui.play.clicked.connect(self.currentMinute)
        self.ui.next_min.clicked.connect(self.nextMinute)
        self.ui.next_day.clicked.connect(self.nextDay)
        self.ui.datetime.dateTimeChanged.connect(self.newDate)

    def deployPopup(self):
        self.Dialog = QtWidgets.QDialog()
        self.popup = Ui_DPL()
        self.popup.setupUi(self.Dialog)
        self.Dialog.setWindowTitle("Deployment settings")
        self.showCurrentSats()
        self.popup.curr_sats_lst.itemClicked.connect(self.selectDeployer)
        self.popup.deploy_now_bt.clicked.connect(self.deploySat)
        self.dpl = Dpl()
        self.Dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.Dialog.exec_()

    def selectDeployer(self, item):
        deployer_name = item.text()
        self.popup.deployer_name_lbl.setText(deployer_name)

    def deploySat(self):
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
        self.ax_cov.append(self.ax.fill([0,0], [0,0], transform=Geodetic(),
                           color='white', alpha=self.cov_alpha)[0])
        self.sat_txt.append(self.ax.text([], [], "", color='yellow', size=8,
                            transform=Geodetic(), ha="center"))
        self.Sats.append(newSat)
        self.sortSats()
        self.Dialog.close()

    def previousDay(self):
        self.dmin -= 1440
        self.format_dt()

    def previousMinute(self):
        self.dmin -= 1
        self.format_dt()

    def currentMinute(self):
        self.dmin = 0
        self.format_dt()

    def nextMinute(self):
        self.dmin += 1
        self.format_dt()

    def nextDay(self):
        self.dmin += 1440
        self.format_dt()

    def format_dt(self):
        self.updateTime()
        self.ui.datetime.setDateTime(self.date)
        self.refreshBackgroundImg()

    def updateTime(self, date=None):
        if (date is None):
            self.date = datetime.utcnow() + timedelta(minutes=self.dmin)
        else:
            self.date = date

    def newDate(self):
        date = self.ui.datetime.dateTime().toPyDateTime()
        self.dmin += (date - self.date).total_seconds()/60.0
        self.date = date
        self.refreshBackgroundImg()

    def updateSatellites(self):
        for Sat in self.Sats:
            Sat.updateOrbitalParameters3(self.date)
            if (self.en_db):
                col = self.db[Sat.name]
                col.insert_one(self.formatDump(Sat))
        

    def setCanvas(self):
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setParent(self)
        self.canvas.draw()
        self.ui.frame_plot.layout().addWidget(self.canvas)

    def setTableConnections(self):
        self.ui.Table.cellClicked.connect(self.changeMainSat)

    def updateTableContent(self):
        rad2deg = 180/pi
        self.ui.Table.setRowCount(len(self.Sats))
        #self.ui.Table.setColumnWidth(0, 50)
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

    def rememberLastRows(self):
        if (len(self.Sats) > 4):
            self.top_index = 0
            self.bottom_index = 5
            self.rememberRow(self.bottom_index)
            if (len(self.Sats) >= 6):
                self.bottom_index = 6
                self.rememberRow(self.bottom_index)
                self.up_bt.grid(row=5, column=17, rowspan=3, sticky="NESW")
                self.down_bt.grid(row=8, column=17, rowspan=3, sticky="NESW")

    def rememberRow(self, r):
        self.root.rowconfigure(r+4, weight=1)
        self.name_bt[r-1].grid(row=r+4, column=0, columnspan=5, sticky="EW")
        self.cat_lbl[r-1].grid(row=r+4, column=5, sticky="W")
        self.lat_lbl[r-1].grid(row=r+4, column=6)
        self.lng_lbl[r-1].grid(row=r+4, column=7)
        self.alt_lbl[r-1].grid(row=r+4, column=8)
        self.spd_lbl[r-1].grid(row=r+4, column=9)
        self.a_lbl[r-1].grid(row=r+4, column=10)
        self.h_lbl[r-1].grid(row=r+4, column=11)
        self.e_lbl[r-1].grid(row=r+4, column=12)
        self.raan_lbl[r-1].grid(row=r+4, column=13)
        self.i_lbl[r-1].grid(row=r+4, column=14)
        self.w_lbl[r-1].grid(row=r+4, column=15)
        self.theta_lbl[r-1].grid(row=r+4, column=16)

    def forgetLastRows(self):
        self.top_index = 0
        self.bottom_index = 4*(len(self.Sats) > 4) + len(self.Sats)*(len(self.Sats) <= 4)
        self.up_bt.grid(row=5, column=17, rowspan=2, sticky="NESW")
        self.down_bt.grid(row=7, column=17, rowspan=2, sticky="NESW")
        if (len(self.Sats) > 4):
            self.forgetRow(self.bottom_index)
            if (len(self.Sats) >= 6):
                self.forgetRow(self.bottom_index + 1)

    def forgetRow(self, r):
        self.name_bt[r].grid_forget()
        self.cat_lbl[r].grid_forget()
        self.lat_lbl[r].grid_forget()
        self.lng_lbl[r].grid_forget()
        self.alt_lbl[r].grid_forget()
        self.spd_lbl[r].grid_forget()
        self.a_lbl[r].grid_forget()
        self.h_lbl[r].grid_forget()
        self.e_lbl[r].grid_forget()
        self.raan_lbl[r].grid_forget()
        self.i_lbl[r].grid_forget()
        self.w_lbl[r].grid_forget()
        self.theta_lbl[r].grid_forget()

    def tableRefresher(self):
        self.updateTableContent()
        self.root.after(500, self.tableRefresher)

    def formatDump(self, Sat):
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

    def saveAs(self):
        file_name = asksaveasfilename()
        if file_name is None:
            return
        self.fig.savefig(file_name)

    def notAvailable(self):
        print("This command is not available yet")

    def enableDB(self):
        self.en_db = True
        self.ui.actionEnable_database.setText("Disable database")
        self.ui.actionEnable_database.triggered.connect(self.disableDB)
        #self.editmenu.entryconfigure(2, label="Disable database",
        #                             command=self.disableDB)

    def disableDB(self):
        self.en_db = False
        self.ui.actionEnable_database.setText("Enable database")
        self.ui.actionEnable_database.triggered.connect(self.enableDB)
        #self.editmenu.entryconfigure(2, label="Enable database",
        #                             command=self.enableDB)

    def removeSAA(self):
        self.saa_alpha = 0
        self.updateCanvas()
        self.ui.actionAdd_south_atlantic_anomaly.setText("Add south atlantic anomaly")
        self.ui.actionAdd_south_atlantic_anomaly.triggered.connect(self.addSAA)
        #self.editmenu.entryconfigure(4, label="Add south atlantic anomaly",
        #                             command=self.addSAA)

    def addSAA(self):
        self.saa_alpha = 0.2
        self.updateCanvas()
        self.ui.actionAdd_south_atlantic_anomaly.setText("Remove south atlantic anomaly")
        self.ui.actionAdd_south_atlantic_anomaly.triggered.connect(self.removeSAA)
        #self.editmenu.entryconfigure(4, label="Remove south atlantic anomaly",
        #                             command=self.removeSAA)

    def removeCoverage(self):
        self.cov_alpha = 0
        self.updateCanvas()
        self.ui.actionRemove_coverage.setText("Add coverage")
        self.ui.actionRemove_coverage.triggered.connect(self.addCoverage)
        #self.editmenu.entryconfigure(5, label="Add coverage",
        #                             command=self.addCoverage)

    def addCoverage(self):
        self.cov_alpha = 0.2
        self.updateCanvas()
        self.ui.actionRemove_coverage.setText("Remove coverage")
        self.ui.actionRemove_coverage.triggered.connect(self.removeCoverage)
        #self.editmenu.entryconfigure(5, label="Remove coverage", 
        #                             command=self.removeCoverage)

    def setSearchBox(self):
        self.srch_lbl = Label(self.popup, text="Search:")
        self.srch_lbl.grid(row=0, column=0)
        self.ent = StringVar()
        self.srch_box = Entry(self.popup, textvariable=self.ent)
        self.srch_box.grid(row=0, column=1)
        self.srch_box.focus()

    #def searchSat(self, event):
    def searchSat(self, text):
        srch = text.upper()#self.srch_box.get().upper()
        self.match = [s for s in self.avail_sats if srch in s]
        self.match.sort(reverse=False) 
        for i in range(0, self.popup.avail_sats_lst.count()):
            #self.avail_sats_lst.delete(0)
            self.popup.avail_sats_lst.clear()#takeItem(i)
        for sat in self.match:
            self.popup.avail_sats_lst.addItem(sat)

    def readSatsFromFile(self, file, lst):
        with open(file, 'r') as f:
            for count, line in enumerate(f):
                if (count % 3 == 0):
                    lst.append(line.strip())

    def readAllSats(self):
        self.argos = []
        self.cubesat = []
        self.dmc = []
        self.goes = []
        self.intelsat = []
        self.iridium = []
        self.iridium_next = []
        self.molniya = []
        self.noaa = []
        self.planet = []
        self.resource = []
        self.sarsat = []
        self.spire = []
        self.tdrss = []
        self.tle_new = []
        self.weather = []
        self.readSatsFromFile("TLE/argos.txt", self.argos)
        self.readSatsFromFile("TLE/cubesat.txt", self.cubesat)
        self.readSatsFromFile("TLE/dmc.txt", self.dmc)
        self.readSatsFromFile("TLE/goes.txt", self.goes)
        self.readSatsFromFile("TLE/intelsat.txt", self.intelsat)
        self.readSatsFromFile("TLE/iridium.txt", self.iridium)
        self.readSatsFromFile("TLE/iridium-NEXT.txt", self.iridium_next)
        self.readSatsFromFile("TLE/molniya.txt", self.molniya)
        self.readSatsFromFile("TLE/noaa.txt", self.noaa)
        self.readSatsFromFile("TLE/planet.txt", self.planet)
        self.readSatsFromFile("TLE/resource.txt", self.resource)
        self.readSatsFromFile("TLE/sarsat.txt", self.sarsat)
        self.readSatsFromFile("TLE/spire.txt", self.spire)
        self.readSatsFromFile("TLE/tdrss.txt", self.tdrss)
        self.readSatsFromFile("TLE/tle-new.txt", self.tle_new)
        self.readSatsFromFile("TLE/weather.txt", self.weather)
        self.avail_sats = self.argos + self.cubesat + self.dmc + self.goes
        self.avail_sats += self.intelsat + self.iridium + self.iridium_next
        self.avail_sats += self.molniya + self.noaa + self.planet
        self.avail_sats += self.resource + self.sarsat + self.spire
        self.avail_sats += self.tdrss + self.tle_new + self.weather
        self.avail_sats.sort()

    def showAvailSats(self):
        #self.avail_sats_lst = Listbox(self.popup, width=28)
        self.readAllSats()
        for sat in self.avail_sats:
            self.popup.avail_sats_lst.addItem(sat)
        #self.avail_sats_lst.grid(row=1, column=0, rowspan=2, columnspan=2)
        #self.popup.columnconfigure(0, weight=1)

    def showCurrentSats(self):
        for sat in self.Sats:
            self.popup.curr_sats_lst.addItem(sat.name)

    def addRemoveButtons(self):
        self.popup.add_sat_bt.clicked.connect(self.addSat)
        self.popup.remove_sat_bt.clicked.connect(self.removeSat)
        '''
        self.add_sat_bt = Button(self.popup, text="→",
                                 command=self.addSat)
        self.add_sat_bt.grid(row=1, column=2)
        self.remove_sat_bt = Button(self.popup, text="←",
                                    command=self.removeSat)
        self.remove_sat_bt.grid(row=2, column=2)
        '''

    def sortSats(self):
        self.Sats.sort(key = lambda s: s.name)

    def addSat(self):
        add_sat = self.popup.avail_sats_lst.currentItem().text()#self.popup.avail_sats_lst.get(self.popup.avail_sats_lst.curselection())
        self.ax_cov.append(self.ax.fill([0,0], [0,0], transform=Geodetic(),
                           color='white', alpha=self.cov_alpha)[0])
        self.sat_txt.append(self.ax.text([], [], "", color='yellow', size=8,
                            transform=Geodetic(), ha="center"))
        if (add_sat in self.argos):
            self.createSatFromFile(add_sat, "TLE/argos.txt", "Argos Data Collection System")
        elif (add_sat in self.cubesat):
            self.createSatFromFile(add_sat, "TLE/cubesat.txt", "CubeSat")
        elif (add_sat in self.dmc):
            self.createSatFromFile(add_sat, "TLE/dmc.txt", "Disaster Monitoring")
        elif (add_sat in self.goes):
            self.createSatFromFile(add_sat, "TLE/goes.txt", "GOES")
        elif (add_sat in self.intelsat):
            self.createSatFromFile(add_sat, "TLE/intelsat.txt", "Intelsat")
        elif (add_sat in self.iridium):
            self.createSatFromFile(add_sat, "TLE/iridium.txt", "Iridium")
        elif (add_sat in self.iridium_next):
            self.createSatFromFile(add_sat, "TLE/iridium-NEXT.txt", "Iridium Next")
        elif (add_sat in self.molniya):
            self.createSatFromFile(add_sat, "TLE/molniya.txt", "Molniya")
        elif (add_sat in self.noaa):
            self.createSatFromFile(add_sat, "TLE/noaa.txt", "NOAA")
        elif (add_sat in self.planet):
            self.createSatFromFile(add_sat, "TLE/planet.txt", "Planet")
        elif (add_sat in self.resource):
            self.createSatFromFile(add_sat, "TLE/resource.txt", "Earth Resources")
        elif (add_sat in self.sarsat):
            self.createSatFromFile(add_sat, "TLE/sarsat.txt", "Search & Rescue")
        elif (add_sat in self.spire):
            self.createSatFromFile(add_sat, "TLE/spire.txt", "Spire")
        elif (add_sat in self.tdrss):
            self.createSatFromFile(add_sat, "TLE/tdrss.txt", "Tracking and Data Relay")
        elif (add_sat in self.tle_new):
            self.createSatFromFile(add_sat, "TLE/tle-new.txt", "Last 30 Days' Launches")
        elif (add_sat in self.weather):
            self.createSatFromFile(add_sat, "TLE/weather.txt", "Weather")
        else:
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat)))
        self.popup.curr_sats_lst.addItem(add_sat)#insert(END, add_sat)
        self.sortSats()
        self.popup.curr_sats_lst.clear()
        self.popup.srch_box.setFocus()
        for i, sat in enumerate(self.Sats):
            self.popup.curr_sats_lst.addItem(sat.name)#insert(i, sat.name)

    def createSatFromFile(self, sat_name, file_name, category):
        newSat = Sat(sat_name, tle=tlefile.read(sat_name, file_name), cat=category)
        newSat.updateOrbitalParameters3(self.date)
        self.Sats.append(newSat)

    def removeSat(self):
        del_sat = self.popup.curr_sats_lst.currentItem().text()
        for i, sat in enumerate(self.Sats):
            if (del_sat == sat.name):
                self.Sats.remove(sat)
                self.popup.curr_sats_lst.takeItem(i)
        self.resetCov()
        
    def addRemoveSat(self):
        '''
        self.popup = Tk()
        self.popup.title("Add/remove satellites")
        #self.popup.call('wm', 'iconphoto', self.popup._w, self.img)
        self.showAvailSats()
        self.setSearchBox()
        self.addRemoveButtons()
        self.showCurrentSats(3)
        self.popup.bind("<Key>", self.searchSat)
        self.popup.protocol("WM_DELETE_WINDOW", self.popup.destroy)
        self.popup.mainloop()
        '''
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
        self.showFullScreen()
        self.ui.actionFullscreen.setText("Exit fullscreen")
        self.ui.actionFullscreen.triggered.connect(self.exitFullscreen)

    def exitFullscreen(self, event=None):
        self.showMaximized()
        self.ui.actionFullscreen.setText("Fullscreen")
        self.ui.actionFullscreen.triggered.connect(self.fullscreen)

    def updateTLEfromNet(self):
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/argos.txt", )
        tlefile.fetch("TLE/argos.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/cubesat.txt", )
        tlefile.fetch("TLE/cubesat.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/dmc.txt", )
        tlefile.fetch("TLE/dmc.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/goes.txt", )
        tlefile.fetch("TLE/goes.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/intelsat.txt", )
        tlefile.fetch("TLE/intelsat.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/iridium.txt", )
        tlefile.fetch("TLE/iridium.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/iridium-NEXT.txt", )
        tlefile.fetch("TLE/iridium-NEXT.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/molniya.txt", )
        tlefile.fetch("TLE/molniya.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/noaa.txt", )
        tlefile.fetch("TLE/noaa.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/planet.txt", )
        tlefile.fetch("TLE/planet.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/resource.txt", )
        tlefile.fetch("TLE/resource.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/sarsat.txt", )
        tlefile.fetch("TLE/sarsat.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/spire.txt", )
        tlefile.fetch("TLE/spire.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/tdrss.txt", )
        tlefile.fetch("TLE/tdrss.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/tle-new.txt", )
        tlefile.fetch("TLE/tle-new.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/visual.txt", )
        tlefile.fetch("TLE/visual.txt")
        tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/weather.txt", )
        tlefile.fetch("TLE/weather.txt")

    def refreshBackgroundImg(self, img=None):
        if (img is None):
            self.map.set_data(self.world_map.fillDarkSideFromPicture(self.date))
        else:
            self.map.set_data(imread(img))
        self.canvas.draw_idle()
        self.fig.canvas.draw_idle()

    def earth(self):
        self.refreshBackgroundImg()
        mu = 5.9722*6.67408*10**13
        for sat in self.Sats:
            sat.changePlanet()

    def mars(self):
        self.refreshBackgroundImg("img/mars_nasa_day.png")
        mu = 0.64171*6.67408*10**13
        for sat in self.Sats:
            sat.changePlanet(M=0.64171*10**24, P_r=3389500, Eq_r=3396200, 
                    Po_r=3376200, J2=0.00196045, P_w=7.08821812*10**(-5))

    def about(self):
        Dialog = QtWidgets.QDialog()
        ui = Ui_About()
        ui.setupUi(Dialog)
        Dialog.setWindowTitle('About Pypredict')
        pixmap = QtGui.QPixmap('img/favicon.png')
        ui.icon.setPixmap(pixmap)
        ui.version.setText("Pypredict {}".format(__version__))
        Dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        Dialog.exec_()

    def setMenu(self):
        self.ui.actionSave_as.triggered.connect(self.saveAs)
        self.ui.actionUpdate_TLE_from_net.triggered.connect(self.updateTLEfromNet)
        self.ui.actionEnable_database.triggered.connect(self.enableDB)
        self.ui.actionAdd_south_atlantic_anomaly.triggered.connect(self.addSAA)
        self.ui.actionRemove_coverage.triggered.connect(self.removeCoverage)
        self.ui.actionAdd_remove_satellites.triggered.connect(self.addRemoveSat)
        self.ui.actionFullscreen.triggered.connect(self.fullscreen)
        self.ui.actionEarth.triggered.connect(self.earth)
        self.ui.actionMars.triggered.connect(self.mars)
        self.ui.actionAbout.triggered.connect(self.about)

    def run(self):
        self.time_timer = QtCore.QTimer()
        self.time_timer.timeout.connect(self.updateTime)
        self.time_timer.start(100)

        self.sats_timer = QtCore.QTimer()
        self.sats_timer.timeout.connect(self.updateSatellites)
        self.sats_timer.start(300)

        self.table_timer = QtCore.QTimer()
        self.table_timer.timeout.connect(self.updateTableContent)
        self.table_timer.start(500)

        self.canvas_timer = QtCore.QTimer()
        self.canvas_timer.timeout.connect(self.updateCanvas)
        self.canvas_timer.start(500)

        self.bg_timer = QtCore.QTimer()
        self.bg_timer.timeout.connect(self.refreshBackgroundImg)
        self.bg_timer.start(60000)
