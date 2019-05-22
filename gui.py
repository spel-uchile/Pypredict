'''
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
'''
__version__ = "2.2.0"

from cartopy.crs import Geodetic, PlateCarree, RotatedPole
from dayNightMap import Map
from dpl import Dpl
from matplotlib.ticker import FixedLocator
from numpy import abs, arange, array, cos, empty, log, pi, sin, tan
from matplotlib.animation import FuncAnimation
from pyorbital import tlefile
from matplotlib.pyplot import imread, subplots, tight_layout
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import datetime, timedelta
from tkinter import Button, Entry, END, Image, Label, Listbox, Menu, StringVar, ttk, Tk
from tkinter.filedialog import asksaveasfilename
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from sat import Sat
from warnings import filterwarnings
import ssl

filterwarnings("ignore", category=RuntimeWarning)
ssl._create_default_https_context = ssl._create_unverified_context

class GUI(object):
    #@profile
    __slots__ = ["Sats", "root", "img", "mainSat", "dt", "mainSat_lats",
                 "mainSat_lngs", "ax_dark_side", "fig", "ax", "ax_tray",
                 "ax_sat", "ax_cov", "sat_txt", "bg", "fg", "active_bg",
                 "night_alpha", "cov_alpha", "top_index", "bottom_index",
                 "name_bt", "cat_lbl", "lat_lbl", "lng_lbl", "alt_lbl",
                 "a_lbl", "e_lbl", "raan_lbl", "i_lbl", "w_lbl",
                 "theta_lbl", "h_lbl", "up_bt", "down_bt", "popup",
                 "srch_lbl", "ent", "srch_box", "match", "avail_sats",
                 "argos", "cubesat", "dmc", "goes", "intelsat",
                 "iridium", "iridium_next", "noaa", "planet", "resource",
                 "sarsat", "spire", "tdrss", "tle_new", "weather",
                 "avail_sats_lst", "curr_lbl", "curr_sats_lst",
                 "add_sat_bt", "remove_sat_bt", "editmenu", "viewmenu",
                 "planetmenu", "map", "ani", "spd_lbl", "dpl_bt",
                 "deploy_now_bt", "loc_bt", "dpl_img", "tdoa_img",
                 "world_map", "dpl_mass_lbl1", "dpl_mass_lbl2",
                 "dpl_spdx_lbl", "dpl_spdy_lbl", "dpl_spdz_lbl",
                 "mass_box1", "mass_box2", "spdx_box", "spdy_box",
                 "spdz_box", "dpl_name_lbl", "name_box", "dpl_cat_lbl",
                 "cat_box", "deployer_lbl", "deployer_name_lbl", "dpl",
                 "prog_name_lbl", "version_lbl", "dev_lbl", "contact_lbl",
                 "updtCnt", "warranty_lbl", "details_lbl", "cov_lat",
                 "cov_lng", "play", "next_min", "next_day", "prev_min",
                 "prev_day", "dmin", "dt_lbl", "molniya"]
    def __init__(self, Sats):
        self.Sats = Sats
        self.sortSats()
        self.root = Tk()
        self.img = Image("photo", file="img/favicon.png") 
        self.root.call('wm', 'iconphoto', self.root._w, self.img)
        self.geometry()
        self.root.title('Pypredict')
        self.setTheme('#404040', 'white', '#303030')
        self.root.protocol("WM_DELETE_WINDOW", exit)
        self.mainSat = self.Sats[0]
        T = int(self.mainSat.getPeriod()*3)      # Total duration in seconds
        self.dt = 1000                           # Step length in miliseconds
        self.updtCnt = 0
        self.dmin = 0
        #self.mainSat_lats, self.mainSat_lngs = self.mainSat.getLocation(T, 50)
        self.mainSat_lats, self.mainSat_lngs = self.mainSat.getTrayectory(T, 50)
        self.world_map = Map("img/earth_nasa_day.png", "img/earth_nasa_night.png")
        self.plotData()
        self.setButtons()
        self.setCanvas()
        self.setMenu()
        self.setTableTitles()
        self.setTableContent()
        self.tableRefresher()
        self.cov_lng = empty(180)
        self.cov_lat = empty(180)
        self.ani = FuncAnimation(self.fig, self.update, self.data_gen(),
                            interval=self.dt, blit=True, repeat=False)
        self.root.bind("<F11>", self.fullscreen)
        self.root.bind("<Escape>", self.exitFullscreen)
        self.run()

    def __call__(self):
        return self

    def sun_pos(self, dt=None):
        """This function computes a rough estimate of the coordinates for
        the point on the surface of the Earth where the Sun is directly
        overhead at the time dt. Precision is down to a few degrees. This
        means that the equinoxes (when the sign of the latitude changes)
        will be off by a few days.

        The function is intended only for visualization. For more precise
        calculations consider for example the PyEphem package.

        Parameters
        ----------
            dt: datetime
                Defaults to datetime.utcnow()

        Returns
        -------
            lat, lng: tuple of floats
                Approximate coordinates of the point where the sun is
                in zenith at the time dt.
        """
        if dt is None:
            dt = datetime.utcnow()

        axial_tilt = 23.4
        ref_solstice = datetime(2016, 6, 21, 22, 22)
        days_per_year = 365.2425
        seconds_per_day = 24*3600.0

        days_since_ref = (dt - ref_solstice).total_seconds()/seconds_per_day
        lat = axial_tilt*cos(2*pi*days_since_ref/days_per_year)
        sec_since_midnight = (dt - datetime(dt.year, dt.month, dt.day)).seconds
        lng = -(sec_since_midnight/seconds_per_day - 0.5)*360
        return lat, lng

    def fill_dark_side(self, time=None):
        """
        Plot a fill on the dark side of the planet (without refraction).

        Parameters
        ----------
            ax : matplotlib axes
                The axes to plot on.
            time : datetime
                The time to calculate terminator for. Defaults to datetime.utcnow()
            **kwargs :
                Passed on to Matplotlib's ax.fill()

        """
        lat, lng = self.sun_pos(time)
        pole_lng = lng
        if lat > 0:
            pole_lat = -90 + lat
            central_rot_lng = 180
        else:
            pole_lat = 90 + lat
            central_rot_lng = 0

        rotated_pole = RotatedPole(pole_latitude=pole_lat,
                                   pole_longitude=pole_lng,
                                   central_rotated_longitude=central_rot_lng)

        x = empty(360)
        y = empty(360)
        x[:180] = -90
        y[:180] = arange(-90, 90.)
        x[180:] = 90
        y[180:] = arange(90, -90., -1)

        self.ax_dark_side.set_xy(array((x,y)).transpose())
        self.ax_dark_side.set_transform(rotated_pole)
        self.ax_dark_side.set_alpha(self.night_alpha)

    def plotData(self):
        screen_height = self.root.winfo_screenheight()
        if (screen_height < 1080):
            self.fig, self.ax = subplots(figsize=(12, 6),
                                    subplot_kw={'projection': PlateCarree()})
        elif (screen_height == 1080):
            self.fig, self.ax = subplots(figsize=(17.5, 8.75),
                                    subplot_kw={'projection': PlateCarree()})
        else:
            self.fig, self.ax = subplots(figsize=(18, 9),
                                    subplot_kw={'projection': PlateCarree()})
        img_extent = (-180, 180, -90, 90)
        #self.map = self.ax.imshow(imread("img/earth_nasa_day.png"), origin='upper',
        #               extent=img_extent, transform=PlateCarree())
        date = datetime.utcnow() + timedelta(minutes=self.dmin)
        self.map = self.ax.imshow(self.world_map.fillDarkSideFromPicture(date),
                origin='upper', extent=img_extent, transform=PlateCarree())
        self.gridAndFormat()
        tight_layout(pad=-0.26)

    def gridAndFormat(self):
        gl = self.ax.gridlines(crs=PlateCarree(), draw_labels=True,
                               linewidth=1, color='gray', alpha=0.5, linestyle='-')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.ylabels_left = True
        gl.xlines = True
        gl.xlocator = FixedLocator([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
        gl.ylocator = FixedLocator([-90, -60, -30, 0, 30, 60, 90])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 9, 'color': 'black'}
        gl.ylabel_style = {'size': 9, 'color': 'white'}
        gl.xpadding = -10
        gl.ypadding = -25
        return gl

    def data_gen(self):
        self.ax_dark_side, = self.ax.fill([0,0], [0,0], color='black',
                                          alpha=self.night_alpha)
        self.ax_tray, = self.ax.plot([], [], color='green',
                                     linewidth=1.5, linestyle='--',
                                     transform=Geodetic())
        self.ax_sat, = self.ax.plot([], [], 'yo', ms=6)
        self.ax_cov = []
        self.sat_txt = []
        sats_lngs = []
        sats_lats = []
        sats_angs = []
        sats_names = []
        for i in range(0, len(self.Sats)):
            self.ax_cov.append(self.ax.fill([0,0], [0,0], transform=Geodetic(),
                               color='white', alpha=self.cov_alpha)[0])
            self.sat_txt.append(self.ax.text([], [], "", color='yellow', size=8,
                                        transform=Geodetic(), ha="center"))
        while True:
            sats_lngs[:] = []
            sats_lats[:] = []
            sats_angs[:] = []
            sats_names[:] = []
            date = datetime.utcnow() + timedelta(minutes=self.dmin)
            for Sat in self.Sats:
                sats_lngs.append(Sat.getLng(date=date))
                sats_lats.append(Sat.getLat())
                sats_angs.append(Sat.getCoverage())
                sats_names.append(Sat.name)
            yield sats_lngs, sats_lats, sats_angs, sats_names,
    #@profile
    def update(self, data):
        sats_lngs, sats_lats, sats_angs, sats_names, = data
        self.fill_dark_side()
        self.ax_tray.set_data(self.mainSat_lngs, self.mainSat_lats)
        self.ax_sat.set_data(sats_lngs, sats_lats)
        for i in range(0, len(self.Sats)):
            lng = sats_lngs[i] - 5*(sats_lngs[i] > 174) + 5*(sats_lngs[i] < -174)
            lat = sats_lats[i]-(1 - 2*(sats_lats[i] < -85))*4
            self.sat_txt[i].set_position(array((lng, lat)))
            self.sat_txt[i].set_text(sats_names[i])
            self.plotCoverage(sats_angs[i], sats_lats[i], sats_lngs[i], i)
        return [self.ax_dark_side] + [self.ax_tray] + [self.ax_sat] + self.sat_txt + self.ax_cov

    def plotCoverage(self, ang, sat_lat, sat_lng, n):
        deg2rad = pi/180
        rad2deg = 180/pi
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
        self.night_alpha = 0#0.7
        self.cov_alpha = 0.2

    def changeMainSat(self, Sat):
        self.mainSat = Sat
        tf = int(self.mainSat.getPeriod()*3)
        #self.mainSat_lats, self.mainSat_lngs = self.mainSat.getLocation(tf, 50, self.dmin)
        self.mainSat_lats, self.mainSat_lngs = self.mainSat.getTrayectory(tf, 50, self.dmin)
   
    def up(self):
        if (self.top_index > 0):
            self.top_index -= 1
            self.bottom_index -= 1

    def down(self):
        if (self.bottom_index < len(self.Sats)):
            self.top_index += 1
            self.bottom_index += 1

    def setButtons(self):
        self.dpl_img = Image("photo", file="img/deploy.png")
        self.dpl_bt = ttk.Button(self.root, text="Simulate deployment",
                style = "BW.TLabel", image=self.dpl_img, compound="bottom",
                command=self.deployPopup)
        self.dpl_bt.grid(row=0, column=0, columnspan=5, sticky="NESW")
        self.tdoa_img = Image("photo", file="img/TDOA.png")
        self.loc_bt = ttk.Button(self.root, text="Simulate localization",
                style="BW.TLabel", image=self.tdoa_img, compound="bottom",
                command=self.notAvailable)
        self.loc_bt.grid(row=1, column=0, columnspan=5, sticky="NESW")
        self.dpl = Dpl()

        s = ttk.Style()
        s.configure('nextPrev.TLabel', font=('TkDefaultFont', 10, 'bold'),
                    foreground=self.fg, background=self.bg, anchor='center')
        s.map("nextPrev.TLabel", background=[("active", self.active_bg)])
        self.prev_day = ttk.Button(self.root, text="|◀◀",
                style = "nextPrev.TLabel", command=self.previousDay)
        self.prev_day.grid(row=2, column=0, sticky="NESW")
        self.prev_min = ttk.Button(self.root, text="|◀",
                style = "nextPrev.TLabel", command=self.previousMinute)
        s.configure('play.TLabel', font=('TkDefaultFont', 16, 'bold'),
                    foreground=self.fg, background=self.bg, anchor='center')
        s.map("play.TLabel", background=[("active", self.active_bg)])
        self.prev_min.grid(row=2, column=1, sticky="NESW")
        self.play = ttk.Button(self.root, text="▶️",
                style = "play.TLabel", command=self.currentMinute)
        self.play.grid(row=2, column=2, sticky="NESW")
        self.next_min = ttk.Button(self.root, text="▶️|",
                style = "nextPrev.TLabel", command=self.nextMinute)
        self.next_min.grid(row=2, column=3, sticky="NESW")
        self.next_day = ttk.Button(self.root, text="▶️▶️|",
                style = "nextPrev.TLabel", command=self.nextDay)
        self.next_day.grid(row=2, column=4, sticky="NESW")
        self.dt_lbl = Label(self.root, text="Δt: 0 minutes", width=19,
                font="TkDefaultFont 10 bold", bg=self.bg, fg=self.fg)
        self.dt_lbl.grid(row=3, column=0, columnspan=5)

    def deployPopup(self):
        self.popup = Tk()
        self.popup.title("Deployment settings")
        self.showCurrentSats(0, rowspan=7)
        self.deployer_lbl = Label(self.popup,
                text="Satellite deployer:")
        self.deployer_lbl.grid(row=0, column=1, sticky="W")
        self.deployer_name_lbl = Label(self.popup,
                text="No sat selected")
        self.deployer_name_lbl.grid(row=0, column=2, sticky="W")
        self.dpl_mass_lbl1 = Label(self.popup,
                text="Deployer's mass [kg]:")
        self.dpl_mass_lbl1.grid(row=1, column=1, sticky="W")
        self.mass_box1 = Entry(self.popup)
        self.mass_box1.grid(row=1, column=2)
        self.dpl_name_lbl = Label(self.popup,
                text="New sat's name:")
        self.dpl_name_lbl.grid(row=2, column=1, sticky="W")
        self.name_box = Entry(self.popup)
        self.name_box.grid(row=2, column=2)
        self.dpl_cat_lbl = Label(self.popup,
                text="New sat's category:")
        self.dpl_cat_lbl.grid(row=3, column=1, sticky="W")
        self.cat_box = Entry(self.popup)
        self.cat_box.grid(row=3, column=2)
        self.dpl_mass_lbl2 = Label(self.popup,
                text="New sat's mass [kg]:")
        self.dpl_mass_lbl2.grid(row=4, column=1, sticky="W")
        self.mass_box2 = Entry(self.popup)
        self.mass_box2.grid(row=4, column=2)
        self.dpl_spdx_lbl = Label(self.popup,
                text="Deployment speed x [m/s]:")
        self.dpl_spdx_lbl.grid(row=5, column=1, sticky="W")
        self.spdx_box = Entry(self.popup)
        self.spdx_box.grid(row=5, column=2)
        self.dpl_spdy_lbl = Label(self.popup,
                text="Deployment speed y [m/s]:")
        self.dpl_spdy_lbl.grid(row=6, column=1, sticky="W")
        self.spdy_box = Entry(self.popup)
        self.spdy_box.grid(row=6, column=2)
        self.dpl_spdz_lbl = Label(self.popup,
                text="Deployment speed z [m/s]:")
        self.dpl_spdz_lbl.grid(row=7, column=1, sticky="W")
        self.spdz_box = Entry(self.popup)
        self.spdz_box.grid(row=7, column=2)
        self.deploy_now_bt = Button(self.popup, text="Deploy",
                                 command=self.deploySat)
        self.deploy_now_bt.grid(row=8, column=0,
                columnspan=3, sticky="NESW")
        self.popup.bind("<Button-1>", self.selectDeployer)
        self.popup.bind("<Return>", self.selectDeployer)
        self.popup.protocol("WM_DELETE_WINDOW", self.popup.destroy)
        self.popup.mainloop()

    def selectDeployer(self, event):
        #select = self.curr_sats_lst.curselection()
        if (event.widget == self.curr_sats_lst):
            select = self.curr_sats_lst.curselection()
            deployer_name = self.curr_sats_lst.get(select)
            self.deployer_name_lbl['text'] = deployer_name
        elif (event.widget == self.deploy_now_bt):
            self.deploySat()

    def deploySat(self):
        deployer_name = self.deployer_name_lbl['text']
        dplyr_mass = float(self.mass_box1.get())
        dplyd_mass = float(self.mass_box2.get())
        spdx = int(self.spdx_box.get())
        spdy = int(self.spdy_box.get())
        spdz = int(self.spdz_box.get())
        name = self.name_box.get()
        cat = self.cat_box.get()
        for sat in self.Sats:
            if (sat.name == deployer_name):
                deployer = sat
                break
        date = datetime.utcnow() + timedelta(minutes=self.dmin)
        newSat = self.dpl.deploy(cat, deployer, dplyr_mass,
                dplyd_mass, name, [spdx, spdy, spdz], date=date)
        self.ax_cov.append(self.ax.fill([0,0], [0,0], transform=Geodetic(),
                           color='white', alpha=self.cov_alpha)[0])
        self.sat_txt.append(self.ax.text([], [], "", color='yellow', size=8,
                            transform=Geodetic(), ha="center"))
        self.Sats.append(newSat)
        self.sortSats()
        self.popup.destroy()

    def previousDay(self):
        self.dmin -= 1440
        self.format_dt()

    def previousMinute(self):
        self.dmin -= 1
        self.format_dt()

    def currentMinute(self):
        self.dmin = 0
        self.dt_lbl['text'] = "Δt: 0 minutes".format(self.dmin)

    def nextMinute(self):
        self.dmin += 1
        self.format_dt()

    def nextDay(self):
        self.dmin += 1440
        self.format_dt()

    def format_dt(self):
        if (abs(self.dmin) > 1439):
            self.dt_lbl['text'] = "Δt: {:+.2f} days".format(self.dmin/1440)
        elif (abs(self.dmin) > 59):
            self.dt_lbl['text'] = "Δt: {:+.2f} hours".format(self.dmin/60)
        else:
            self.dt_lbl['text'] = "Δt: {:+d} minutes".format(self.dmin)

    def setCanvas(self):
        canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        #self.root.rowconfigure(0, weight=1)
        #self.root.columnconfigure(0, weight=1)
        canvas.get_tk_widget().grid(row=0, column=5, rowspan=4,
                columnspan=13, sticky="NES")
        canvas.draw()

    def setTableTitles(self):
        self.root.rowconfigure(4, weight=1)
        Label(self.root, text="Satellite", font="TkDefaultFont 10 bold", 
                bg=self.bg, fg=self.fg, width=19, anchor='w').grid(row=4,
                        column=0, columnspan=5, sticky="W")
        Label(self.root, text="Category", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg, width=20, anchor='w').grid(row=4,
                        column=5, sticky="W")
        Label(self.root, text="Latitude", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg).grid(row=4, column=6)
        Label(self.root, text="Longitude", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg).grid(row=4, column=7)
        Label(self.root, text="Alt. [km]", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg).grid(row=4, column=8)
        Label(self.root, text="Spd. [m/s]", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg).grid(row=4, column=9)
        Label(self.root, text="a [km]", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg).grid(row=4, column=10)
        Label(self.root, text="h [km²/s]", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg).grid(row=4, column=11)
        Label(self.root, text="e", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg).grid(row=4, column=12)
        Label(self.root, text="Ω", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg).grid(row=4, column=13)
        Label(self.root, text="i", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg).grid(row=4, column=14)
        Label(self.root, text="ω", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg).grid(row=4, column=15)
        Label(self.root, text="T. Anom.", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg).grid(row=4, column=16)

    def setTableContent(self):
        self.name_bt = []
        self.cat_lbl = []
        self.lat_lbl = []
        self.lng_lbl = []
        self.alt_lbl = []
        self.spd_lbl = []
        self.a_lbl = []
        self.h_lbl = []
        self.e_lbl = []
        self.raan_lbl = []
        self.i_lbl = []
        self.w_lbl = []
        self.theta_lbl = []
        self.top_index = 0
        self.bottom_index = 4*(len(self.Sats) > 4) + len(self.Sats)*(len(self.Sats) <= 4)
        for i in range(0, 17):
            self.root.columnconfigure(i, weight=1)
        for i in range(0, 6):
            self.name_bt.append(ttk.Button(self.root, style = "BW.TLabel"))
            self.cat_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=24, anchor='w'))
            self.lat_lbl.append(Label(self.root, bg=self.bg, fg=self.fg,  width=9, anchor='e'))
            self.lng_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=10, anchor='e'))
            self.alt_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=9, anchor='e'))
            self.spd_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=10, anchor='e'))
            self.a_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=7, anchor='e'))
            self.h_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=8, anchor='e'))
            self.e_lbl.append(Label(self.root, bg=self.bg, fg=self.fg))
            self.raan_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=7, anchor='e'))
            self.i_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=7, anchor='e'))
            self.w_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=7, anchor='e'))
            self.theta_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=7, anchor='e'))
        self.updateTableContent()
        for i in range(self.top_index, self.bottom_index):
            #self.root.rowconfigure(i+5, weight=1)
            self.rememberRow(i+1)
        self.up_bt = ttk.Button(self.root, text="▲", style = "BW.TLabel", command=self.up)
        self.up_bt.grid(row=5, column=17, rowspan=2, sticky="NESW")
        self.down_bt = ttk.Button(self.root, text="▼", style = "BW.TLabel", command=self.down)
        self.down_bt.grid(row=7, column=17, rowspan=2, sticky="NESW")
        self.forgetLastRows()

    def updateTableContent(self):
        rad2deg = 180/pi
        self.updtCnt += 1
        if (self.updtCnt > 120):
            self.refreshBackgroundImg()
            self.updtCnt = 0
        date = datetime.utcnow() + timedelta(minutes=self.dmin)
        for Sat in self.Sats:
            #Sat.updateWithDragEffect()
            #Sat.updateOrbitalParameters2()
            Sat.updateOrbitalParameters3(date)
        i = 0
        for Sat in self.Sats[self.top_index:self.bottom_index]:
            self.name_bt[i]['text'] = Sat.name
            self.name_bt[i]['command'] = lambda Sat=Sat: self.changeMainSat(Sat)
            self.cat_lbl[i]['text'] = Sat.getCategory()
            self.lat_lbl[i]['text'] = "{:0.4f}{}".format(Sat.getLat(), "°")
            self.lng_lbl[i]['text'] = "{:0.4f}{}".format(Sat.getLng(date=date), "°")
            self.alt_lbl[i]['text'] = "{:0.1f}".format((Sat.getAlt()*0.001))
            self.spd_lbl[i]['text'] = "{}".format(int(Sat.getSpeed()))
            self.a_lbl[i]['text'] = "{:0.1f}".format((Sat.getSemiMajorAxis()*0.001))
            self.h_lbl[i]['text'] = "{:0.1f}".format(Sat.getSpecAngMomentum()*0.000001)
            self.e_lbl[i]['text'] = "{:0.4f}".format(Sat.getEccentricity())
            self.raan_lbl[i]['text'] = "{:0.2f}{}".format((Sat.getRAAN()*rad2deg), "°")
            self.i_lbl[i]['text'] = "{:0.2f}{}".format((Sat.getInclination()*rad2deg), "°")
            self.w_lbl[i]['text'] = "{:0.2f}{}".format((Sat.getArgPerigee()*rad2deg), "°")
            self.theta_lbl[i]['text'] = "{:0.2f}{}".format((Sat.getAnomaly()*rad2deg), "°")
            i += 1

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

    def saveAs(self):
        file_name = asksaveasfilename()
        if file_name is None:
            return
        self.fig.savefig(file_name)

    def notAvailable(self):
        print("This command is not available yet")

    def geometry(self):
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        screen_resolution = "{}{}{}".format(str(screen_width),
                                            'x',
                                            str(screen_height))
        self.root.geometry(screen_resolution)

    def removeNight(self):
        self.night_alpha = 0
        self.editmenu.entryconfigure(3, label="Add day/night terminator",
                                     command=self.addNight)

    def addNight(self):
        self.night_alpha = 0.7
        self.editmenu.entryconfigure(3, label="Remove day/night terminator",
                                     command=self.removeNight)

    def removeCoverage(self):
        self.cov_alpha = 0
        self.editmenu.entryconfigure(4, label="Add coverage",
                                     command=self.addCoverage)

    def addCoverage(self):
        self.cov_alpha = 0.2
        self.editmenu.entryconfigure(4, label="Remove coverage", 
                                     command=self.removeCoverage)

    def setSearchBox(self):
        self.srch_lbl = Label(self.popup, text="Search:")
        self.srch_lbl.grid(row=0, column=0)
        self.ent = StringVar()
        self.srch_box = Entry(self.popup, textvariable=self.ent)
        self.srch_box.grid(row=0, column=1)
        self.srch_box.focus()

    def searchSat(self, event):
        srch = self.srch_box.get().upper()
        self.match = [s for s in self.avail_sats if srch in s]
        self.match.sort(reverse=True)
        for i in range(0, self.avail_sats_lst.size()):
            self.avail_sats_lst.delete(0)
        for sat in self.match:
            self.avail_sats_lst.insert(0, sat)

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
        self.avail_sats_lst = Listbox(self.popup, width=28)
        self.readAllSats()
        for sat in self.avail_sats:
            self.avail_sats_lst.insert(END, sat)
        self.avail_sats_lst.grid(row=1, column=0, rowspan=2, columnspan=2)
        self.popup.columnconfigure(0, weight=1)

    def showCurrentSats(self, col, rowspan=2):
        self.curr_lbl = Label(self.popup, text="Current satellites:", anchor='w')
        self.curr_lbl.grid(row=0, column=col, sticky="W")
        self.curr_sats_lst = Listbox(self.popup, width=28)
        for sat in self.Sats:
            self.curr_sats_lst.insert(END, sat.name)
        self.curr_sats_lst.grid(row=1, column=col, rowspan=rowspan)

    def addRemoveButtons(self):
        self.add_sat_bt = Button(self.popup, text="→",
                                 command=self.addSat)
        self.add_sat_bt.grid(row=1, column=2)
        self.remove_sat_bt = Button(self.popup, text="←",
                                    command=self.removeSat)
        self.remove_sat_bt.grid(row=2, column=2)

    def sortSats(self):
        self.Sats.sort(key = lambda s: s.name)

    def addSat(self):
        add_sat = self.avail_sats_lst.get(self.avail_sats_lst.curselection())
        self.ax_cov.append(self.ax.fill([0,0], [0,0], transform=Geodetic(),
                           color='white', alpha=self.cov_alpha)[0])
        self.sat_txt.append(self.ax.text([], [], "", color='yellow', size=8,
                            transform=Geodetic(), ha="center"))
        if (add_sat in self.argos):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/argos.txt"), cat="Argos Data Collection System"))
        elif (add_sat in self.cubesat):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/cubesat.txt"), cat="CubeSat"))
        elif (add_sat in self.dmc):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/dmc.txt"), cat="Disaster Monitoring"))
        elif (add_sat in self.goes):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/goes.txt"), cat="GOES"))
        elif (add_sat in self.intelsat):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/intelsat.txt"), cat="Intelsat"))
        elif (add_sat in self.iridium):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/iridium.txt"), cat="Iridium"))
        elif (add_sat in self.iridium_next):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/iridium-NEXT.txt"), cat="Iridium Next"))
        elif (add_sat in self.molniya):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/molniya.txt"), cat="Molniya"))
        elif (add_sat in self.noaa):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/noaa.txt"), cat="NOAA"))
        elif (add_sat in self.planet):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/planet.txt"), cat="Planet"))
        elif (add_sat in self.resource):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/resource.txt"), cat="Earth Resources"))
        elif (add_sat in self.sarsat):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/sarsat.txt"), cat="Search & Rescue"))
        elif (add_sat in self.spire):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/spire.txt"), cat="Spire"))
        elif (add_sat in self.tdrss):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/tdrss.txt"), cat="Tracking and Data Relay"))
        elif (add_sat in self.tle_new):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/tle-new.txt"), cat="Last 30 Days' Launches"))
        elif (add_sat in self.weather):
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat, "TLE/weather.txt"), cat="Weather"))
        else:
            self.Sats.append(Sat(add_sat, tle=tlefile.read(add_sat)))
        self.curr_sats_lst.insert(END, add_sat)
        self.sortSats()
        self.srch_box.focus()
        for i, sat in enumerate(self.Sats):
            self.curr_sats_lst.delete(i)
            self.curr_sats_lst.insert(i, sat.name)
        

    def removeSat(self):
        del_sat = self.curr_sats_lst.get(self.curr_sats_lst.curselection())
        self.ax_cov.pop()
        self.sat_txt.pop()
        for i, sat in enumerate(self.Sats):
            if (del_sat == sat.name):
                self.Sats.remove(sat)
                self.curr_sats_lst.delete(i)
        
    def addRemoveSat(self):
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

    def fullscreen(self, event=None):
        self.root.attributes("-fullscreen", True)
        self.viewmenu.entryconfigure(0, label="Exit fullscreen", 
                                     command=self.exitFullscreen)
        self.rememberLastRows()

    def exitFullscreen(self, event=None):
        self.root.attributes("-fullscreen", False)
        self.viewmenu.entryconfigure(0, label="Fullscreen",
                                     command=self.fullscreen)
        self.forgetLastRows()

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
        self.ax.clear()
        self.plotData()
        self.setCanvas()
        if (img is None):
            date = datetime.utcnow() + timedelta(minutes=self.dmin)
            self.map.set_data(self.world_map.fillDarkSideFromPicture(date))
        else:
            self.map.set_data(imread(img))
        self.ani._stop()
        self.ani._blit_clear(self.ani._drawn_artists, self.ani._blit_cache)
        self.fig.canvas.draw_idle()
        self.ani.__init__(self.fig, self.update, self.data_gen(),
                            interval=self.dt, blit=True, repeat=False)

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
        self.popup = Tk()
        self.popup.title("About Pypredict")
        self.prog_name_lbl = Label(self.popup,
                text="Pypredict", font="TkDefaultFont 10 bold")
        self.prog_name_lbl.grid(row=0, column=0, sticky="EW")
        self.version_lbl = Label(self.popup,
                text=__version__)
        self.version_lbl.grid(row=1, column=0, sticky="EW")
        self.dev_lbl = Label(self.popup,
                text="\nCopyright (C) 2018-2019, Matías Vidal Valladares.")
        self.dev_lbl.grid(row=2, column=0, sticky="EW")
        self.contact_lbl = Label(self.popup,
                text="E-mail: matias.vidal.v@gmail.com")
        self.contact_lbl.grid(row=3, column=0, sticky="EW")
        self.warranty_lbl = Label(self.popup,
                text="\nThis program comes with ABSOLUTELY NO WARRANTY.")
        self.warranty_lbl.grid(row=4, column=0, sticky="EW")
        self.details_lbl = Label(self.popup,
                text="See the GNU General Public License, version 3 or later for details.")
        self.details_lbl.grid(row=5, column=0, sticky="EW")
        self.popup.columnconfigure(0,weight=1)
        self.popup.columnconfigure(1,weight=1)
        self.popup.columnconfigure(2,weight=1)
        self.popup.columnconfigure(3,weight=1)
        self.popup.protocol("WM_DELETE_WINDOW", self.popup.destroy)
        self.popup.mainloop()

    def setMenu(self):
        menubar = Menu(self.root, bg=self.bg, fg=self.fg, activeforeground=self.fg,
                        activebackground=self.active_bg)

        # Create a pulldown menu, and add it to the menu bar
        filemenu = Menu(menubar, tearoff=0, bg=self.active_bg, fg=self.fg,
                        activeforeground=self.fg, activebackground=self.bg)
        filemenu.add_command(label="Open", command=self.notAvailable)
        filemenu.add_command(label="Save as", command=self.saveAs)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=exit)
        menubar.add_cascade(label="File", menu=filemenu)

        # Create more pulldown menus
        self.editmenu = Menu(menubar, tearoff=0, bg=self.active_bg, fg=self.fg,
                             activeforeground=self.fg, activebackground=self.bg)
        self.editmenu.add_command(label="Update TLE from net", command=self.updateTLEfromNet)
        self.editmenu.add_command(label="Update TLE from file", command=self.notAvailable)
        self.editmenu.add_separator()
        self.editmenu.add_command(label="Add day/night terminator", command=self.addNight)
        self.editmenu.add_command(label="Remove coverage", command=self.removeCoverage)
        self.editmenu.add_command(label="Add/remove satellites", command=self.addRemoveSat)
        menubar.add_cascade(label="Edit", menu=self.editmenu)
        
        self.viewmenu = Menu(menubar, tearoff=0, bg=self.active_bg, fg=self.fg,
                             activeforeground=self.fg, activebackground=self.bg)
        self.viewmenu.add_command(label="Fullscreen", command=self.fullscreen)
        #self.viewmenu.add_command(label="Change planet", command=self.notAvailable)
        self.planetmenu = Menu(self.viewmenu, tearoff=0, bg=self.active_bg, fg=self.fg,
                               activeforeground=self.fg, activebackground=self.bg)
        self.planetmenu.add_command(label="Earth", command=self.earth)
        self.planetmenu.add_command(label="Mars", command=self.mars)
        self.viewmenu.add_cascade(label="Change planet", menu=self.planetmenu)
        menubar.add_cascade(label="View", menu=self.viewmenu)

        helpmenu = Menu(menubar, tearoff=0, bg=self.active_bg, fg=self.fg,
                        activeforeground=self.fg, activebackground=self.bg)
        helpmenu.add_command(label="About Pypredict", command=self.about)
        menubar.add_cascade(label="Help", menu=helpmenu)

        # Display the menu
        self.root.config(menu=menubar)

    def run(self):
        self.root.mainloop()

