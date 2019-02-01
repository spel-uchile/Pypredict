from cartopy.crs import Geodetic, PlateCarree, RotatedPole
from matplotlib.ticker import FixedLocator
from numpy import abs, arange, array, cos, empty, log, pi, sin, tan
from matplotlib.animation import FuncAnimation
from pyorbital import tlefile
from matplotlib.pyplot import imread, subplots, tight_layout
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import datetime
from tkinter import Button, Entry, END, Image, Label, Listbox, Menu, StringVar, ttk, Tk
from tkinter.filedialog import asksaveasfilename
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from sat import Sat
from warnings import filterwarnings
import ssl

filterwarnings("ignore", category=RuntimeWarning)
ssl._create_default_https_context = ssl._create_unverified_context

class GUI():
    #@profile
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
        tf = int(self.mainSat.getPeriod()*3)           # Total duration in seconds
        self.dt = 1                                    # Step's length in seconds
        self.mainSat_lats, self.mainSat_lngs = self.mainSat.getLocation(tf, self.dt)
        self.plotData()
        self.setCanvas()
        self.setMenu()
        self.setTableTitles()
        self.setTableContent()
        self.tableRefresher()
        ani = FuncAnimation(self.fig, self.update, self.data_gen(),
                            interval=self.dt*964, blit=True, repeat=False)
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
        else:
            self.fig, self.ax = subplots(figsize=(18, 9),
                                    subplot_kw={'projection': PlateCarree()})
        img_extent = (-180, 180, -90, 90)
        self.map = self.ax.imshow(imread("img/earth_nasa_day.png"), origin='upper',
                       extent=img_extent, transform=PlateCarree())
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
            for Sat in self.Sats:
                Sat.updateOrbitalParameters()
                sats_lngs.append(Sat.getLng())
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
        lng = empty(360)
        lat = empty(360)
        for i in range(0, 360):
            theta = i*deg2rad
            dlat = ang*deg2rad * cos(theta)
            lat[359 - i] = sat_lat*deg2rad + dlat
            dpsi = log(tan(lat[359 - i]/2 + pi/4)/tan(sat_lat*deg2rad/2 + pi/4))
            if (abs(dpsi) > 10e-12):
                q = dlat / dpsi
            else:
                q = cos(sat_lat*deg2rad)
            dlng = ang*deg2rad*sin(theta)/q
            lng[359 - i] = sat_lng*deg2rad + dlng
            if (abs(lat[359 - i]) > (6*deg2rad + pi - ang*deg2rad - abs(sat_lat*deg2rad))):
                lat[359 - i] = (2*(sat_lat > 0) - 1)*(6*deg2rad + pi - ang*deg2rad - abs(sat_lat*deg2rad))
                lng[359 - i] = sat_lng*deg2rad + ((sat_lng*deg2rad - pi) > 0)*(-pi) + ((sat_lng*deg2rad - pi) <= 0)*pi
            lat[359 - i] = lat[359 - i]*rad2deg
            lng[359 - i] = lng[359 - i]*rad2deg
        self.ax_cov[n].set_alpha(self.cov_alpha)    
        self.ax_cov[n].set_xy(array((lng, lat)).transpose())

    def setTheme(self, bg, fg, active_bg):
        self.bg = bg
        self.fg = fg
        self.active_bg = active_bg
        self.root.configure(background=bg)
        style = ttk.Style()
        style.configure("BW.TLabel", foreground=self.fg, background=self.bg)
        style.map("BW.TLabel", background=[("active", self.active_bg)])
        self.night_alpha = 0.7
        self.cov_alpha = 0.2

    def changeMainSat(self, Sat):
        self.mainSat = Sat
        tf = int(self.mainSat.getPeriod()*3)
        self.mainSat_lats, self.mainSat_lngs = self.mainSat.getLocation(tf, self.dt)
   
    def up(self):
        if (self.top_index > 0):
            self.top_index -= 1
            self.bottom_index -= 1

    def down(self):
        if (self.bottom_index < len(self.Sats)):
            self.top_index += 1
            self.bottom_index += 1

    def setCanvas(self):
        canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        #self.root.rowconfigure(0, weight=1)
        #self.root.columnconfigure(0, weight=1)
        canvas.get_tk_widget().grid(row=0, column=0, columnspan=13, sticky="NES")
        canvas.draw()

    def setTableTitles(self):
        #self.root.rowconfigure(1, weight=1)
        Label(self.root, text="Satellite", font="TkDefaultFont 10 bold", 
                bg=self.bg, fg=self.fg, width=18, anchor='w').grid(row=1, column=0, sticky="W")
        Label(self.root, text="Category", font="TkDefaultFont 10 bold",
                bg=self.bg, fg=self.fg, width=24, anchor='w').grid(row=1, column=1, sticky="W")
        Label(self.root, text="Latitude", font="TkDefaultFont 10 bold", bg=self.bg, fg=self.fg).grid(row=1, column=2)
        Label(self.root, text="Longitude", font="TkDefaultFont 10 bold", bg=self.bg, fg=self.fg).grid(row=1, column=3)
        Label(self.root, text="Alt. [km]", font="TkDefaultFont 10 bold", bg=self.bg, fg=self.fg).grid(row=1, column=4)
        Label(self.root, text="a [km]", font="TkDefaultFont 10 bold", bg=self.bg, fg=self.fg).grid(row=1, column=5)
        Label(self.root, text="e", font="TkDefaultFont 10 bold", bg=self.bg, fg=self.fg).grid(row=1, column=6)
        Label(self.root, text="Ω", font="TkDefaultFont 10 bold", bg=self.bg, fg=self.fg).grid(row=1, column=7)
        Label(self.root, text="i", font="TkDefaultFont 10 bold", bg=self.bg, fg=self.fg).grid(row=1, column=8)
        Label(self.root, text="ω", font="TkDefaultFont 10 bold", bg=self.bg, fg=self.fg).grid(row=1, column=9)
        Label(self.root, text="T. Anom.", font="TkDefaultFont 10 bold", bg=self.bg, fg=self.fg).grid(row=1, column=10)
        Label(self.root, text="M. Anom.", font="TkDefaultFont 10 bold", bg=self.bg, fg=self.fg).grid(row=1, column=11)

    def setTableContent(self):
        rad2deg = 180/pi
        self.name_bt = []
        self.cat_lbl = []
        self.lat_lbl = []
        self.lng_lbl = []
        self.alt_lbl = []
        self.a_lbl = []
        self.e_lbl = []
        self.raan_lbl = []
        self.i_lbl = []
        self.w_lbl = []
        self.theta_lbl = []
        self.ma_lbl = []
        self.top_index = 0
        self.bottom_index = 4*(len(self.Sats) > 4) + len(self.Sats)*(len(self.Sats) <= 4)
        for i in range(0, 12):
            self.root.columnconfigure(i, weight=1)
        for i in range(0, 6):
            self.name_bt.append(ttk.Button(self.root, style = "BW.TLabel"))
            self.cat_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=24, anchor='w'))
            self.lat_lbl.append(Label(self.root, bg=self.bg, fg=self.fg,  width=9, anchor='e'))
            self.lng_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=10, anchor='e'))
            self.alt_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=9, anchor='e'))
            self.a_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=7, anchor='e'))
            self.e_lbl.append(Label(self.root, bg=self.bg, fg=self.fg))
            self.raan_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=7, anchor='e'))
            self.i_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=7, anchor='e'))
            self.w_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=7, anchor='e'))
            self.theta_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=7, anchor='e'))
            self.ma_lbl.append(Label(self.root, bg=self.bg, fg=self.fg, width=7, anchor='e'))
        self.updateTableContent()
        for i in range(self.top_index, self.bottom_index):
            self.root.rowconfigure(i+2, weight=1)
            self.rememberRow(i+1)
        self.up_bt = ttk.Button(self.root, text="▲", style = "BW.TLabel", command=self.up)
        self.up_bt.grid(row=2, column=11, rowspan=2, sticky="NESW")
        self.down_bt = ttk.Button(self.root, text="▼", style = "BW.TLabel", command=self.down)
        self.down_bt.grid(row=4, column=11, rowspan=2, sticky="NESW")
        self.forgetLastRows()

    def updateTableContent(self):
        rad2deg = 180/pi
        i = 0
        for Sat in self.Sats[self.top_index:self.bottom_index]:
            Sat.updateOrbitalParameters()
            self.name_bt[i]['text'] = Sat.name
            self.name_bt[i]['command'] = lambda Sat=Sat: self.changeMainSat(Sat)
            self.cat_lbl[i]['text'] = Sat.getCategory()
            self.lat_lbl[i]['text'] = "{:0.4f}{}".format(Sat.getLat(), "°")
            self.lng_lbl[i]['text'] = "{:0.4f}{}".format(Sat.getLng(), "°")
            self.alt_lbl[i]['text'] = "{:0.1f}".format((Sat.getAlt()/1000))
            self.a_lbl[i]['text'] = "{:0.2f}".format((Sat.getSemiMajorAxis()/1000))
            self.e_lbl[i]['text'] = "{:0.4f}".format(Sat.getEccentricity())
            self.raan_lbl[i]['text'] = "{:0.2f}{}".format((Sat.getRAAN()*rad2deg), "°")
            self.i_lbl[i]['text'] = "{:0.2f}{}".format((Sat.getInclination()*rad2deg), "°")
            self.w_lbl[i]['text'] = "{:0.2f}{}".format((Sat.getArgPerigee()*rad2deg), "°")
            self.theta_lbl[i]['text'] = "{:0.2f}{}".format((Sat.getAnomaly()*rad2deg), "°")
            self.ma_lbl[i]['text'] = "{:0.2f}{}".format((Sat.getMeanAnomaly()*rad2deg), "°")
            i += 1

    def rememberLastRows(self):
        if (len(self.Sats) > 4):
            self.top_index = 0
            self.bottom_index = 5
            self.rememberRow(self.bottom_index)
            if (len(self.Sats) >= 6):
                self.bottom_index = 6
                self.rememberRow(self.bottom_index)
                self.up_bt.grid(row=2, column=12, rowspan=3, sticky="NESW")
                self.down_bt.grid(row=5, column=12, rowspan=3, sticky="NESW")

    def rememberRow(self, r):
        self.root.rowconfigure(r+1, weight=1)
        self.name_bt[r-1].grid(row=r+1, column=0, sticky="W")
        self.cat_lbl[r-1].grid(row=r+1, column=1, sticky="W")
        self.lat_lbl[r-1].grid(row=r+1, column=2)
        self.lng_lbl[r-1].grid(row=r+1, column=3)
        self.alt_lbl[r-1].grid(row=r+1, column=4)
        self.a_lbl[r-1].grid(row=r+1, column=5)
        self.e_lbl[r-1].grid(row=r+1, column=6)
        self.raan_lbl[r-1].grid(row=r+1, column=7)
        self.i_lbl[r-1].grid(row=r+1, column=8)
        self.w_lbl[r-1].grid(row=r+1, column=9)
        self.theta_lbl[r-1].grid(row=r+1, column=10)
        self.ma_lbl[r-1].grid(row=r+1, column=11)

    def forgetLastRows(self):
        self.top_index = 0
        self.bottom_index = 4*(len(self.Sats) > 4) + len(self.Sats)*(len(self.Sats) <= 4)
        self.up_bt.grid(row=2, column=12, rowspan=2, sticky="NESW")
        self.down_bt.grid(row=4, column=12, rowspan=2, sticky="NESW")
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
        self.a_lbl[r].grid_forget()
        self.e_lbl[r].grid_forget()
        self.raan_lbl[r].grid_forget()
        self.i_lbl[r].grid_forget()
        self.w_lbl[r].grid_forget()
        self.theta_lbl[r].grid_forget()
        self.ma_lbl[r].grid_forget()

    def setTableContent3(self, root, Sats):
        rad2deg = 180/pi
        self.scrollbar = Scrollbar(root)
        self.scrollbar.grid(row=2, column=10)
        self.table_lst = Listbox(root, yscrollcommand = self.scrollbar.set, width=100)
        for i in range(0, 10):
            root.columnconfigure(i, weight=1)
        for Sat in Sats:
            Sat.updateOrbitalParameters()
            data = '{} {:f}{} {:f}{} {:0.1f} {:0.2f} {:0.4f} {:0.2f}{} {:0.2f}{} {:0.2f}{} {:0.2f}{}'.format(Sat.name,
                                                                                           Sat.getLat(), '°',
                                                                                           Sat.getLng(), '°',
                                                                                           (Sat.getAlt()/1000),
                                                                                           (Sat.getSemiMajorAxis()/1000),
                                                                                           Sat.getEccentricity(),
                                                                                           (Sat.getRAAN()*rad2deg), '°',
                                                                                           (Sat.getInclination()*rad2deg), '°',
                                                                                           (Sat.getArgPerigee()*rad2deg), '°',
                                                                                           (Sat.getAnomaly()*rad2deg % 360), '°')
            
            self.table_lst.insert(END, data)
            i += 1
        self.table_lst.grid(row=2, column=1, columnspan=10)
    
    def updateTableContent3(self, root, Sats):
        rad2deg = 180/pi
        i = 0
        for Sat in Sats:
            Sat.updateOrbitalParameters()
            data = '{} {:f}{} {:f}{} {:0.1f} {:0.2f} {:0.4f} {:0.2f}{} {:0.2f}{} {:0.2f}{} {:0.2f}{}'.format(Sat.name,
                                                                                           Sat.getLat(), '°',
                                                                                           Sat.getLng(), '°',
                                                                                           (Sat.getAlt()/1000),
                                                                                           (Sat.getSemiMajorAxis()/1000),
                                                                                           Sat.getEccentricity(),
                                                                                           (Sat.getRAAN()*rad2deg), '°',
                                                                                           (Sat.getInclination()*rad2deg), '°',
                                                                                           (Sat.getArgPerigee()*rad2deg), '°',
                                                                                           (Sat.getAnomaly()*rad2deg % 360), '°')
            self.table_lst.insert(i, data)
            self.table_lst.delete(i+1)
            i += 1

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
        self.avail_sats += self.noaa + self.planet + self.resource
        self.avail_sats += self.sarsat + self.spire + self.tdrss
        self.avail_sats += self.tle_new + self.weather
        self.avail_sats.sort()

    def showAvailSats(self):
        self.avail_sats_lst = Listbox(self.popup, width=28)
        self.readAllSats()
        for sat in self.avail_sats:
            self.avail_sats_lst.insert(END, sat)
        self.avail_sats_lst.grid(row=1, column=0, rowspan=2, columnspan=2)
        self.popup.columnconfigure(0, weight=1)

    def showCurrentSats(self):
        self.curr_lbl = Label(self.popup, text="Current satellites:", anchor='w')
        self.curr_lbl.grid(row=0, column=3, sticky="W")
        self.curr_sats_lst = Listbox(self.popup, width=28)
        for Sat in self.Sats:
            self.curr_sats_lst.insert(END, Sat.name)
        self.curr_sats_lst.grid(row=1, column=3, rowspan=2)

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
        self.srch_box.focus()
        self.sortSats()
        

    def removeSat(self):
        del_sat = self.curr_sats_lst.get(self.curr_sats_lst.curselection())
        self.ax_cov.pop()
        self.sat_txt.pop()
        for i, Sat in enumerate(self.Sats):
            if (del_sat == Sat.name):
                self.Sats.remove(Sat)
                self.curr_sats_lst.delete(i)
        
    def addRemoveSat(self):
        self.popup = Tk()
        self.popup.title("Add/remove satellites")
        #self.popup.call('wm', 'iconphoto', self.popup._w, self.img)
        self.showAvailSats()
        self.setSearchBox()
        self.addRemoveButtons()
        self.showCurrentSats()
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

    def earth(self):
        img_extent = (-180, 180, -90, 90)
        self.map.set_data(imread("img/earth_nasa_day.png"))
        self.fig.canvas.draw_idle()
        mu = 5.9722*6.67408*10**13
        for sat in self.Sats:
            sat.changePlanet()

    def mars(self):
        img_extent = (-180, 180, -90, 90) 
        self.map.set_data(imread("img/mars_nasa_day.png"))
        self.fig.canvas.draw_idle()
        mu = 0.64171*6.67408*10**13
        for sat in self.Sats:
            sat.changePlanet(M=0.64171*10**24, P_r=3389500, Eq_r=3396200, 
                    Po_r=3376200, J2=0.00196045, P_w=7.08821812*10**(-5))

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
        self.editmenu.add_command(label="Remove day/night terminator", command=self.removeNight)
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
        helpmenu.add_command(label="About Pypredict", command=self.notAvailable)
        menubar.add_cascade(label="Help", menu=helpmenu)

        # Display the menu
        self.root.config(menu=menubar)

    def run(self):
        self.root.mainloop()

