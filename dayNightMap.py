from cartopy.crs import Geodetic, PlateCarree, RotatedPole
from PIL import Image, ImageDraw
from datetime import datetime
from numpy import arange, arcsin, arctan2, argmax, argmin, cos, empty, log, pi, sin, tan
from warnings import filterwarnings

filterwarnings("ignore", category=RuntimeWarning)

class Map(object):
    __slots__ = ["day", "night", "lat", "lng", "xy", "sun_lat", "sun_lng", "dark_lat", "dark_lng"]
    def __init__(self, day, night):
        self.day = day
        self.night = night
        self.lat = empty(360)
        self.lng = empty(360)
        self.xy = []

    def __call__(self):
        return self

    def getSunCoordinates(self, date=None):
        if date is None:
            date = datetime.utcnow()
        D = date.day + date.hour/24 + date.minute/(24*60) + date.second/(24*3600)
        M = date.month
        Y = date.year
        JD = 367*Y - int(7/4*(Y + int((M + 9)/12))) + int(275*M/9) + D + 1721013.5
        T0 = (JD - 2451545.0)/36525
        axial_tilt = 23.43929111 + (-46.8150*T0 - 0.00059*T0*T0 + 0.001813*T0*T0*T0)/3600.0
        ref_solstice = datetime(2016, 6, 21, 22, 22)
        days_per_year = 365.2425
        seconds_per_day = 24*3600.0
        days_since_ref = (date - ref_solstice).total_seconds()/seconds_per_day
        self.sun_lat = axial_tilt*cos(2*pi*days_since_ref/days_per_year)
        sec_since_midnight = (date - datetime(date.year, date.month, date.day)).seconds
        self.sun_lng = -(sec_since_midnight/seconds_per_day - 0.5)*360
        self.dark_lat = -self.sun_lat
        if (self.sun_lng > 0):
            self.dark_lng = (((self.sun_lng - 180) + 180) % 360) - 180
        else:
            self.dark_lng = (((self.sun_lng + 180) + 180) % 360) - 180

    def getCoverageCoordinates(self, ang, center_lat, center_lng):
        pole_lng = 0#center_lng
        if center_lat > 0:
            pole_lat = -90 + center_lat
            central_rot_lng = 180
        else:
            pole_lat = 90 + center_lat
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

        z = empty(360)
        transformed = rotated_pole.transform_points(PlateCarree(), x, y)
        for i in range(0, 360):
            self.lng[i], self.lat[i], z[i] = transformed[i]
            self.lat[i] += 90
            self.lng[i] = (self.lng[i] + 180 + center_lng) % 360

    def coord2pixels(self, ang, center_lat, center_lng):
        self.getCoverageCoordinates(ang, center_lat, center_lng)
        self.xy[:] = []
        for i in range(0, len(self.lat)):
            self.xy.append((int(self.lng[i]*2200/360),
                            int(self.lat[i]*1100/180)))
        if (max(self.lng) is not 360):
            self.xy[argmax(self.lng)] = ((2200,
                int(self.lat[argmax(self.lng)]*1100/180)))
        if (min(self.lng) is not 0):
            self.xy[argmin(self.lng)] = ((0,
                int(self.lat[argmin(self.lng)]*1100/180)))

    def fillDarkSideFromPicture(self, ang, date=None):
        self.getSunCoordinates(date)
        self.coord2pixels(ang, self.dark_lat, self.dark_lng)
        self.xy = sorted(self.xy, key=lambda x: x[0])
        self.xy.append((2200, 0))
        self.xy.append((0, 0))
        img1 = Image.open(self.day).convert('RGBA')
        #drw1 = ImageDraw.Draw(img1, 'RGBA')
        #drw1.polygon(self.xy, fill=(255, 255, 255, 0))
        #del drw1

        self.xy[len(self.xy) - 1] = (0, 1100)
        self.xy[len(self.xy) - 2] = (2200, 1100)
        img2 = Image.open(self.night).convert('RGBA')
        drw2 = ImageDraw.Draw(img2, 'RGBA')
        drw2.polygon(self.xy, fill=(255, 255, 255, 0))
        del drw2

        composite = Image.alpha_composite(img1, img2)
        img1.close()
        img2.close()
        return composite
