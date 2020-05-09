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
from cartopy.crs import PlateCarree, RotatedPole
from datetime import datetime
from numpy import arange, argmax, argmin, cos, empty, pi
from PIL import Image, ImageDraw

class Map(object):
    __slots__ = ["day", "night", "lat", "lng", "xy",
                 "sun_lat", "sun_lng", "dark_lat",
                 "dark_lng", "x", "y", "day_img", "night_img"]

    def __init__(self, day, night):
        """
        Initialize Map class. This class creates a composite
        of the world map in daylight and in night, given the
        current position of the sun.

        Parameters
        ----------
        day   : string
                Path to the picture in png of the world at day
        night : string
                Path to the picture in png of the world at
                night
        """
        self.day = day
        self.night = night
        self.day_img = Image.open(self.day).convert('RGBA')
        self.night_img = Image.open(self.night).convert('RGBA')
        self.lat = empty(360)
        self.lng = empty(360)
        self.x = empty(360)
        self.y = empty(360)
        self.xy = []

    def __call__(self):
        return self

    def getSunCoordinates(self, date=None):
        """
        Calculates the latitude and longitude of the sun
        and of the dark side.

        Parameters
        ----------
        date : datetime.utcnow(), optional
               Date to calculate the sun's position
        """
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

    def getCoverageCoordinates(self, center_lat, center_lng):
        """
        Returns the coordinates of the day/night terminator.

        center_lat : float
                     Center latitude of the dark side in degrees
        center_lng : float
                     Center longitude of the dark side in degrees
        """
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

        self.x[:180] = -90
        self.y[:180] = arange(-90, 90.)
        self.x[180:] = 90
        self.y[180:] = arange(90, -90., -1)

        transformed = rotated_pole.transform_points(PlateCarree(),
                                                    self.x, self.y)
        self.lng = transformed[:,0]
        self.lat = transformed[:,1]
        if (center_lat > 0):
            self.lat += 90
        else:
            self.lat = 90 - self.lat
        self.lng = (self.lng + 180 + center_lng) % 360

    def coord2pixels(self, center_lat, center_lng):
        """
        Transforms the coordinates to pixels.

        Parameters
        ----------
        center_lat : float
                     Center latitude of the dark side in degrees
        center_lng : float
                     Center longitude of the dark side in degrees
        """
        self.getCoverageCoordinates(center_lat, center_lng)
        lng2pix = 2200/360
        lat2pix = 1100/180
        self.lng = self.lng*lng2pix
        self.lat = self.lat*lat2pix
        self.xy = list(zip(self.lng.astype(int),
                           self.lat.astype(int)))
        if (max(self.lng) is not 2200):
            self.xy[argmax(self.lng)] = ((2200,
                self.lat[argmax(self.lng)]))
        if (min(self.lng) is not 0):
            self.xy[argmin(self.lng)] = ((0,
                self.lat[argmin(self.lng)]))

    def fillDarkSideFromPicture(self, date=None):
        """
        Returns a composite of the images of the Earth at
        day and at night, considering the sun's position.

        Parameters
        ----------
        date : datetime.utcnow(), optional
               Date to calculate the sun's position
        """
        self.getSunCoordinates(date)
        self.coord2pixels(self.dark_lat, self.dark_lng)
        self.xy = sorted(self.xy, key=lambda x: x[0])
        self.xy.append((2200, 0))
        self.xy.append((0, 0))
        if (self.dark_lat > 0):
            self.xy[len(self.xy) - 1] = (0, 1100)
            self.xy[len(self.xy) - 2] = (2200, 1100)
        temp_night = self.night_img.copy()
        drw2 = ImageDraw.Draw(temp_night, 'RGBA')
        drw2.polygon(self.xy, fill=(255, 255, 255, 0))
        del drw2
        composite = Image.alpha_composite(self.day_img, temp_night)
        temp_night.close()
        return composite
