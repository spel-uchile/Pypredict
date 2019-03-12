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
class Node(object):
    __slots__ = ["name", "lat", "lng", "alt", "freq"]
    def __init__(self, name="", lat=0, lng=0, alt=0, freq=437225000):
        self.name = name
        self.lat = lat
        self.lng = lng
        self.alt = alt
        self.freq = freq

    def __call__(self):
        return self

    def __str__(self):
        return self.name

    def setLat(self, lat):
        self.lat = lat

    def setLng(self, lng):
        self.lng = lng

    def setAlt(self, alt):
        self.alt = alt

    def setFreq(self, freq):
        self.freq = freq

    def getLat(self):
        return self.lat

    def getLng(self):
        return self.lng

    def getAlt(self):
        return self.alt

    def getFreq(self):
        return self.freq
