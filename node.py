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
class Node(object):
    __slots__ = ["name", "lat", "lng", "alt", "freq"]
    def __init__(self, name="", lat=0, lng=0, alt=0, freq=437225000):
        """
        Parameters
        ----------
        name : str, optional
            The node's name
        lat : float, optional
            The node's, latitude (default is 0)
        lng : float, optional
            The node's longitude (default is 0)
        alt : float, optional
            The node's altitude in meters (default is 0)
        freq : int, optional
            The node's frequency in Hz (default is 437225000)
        """
        self.name = name
        self.lat = lat
        self.lng = lng
        self.alt = alt
        self.freq = freq

    def __call__(self):
        return self

    def __str__(self):
        """Returns the node's name."""
        return self.name

    def setLat(self, lat):
        """Sets the node's latitude."""
        self.lat = lat

    def setLng(self, lng):
        """Sets the node's longitude."""
        self.lng = lng

    def setAlt(self, alt):
        """Sets the node's altitude."""
        self.alt = alt

    def setFreq(self, freq):
        """Sets the node's frequency."""
        self.freq = freq

    def getLat(self):
        """Returns the node's latitude."""
        return self.lat

    def getLng(self):
        """Returns the node's longitude."""
        return self.lng

    def getAlt(self):
        """Returns the node's altitude."""
        return self.alt

    def getFreq(self):
        """Returns the node's frequency."""
        return self.freq
