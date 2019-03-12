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
from pyorbital import tlefile
from sat import Sat
from gui import GUI

SUCHAI = Sat(name="SUCHAI", tle=tlefile.read("SUCHAI", "TLE/cubesat.txt"), cat="CubeSat")
HODOYOSHI3 = Sat(name="HODOYOSHI-3", tle=tlefile.read("HODOYOSHI-3", "TLE/resource.txt"), cat="Earth Resources")
HODOYOSHI4 = Sat(name="HODOYOSHI-4", tle=tlefile.read("HODOYOSHI-4", "TLE/resource.txt"), cat="Earth Resources")
CUBESATXI_IV = Sat(name="CUBESAT XI-IV (CO-57)", tle=tlefile.read("CUBESAT XI-IV (CO-57)", "TLE/cubesat.txt"), cat="CubeSat")

Sats = [SUCHAI, HODOYOSHI3, HODOYOSHI4, CUBESATXI_IV]
GUI = GUI(Sats=Sats)
