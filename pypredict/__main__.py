#!/usr/bin/python3
"""
                                Pypredict
    Orbit prediction software. Displays the satellites' position and
    orbital parameters in real time. Simulates satellite localization
    and deployment.
    
    Copyright (C) 2018-2022, Matías Vidal Valladares, matvidal.
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
from pkg_resources import resource_filename
from pypredict.sat import Sat
from PyQt5 import QtWidgets
from pypredict.app import ApplicationWindow
import sys

data_path = resource_filename("pypredict","data/")

SUCHAI2 = Sat(name="SUCHAI-2", tlepath="{}active.txt".format(data_path), cat="CubeSat")
SUCHAI3 = Sat(name="SUCHAI-3", tlepath="{}active.txt".format(data_path), cat="CubeSat")
PLANTSAT = Sat(name="PLANTSAT", tlepath="{}active.txt".format(data_path), cat="CubeSat")
ISS = Sat(name="ISS (ZARYA)", tlepath="{}tdrss.txt".format(data_path), cat="Tracking and Data Relay")
Sats = [SUCHAI2, SUCHAI3, PLANTSAT, ISS]

def main(args=None):
    if args is None:
        args = sys.argv[1:]
    app = QtWidgets.QApplication(sys.argv)
    application = ApplicationWindow(Sats=Sats)
    application.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
