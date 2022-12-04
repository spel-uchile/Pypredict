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
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT
from pkg_resources import resource_filename

class Toolbar(NavigationToolbar2QT):

    def __init__(self, canvas, applicationWindow):
        img_path = resource_filename("pypredict","img/")
        self.toolitems = (
            ('Home', 'Reset original view', '{}home'.format(img_path), 'home'),
            ('Back', 'Back to previous view', '{}back'.format(img_path), 'back'),
            ('Forward', 'Forward to next view', '{}forward'.format(img_path), 'forward'),
            (None, None, None, None),
            # ('Subplots', 'Configure subplots', 'subplots', 'configure_subplots'),
            ('Pan', 'Pan axes with left mouse, zoom with right', '{}move'.format(img_path), 'pan'),
            ('Zoom', 'Zoom to rectangle', '{}zoom'.format(img_path), 'zoom'),
            (None, None, None, None),
            # ('Save', 'Save the figure', 'filesave', 'save_figure'),
        )
        super(Toolbar, self).__init__(canvas, applicationWindow)
