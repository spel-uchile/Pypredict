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
from matplotlib import path
from numpy import array

class SAA(object):
    __slots__ = ["lat", "lng", "path", "vertices"]
    def __init__(self):
        self.setPath()

    def setPath(self):
        self.lng = array([ 24.096,  24.05 ,  23.956,  19.85 ,  19.71 ,  15.49 ,  15.319,
                           11.064,  10.854,  10.791,   6.421,   6.337,   1.875,   1.755,
                           -2.539,  -2.688,  -6.856,  -7.026, -11.244, -11.451, -11.511,
                          -15.652, -15.736, -19.872, -19.96 , -23.996, -24.116, -28.322,
                          -28.456, -32.428, -32.607, -32.645, -37.16 , -37.236, -41.697,
                          -41.812, -46.403, -46.548, -51.347, -51.528, -56.219, -56.431,
                          -56.488, -61.367, -61.466, -66.373, -66.495, -71.29 , -71.45 ,
                          -76.601, -76.782, -82.773, -83.075, -83.125, -83.125, -83.075,
                          -73.511, -73.039, -66.349, -66.047, -60.667, -60.412, -55.567,
                          -55.415, -50.726, -50.539, -45.782, -45.679, -40.977, -40.786,
                          -36.628, -36.457, -32.083, -31.95 , -27.794, -27.705, -23.552,
                          -23.484, -19.131, -19.085, -18.897, -14.844, -14.674, -10.473,
                          -10.353,  -6.228,  -6.161,  -1.858,  -1.775,   2.668,   2.718,
                            2.921,   6.932,   7.088,  11.425,  11.57 ,  16.01 ,  16.162,
                           20.183,  20.295])
        self.lat = array([-29.183, -28.037, -27.856, -25.883, -25.766, -24.429, -24.306,
                          -23.277, -23.191, -23.154, -21.941, -21.856, -21.23 , -21.169,
                          -19.926, -19.853, -18.171, -18.046, -16.763, -16.668, -16.63 ,
                          -14.202, -14.126, -11.829, -11.613,  -8.967,  -8.759,  -7.137,
                           -6.841,  -4.202,  -3.975,  -3.832,  -3.19 ,  -3.084,  -2.362,
                           -2.296,  -2.427,  -2.354,  -3.747,  -3.703,  -4.729,  -4.68 ,
                           -4.645,  -5.859,  -5.882,  -7.422,  -7.4  ,  -8.56 ,  -8.57 ,
                          -11.701, -11.665, -18.999, -19.382, -19.325, -19.325, -19.382,
                          -37.638, -39.078, -45.589, -45.885, -48.384, -48.613, -49.416,
                          -49.446, -49.716, -49.92 , -50.265, -50.333, -50.315, -50.276,
                          -49.349, -49.34 , -48.872, -48.846, -47.652, -47.59 , -46.227,
                          -46.194, -45.266, -45.263, -45.202, -43.788, -43.77 , -42.66 ,
                          -42.581, -41.051, -40.897, -39.792, -39.802, -39.038, -39.041,
                          -39.019, -37.222, -37.142, -36.369, -36.366, -35.83 , -35.976,
                          -33.654, -33.772])
        self.vertices = array((self.lng, self.lat)).transpose()
        self.path = path.Path(self.vertices)

    def contains(self, lat, lng):
        if (lng > 24.096 or lng < -83.125 or lat > -2.296 or lat < -50.333):  
            inside = False
        elif (lng > -66.5 and lng < -19.8 and lat > -45.4 and lat < -11.85):
            inside = True
        else:
            inside = self.path.contains_point((lng, lat))
        return inside
