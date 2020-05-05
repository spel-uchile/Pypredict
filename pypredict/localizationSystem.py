'''
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
'''
from pypredict.ekf import EKF
from pypredict.tdoa import TDOA
from numpy import cos, matrix, pi, random, sin, transpose

class Loc(object):

    def __init__(self):
        self.ekf = EKF()
        self.tdoa = TDOA()
        self.z = matrix([[0], [0], [0], [0], [0], [0]])

    def __call__(self):
        return self

    def __str__(self):
        return "This localization system combines a model of an orbit obtained at deploy time with TDOA measurements. This is done by using an Extended Kalman Filter (EKF)."

    def initTransformationData(self, i, RAAN, theta, w, r):
        """
        Initializes the rotation matrices and a translation vector.

        Parameters
        ----------
        i: float
           The inclination of the main satellite, the one which is
           at the center of the orbital frame, in radians.
        RAAN: float
              The right ascention of the ascending node of the
              main satellite in radians.
        theta: float
               The true anomaly of the main satellite in radians.
        w: float
           The argument of the perigee of the main satellite in
           radians.
        r: matrix
           The position vector of the main satellite in meters.
        """
        sin_RAAN = sin(RAAN)
        cos_RAAN = cos(RAAN)
        sin_i = sin(i)
        cos_i = cos(i)
        sin_w = sin(w)
        cos_w = cos(w)
        self.Q_ip = matrix([[-sin_RAAN*cos_i*sin_w + cos_RAAN*cos_w,
                             cos_RAAN*cos_i*sin_w + sin_RAAN*cos_w,
                             sin_i*sin_w],
                            [-sin_RAAN*cos_i*cos_w - cos_RAAN*sin_w,
                             cos_RAAN*cos_i*cos_w - sin_RAAN*sin_w,
                             sin_i*cos_w],
                            [sin_RAAN*sin_i, -cos_RAAN*sin_i, cos_i]])
        phi = theta + pi/2
        self.Q_po = matrix([[cos(phi), sin(phi), 0],
                            [0, 0, -1],
                            [-sin(phi), cos(phi), 0]])
        self.translation = self.Q_po*self.Q_ip*r

    def transformPosition(self, r):
        """
        Transforms the position vector r from the inertial frame to
        the perifocal frame, and then from the perifocal frame to the
        orbital frame.

        Parameters
        ----------
        r: matrix
           Position vector in the inertial frame.
        """
        return self.Q_po*self.Q_ip*r - self.translation

    def transformVelocity(self, v):
        """
        Transforms the velocity from the inertial frame to the
        perifocal frame, and then from the perifocal frame to the
        orbital frame.

        Parameters
        ----------
        v: matrix
           Velocity in the inertial frame.
        """
        v = self.Q_po*self.Q_ip*v
        v_orb = v - self.mainSat_v
        return v_orb

    def orbital2inertial(self, r):
        """
        Transforms the vector r from the orbital frame to the
        perifocal frame, and then from the perifocal frame to the
        inertial frame.

        Parameters
        ----------
        r: matrix
           Position vector in the orbital frame.
        """
        r_in = transpose(self.Q_ip)*transpose(self.Q_po)*(r + self.translation)
        return r_in
    
    def estimateLocation(self, i, RAAN, theta, w, v, dt, n, r0_i, r1_i, r2_i, r3_i, r4_i, x, z=None):
        self.initTransformationData(i, RAAN, theta, w, r1_i)
        if (z is None):
            r0_o = self.transformPosition(r0_i)
            r1_o = self.transformPosition(r1_i)
            r2_o = self.transformPosition(r2_i)
            r3_o = self.transformPosition(r3_i)
            r4_o = self.transformPosition(r4_i)
            self.z[0:3] = self.tdoa.calculateLocation(r0_o, r1_o,
                                                      r2_o, r3_o, r4_o)
        else:
            self.z[0:3] = self.transformPosition(z)
        # Adding noise
        self.z[0] = random.normal(self.z[0], 3.2)
        self.z[1] = random.normal(self.z[1], 3.2)
        self.z[2] = random.normal(self.z[2], 3.2)
        # End
        self.mainSat_v = self.Q_po*self.Q_ip*v
        x[0:3] = self.transformPosition(x[0:3])
        x[3:6] = self.transformVelocity(x[3:6])
        self.ekf.newIteration(dt, n, x, self.z[0:3])

    def getEstimatedPos(self):
        """
        Returns the estimated position vector in the inertial frame.
        """
        r_orb = self.ekf.getPosition()
        return self.orbital2inertial(r_orb)

    def getEstimatedVel(self):
        """
        Returns the estimated velocity vector in the inertial frame.
        """
        v_orb = self.ekf.getVelocity() + self.mainSat_v
        v_in = transpose(self.Q_ip)*transpose(self.Q_po)*v_orb
        return v_in
