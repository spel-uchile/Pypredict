"""
                                Pypredict
    Orbit prediction software. Displays the satellites' position and
    orbital parameters in real time. Simulates satellite localization
    and deployment.
    
    Copyright (C) 2018-2021, Matías Vidal Valladares, matvidal.
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

from AOAwithTDOA import Locate
from pypredict.dpl import Dpl
from pypredict.ekf import EKF
from pypredict.sat import Sat
from pypredict.localizationSystem import Loc
from numpy import arange, arctan2, arccos, array, cos, linalg, matrix, mean, pi, random, savetxt, sin, std, sqrt, zeros, transpose
from datetime import datetime, timedelta
from pkg_resources import resource_filename


class Loc(object):

    def __init__(self):
        self.ekf = EKF()
        self.loc = Locate()
        self.z = matrix([[0], [0], [0], [0], [0], [0]])

    def __call__(self):
        return self

    def __str__(self):
        return "This localization system combines a model of an orbit obtained at deploy time with TDOA and AOA measurements. This is done by using an Extended Kalman Filter (EKF)."

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

    def generate_noise(self, std_RD, std_AOA, std_GNSS, L):
        e = self.loc.get_error_vector(std_RD, std_AOA, L)
        GNSS_noise = self.loc.get_GNSS_noise(std_GNSS, L)
        return e, GNSS_noise

    def estimate_before_ekf(self, e, GNSS_noise, s1, s2, u, Q):
        noisy_s1, noisy_s2 = self.loc.add_GNSS_error(s1, s2, GNSS_noise)
        k_w_GNSS_error = self.loc.get_real_vector(u, noisy_s1, noisy_s2)
        u_hat = self.loc.estimate(noisy_s1, noisy_s2, k_w_GNSS_error, e, Q)
        return noisy_s1, noisy_s2, u_hat

    def estimate_with_ekf(self, v, dt, n, rc, x, s1, s2, z):
        self.mainSat_v = self.Q_po*self.Q_ip*v
        x[0:3] = self.transformPosition(x[0:3])
        x[3:6] = self.transformVelocity(x[3:6])
        self.ekf.newIteration(dt, n, rc, x, s1, s2, z)

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

    def test_case(self):
        date0 = datetime(2020, 11, 17, 00, 13, 33, 0)
        dep_date = date0 + timedelta(minutes=10)
        v = [0, 1, 0]
        x = matrix([[0], [0], [0], [0], [0], [0]])
        std_ADS = 36/3600 # STT of 36 arcseconds.
        std_ACS = 0.06    # RW of 0.06°.
        std_RD = 10.0
        std_AOA = 1.0*pi/180.0
        std_GNSS = 10.0
        dpl = Dpl()
        s1_mass = 3.2
        u_mass = 0.08
        L = 5000
        minutes = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
        rmse = zeros(len(minutes))
        bias = zeros(len(minutes))
        div = 32
        finer_minutes = arange(0, div*100+1, dtype="float")/div
        rcrb = zeros(len(finer_minutes))
        dist_u_s1 = zeros(len(finer_minutes))
        dist_u_s2 = zeros(len(finer_minutes))
        alpha_max = zeros(len(finer_minutes))
        alpha_mean = zeros(len(finer_minutes))
        alpha_min = zeros(len(finer_minutes)) + 180
        data_path = resource_filename("pypredict","data/")
        sat_s2 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        L0 = 1
        dep_noise = self.loc.get_deployment_noise(std_ADS, std_ACS, L0)
        studied_date = dep_date + timedelta(days=3)
        sat_s1 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s1.updateOrbitalParameters(dep_date)
        i = 0
        vel = self.loc.noisy_dep_velocity(v, dep_noise[:,i])
        sat_u = dpl.deploy("Femto", sat_s1, s1_mass, u_mass, "FE1", vel, dep_date)
        sat_u.updateOrbitalParameters(studied_date)
        sat_s1.updateOrbitalParameters(studied_date)
        sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
        e, GNSS_noise = self.generate_noise(std_RD, std_AOA, std_GNSS, L)
        u = sat_u.getXYZ()
        s1 = sat_s1.getXYZ()
        s2 = sat_s2.getXYZ()
        Q = self.loc.get_Q(std_RD, std_AOA)
        noisy_s1, noisy_s2, u_hat = self.estimate_before_ekf(e[:,1], GNSS_noise[:,1], s1, s2, u, Q)
        self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, s1)
        dt = 10
        x[0:3] = s1
        x[3:6] = sat_s1.getInertialVel()
        rc = sat_s1.a*sqrt(1 - sat_s1.e**2)
        self.estimate_with_ekf(sat_s1.getInertialVel(), dt, sat_s1.n, rc, x, noisy_s1, noisy_s2, u_hat)

