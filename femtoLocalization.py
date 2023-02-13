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
from matplotlib.pyplot import axis, figure, plot, savefig, setp, show, subplots, subplots_adjust, tight_layout, xticks
from datetime import datetime, timedelta
from pkg_resources import resource_filename
import csv


class Loc(object):

    def __init__(self):
        self.ekf = EKF()
        self.loc = Locate()
        self.z = matrix([[0.0], [0.0], [0.0], [0.0], [0.0], [0.0]])
        self.Q_ip = matrix([[0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0]])
        self.Q_po = matrix([[0.0, 0.0,  0.0],
                            [0.0, 0.0, -1.0],
                            [0.0, 0.0,  0.0]])
        self.dpl = Dpl()
        #self.test_case()

    def __call__(self):
        return self

    def __str__(self):
        return "This localization system combines a model of an orbit obtained at deploy time with TDOA and AOA measurements. This is done by using an Extended Kalman Filter (EKF)."

    def initTransformationData(self, i, RAAN, theta, w, r):
        """
        Initializes the rotation matrices and a translation vector.

        Parameters
        ----------
        i     : float
                The inclination of the main satellite, the one which is
                at the center of the orbital frame, in radians.
        RAAN  : float
                The right ascention of the ascending node of the
                main satellite in radians.
        theta : float
                The true anomaly of the main satellite in radians.
        w     : float
                The argument of the perigee of the main satellite in
                radians.
        r     : matrix
                The position vector of the main satellite in meters.
        """
        sin_RAAN = sin(RAAN)
        cos_RAAN = cos(RAAN)
        sin_i = sin(i)
        cos_i = cos(i)
        sin_w = sin(w)
        cos_w = cos(w)
        #self.Q_ip = matrix([[-sin_RAAN*cos_i*sin_w + cos_RAAN*cos_w,
        #                     cos_RAAN*cos_i*sin_w + sin_RAAN*cos_w,
        #                     sin_i*sin_w],
        #                    [-sin_RAAN*cos_i*cos_w - cos_RAAN*sin_w,
        #                     cos_RAAN*cos_i*cos_w - sin_RAAN*sin_w,
        #                     sin_i*cos_w],
        #                    [sin_RAAN*sin_i, -cos_RAAN*sin_i, cos_i]])
        #print("A: {}\n".format(self.Q_ip))
        self.Q_ip[0,0] = -sin_RAAN*cos_i*sin_w + cos_RAAN*cos_w
        self.Q_ip[0,1] = cos_RAAN*cos_i*sin_w + sin_RAAN*cos_w
        self.Q_ip[0,2] = sin_i*sin_w
        self.Q_ip[1,0] = -sin_RAAN*cos_i*cos_w - cos_RAAN*sin_w
        self.Q_ip[1,1] = cos_RAAN*cos_i*cos_w - sin_RAAN*sin_w
        self.Q_ip[1,2] = sin_i*cos_w
        self.Q_ip[2,0] = sin_RAAN*sin_i
        self.Q_ip[2,1] = -cos_RAAN*sin_i
        self.Q_ip[2,2] = cos_i
        #print("B: {}\n".format(self.Q_ip))
        phi = theta + pi*0.5
        cos_phi = cos(phi)
        sin_phi = sin(phi)
        #self.Q_po = matrix([[cos_phi, sin_phi, 0],
        #                    [0, 0, -1],
        #                    [-sin_phi, cos_phi, 0]])
        self.Q_po[0,0] = cos_phi
        self.Q_po[0,1] = sin_phi
        self.Q_po[2,0] = -sin_phi
        self.Q_po[2,1] = cos_phi
        self.Q_poQ_ip = self.Q_po*self.Q_ip
        self.translation = self.Q_poQ_ip*r

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
        return self.Q_poQ_ip*r - self.translation

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
        v = self.Q_poQ_ip*v
        return v - self.mainSat_v

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
        #r_in = transpose(self.Q_ip)*transpose(self.Q_po)*(r + self.translation)
        r_in = (self.Q_po*self.Q_ip).T*(r + self.translation)
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
        self.mainSat_v = self.Q_poQ_ip*v
        #x[0:3] = self.transformPosition(x[0:3])
        #x[3:6] = self.transformVelocity(x[3:6])
        lvlh_s1 = self.transformPosition(s1)
        lvlh_s2 = self.transformPosition(s2)
        lvlh_z = self.transformPosition(z)
        #print("z: {}".format(lvlh_z))
        self.ekf.newIteration(dt, n, rc, x, lvlh_s1, lvlh_s2, lvlh_z)

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
        #v_in = transpose(self.Q_ip)*transpose(self.Q_po)*v_orb
        v_in = (self.Q_po*self.Q_ip).T*v_orb
        return v_in

    def test_case(self):
        start = datetime.utcnow()
        date0 = datetime(2020, 11, 17, 00, 13, 33, 0)
        dep_date = date0 + timedelta(minutes=10)
        x = matrix([[0], [0], [0], [0], [0], [0]])
        std_RD = 10.0
        std_AOA = 1.0*pi/180.0
        std_GNSS = 10.0
        Q = self.loc.get_Q(std_RD, std_AOA)
        L = 2#0#00 #Con 50 -> +9 hrs.
        L0 = 2#000
        minutes = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
        acum_ekf = zeros((len(minutes), 3, 1))
        real_positions = zeros((len(minutes), 3, 1))
        rmse = zeros(len(minutes))
        bias = zeros(len(minutes))
        meas_error = zeros(len(minutes))
        div = 32
        finer_minutes = arange(0, div*100+1, dtype="float")/div
        dt = 2#5
        ekf_time = arange(0, int(100*60/dt)+1, dtype="float")
        rcrb = zeros(len(finer_minutes))
        ekf_error = zeros(len(ekf_time))
        dist_u_s1 = zeros(len(finer_minutes))
        dist_u_s2 = zeros(len(finer_minutes))
        alpha_max = zeros(len(finer_minutes))
        alpha_mean = zeros(len(finer_minutes))
        alpha_min = zeros(len(finer_minutes)) + 180
        data_path = resource_filename("pypredict","data/")
        sat_u = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s1 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2.setTLE("1 44814U 19081L   20321.73053029  .00001305  00000-0  63025-4 0  9996",
                      "2 44814  97.4788  21.6285 0013387  80.2501 280.0246 15.20374749 54001")
        with open("u_TLEs.txt", 'r') as u_TLEs:
            sat_u_TLEs = u_TLEs.readlines()
        with open("s1_TLEs.txt", 'r') as s1_TLEs:
            sat_s1_TLEs = s1_TLEs.readlines()
        for i in range(L0):
            print("Deployment {}/{}".format(i+1, L0))
            #studied_date = dep_date + timedelta(days=3)
            #j = 0
            #e, GNSS_noise = self.generate_noise(std_RD, std_AOA, std_GNSS, L)
            sat_u.setTLE(sat_u_TLEs[1 + i*3], sat_u_TLEs[2 + i*3])
            sat_s1.setTLE(sat_s1_TLEs[1 + i*3], sat_s1_TLEs[2 + i*3])
            acum_ekf = acum_ekf*0
            self.ekf.initMatrices()
            #sat_u.updateOrbitalParameters(studied_date)
            #sat_s1.updateOrbitalParameters(studied_date)
            #sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
            #u = sat_u.getXYZ()
            #s1 = sat_s1.getXYZ()
            #s2 = sat_s2.getXYZ()
            #k = self.loc.get_real_vector(u, s1, s2)
            #MSE = self.loc.get_MSE(u, s1, s2, k, Q)
            #rcrb[j] += MSE[0,0] + MSE[1,1] + MSE[2,2]
            #noisy_s1, noisy_s2, u_hat = self.estimate_before_ekf(e[:,j], GNSS_noise[:,j], s1, s2, u, Q)
            #self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, s1)
            #x[0:3] = self.transformPosition(s1)#self.transformPosition(u_hat.copy()) #noisy_s1
            #x[3:6] = 0 #sat_s1self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, noisy_s1).getInertialVel()
            #meas_error[0] += linalg.norm(u - u_hat)
            #ekf_error[0] += linalg.norm(u - u_hat)
            #print("Measurement error: {}".format(meas_error[0]))
            for j in range(L):
                studied_date = dep_date + timedelta(days=3)
                sat_s1.updateOrbitalParameters(studied_date - timedelta(seconds=dt))
                s1 = sat_s1.getXYZ().copy()
                self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, s1)
                x[0:3] = self.transformPosition(s1)#self.transformPosition(u_hat.copy()) #noisy_s1
                x[3:6] = 0 #sat_s1.getInertialVel()
                rc = sat_s1.a*sqrt(1 - sat_s1.e**2)
                self.ekf.Predict(dt, sat_s1.n, rc, x)
                x = self.ekf.x_pred.copy()
                e, GNSS_noise = self.generate_noise(std_RD, std_AOA, std_GNSS, len(minutes))
                for t in range(int(100*60/dt)+1):
                    #studied_date += timedelta(seconds=dt)
                    sat_u.updateOrbitalParameters(studied_date)
                    sat_s1.updateOrbitalParameters(studied_date)
                    s1 = sat_s1.getXYZ().copy()
                    u = sat_u.getXYZ().copy()
                    self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, s1)
                    lvlh_u = self.transformPosition(u).copy()
                    rc = sat_s1.a*sqrt(1 - sat_s1.e**2)
                    #print(rc)
                    #self.ekf.Predict(dt, sat_s1.n, rc, x)
                    #x = self.ekf.x_pred.copy()
                    #if (j == 0):
                        #sat_u.updateOrbitalParameters(studied_date)
                        #sat_s1.updateOrbitalParameters(studied_date)
                        #sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
                        #u = sat_u.getXYZ()
                        #s1 = sat_s1.getXYZ()
                        #s2 = sat_s2.getXYZ()
                        #k = self.loc.get_real_vector(u, s1, s2)
                        #MSE = self.loc.get_MSE(u, s1, s2, k, Q)
                        #rcrb[t] += MSE[0,0] + MSE[1,1] + MSE[2,2]
                    if (t%(300/dt) == 0):
                        #ekf_estimations = []
                        #acumulated_meas = 0
                        #acumulated_ekf = 0
                        #sat_u.updateOrbitalParameters(studied_date)
                        #sat_s1.updateOrbitalParameters(studied_date)
                        sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
                        #u = sat_u.getXYZ()
                        #s1 = sat_s1.getXYZ().copy()
                        s2 = sat_s2.getXYZ().copy()
                        index = int(t/(300/dt))
                        #rc = sat_s1.a*sqrt(1 - sat_s1.e**2)
                        noisy_s1, noisy_s2, u_hat = self.estimate_before_ekf(e[:,index], GNSS_noise[:,index], s1, s2, u, Q)
                        #self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, noisy_s1)
                        #self.estimate_with_ekf(sat_s1.getInertialVel(), dt, sat_s1.n, rc, x, noisy_s1, noisy_s2, u_hat)
                        #self.mainSat_v = self.Q_poQ_ip*sat_s1.getInertialVel()
                        lvlh_s1 = self.transformPosition(noisy_s1)
                        lvlh_s2 = self.transformPosition(noisy_s2)
                        lvlh_z = self.transformPosition(u_hat)
                        self.ekf.Update(lvlh_s1, lvlh_s2, lvlh_z)
                        #acumulated_meas += linalg.norm(u - u_hat)
                        #print("Measurement error: {}".format(meas_error[j]))
                        #acumulated_ekf += linalg.norm(self.transformPosition(u) - self.ekf.x_up[0:3])
                        #print("{}.- EKF error: {}".format(j, ekf_error[j-1]))
                        x = self.ekf.x_up.copy()
                        #ekf_estimations.append(x[0:3])
                        #lvlh_u = self.transformPosition(u)
                        meas_error[index] += linalg.norm(u - u_hat)
                        #ekf_error[index] += linalg.norm(lvlh_u - x[0:3])
                        rmse[index] += linalg.norm(lvlh_u - x[0:3])**2
                        acum_ekf[index] += x[0:3]
                        if (j == 0):
                            real_positions[index] = lvlh_u
                        #self.loc.RMSE(self.transformPosition(u), ekf_estimations, L)
                        #bias[int(t/30)] += self.loc.Bias(self.transformPosition(u), ekf_estimations, L)
                        #print("RMSE: {}, Bias: {}".format(rmse[int(j/30)], bias[int(j/30)]))
                        #x = self.ekf.x_up.copy()
                    ekf_error[t] += linalg.norm(lvlh_u - x[0:3])
                    studied_date += timedelta(seconds=dt)
                    self.ekf.Predict(dt, sat_s1.n, rc, x)
                    x = self.ekf.x_pred.copy()
                #meas_error[int(t/30)] += acumulated_meas/L
                #ekf_error[int(t/30)] += acumulated_ekf/L
                #rmse[int(t/30)] += self.loc.RMSE(self.transformPosition(u), ekf_estimations, L)
            #bias[int(t/30)] += self.loc.Bias(self.transformPosition(u), ekf_estimations, L)
            for en, acum in enumerate(acum_ekf): bias[en] += linalg.norm(acum/L - real_positions[en])**2
            #print(bias)
        meas_error = meas_error/L0/L
        ekf_error = ekf_error/L0/L
        rmse = sqrt(rmse/L0/L)
        bias = sqrt(bias/L0)
        with open("rcrb.csv", newline='') as f:
            reader = csv.reader(f, delimiter=",")
            next(reader)
            data = array(list(reader))#[:,0]
        for i in range(len(data)):
            rcrb[i] = float(data[i,0])
        #rcrb = sqrt(rcrb/L0)
        fig, ax = subplots(1, 1, sharey=False)
        ax.semilogy(minutes, rmse, 'o', linewidth=2.0, markersize=12,# clip_on=False,
                    fillstyle="none", label="RMSE")
        ax.semilogy(finer_minutes, rcrb, '-', linewidth=2.0, markersize=12,
                    label="Root CRB", color="tab:orange", alpha=0.7)
        ax.semilogy(minutes, bias, '+', linewidth=2.0, markersize=12,# clip_on=False,
                label="Bias", color="tab:green")
        ax.semilogy(minutes, meas_error, '-', linewidth=1.8, markersize=12,
                    label="Measurement error", color="tab:purple")
        ax.semilogy(ekf_time*dt/60, ekf_error, '-', linewidth=1.8, markersize=12,
                    label="EKF error", color="tab:red")
        ax.grid()
        ax.legend(fontsize=14, loc="upper left")
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        tight_layout()
        finish = datetime.utcnow()
        print("Elapsed time: {}".format(finish-start))
        show()

    def test_case2(self):
        std_RD = 10.0
        std_AOA = 1.0*pi/180.0
        std_GNSS = 10.0
        Q = self.loc.get_Q(std_RD, std_AOA)
        e_rcrb = zeros(42)
        bias = zeros(42)
        acum_est = zeros([42, 3, 1])
        real_u = zeros([42, 3, 1])
        real_s1 = zeros([42, 3, 1])
        real_s2 = zeros([42, 3, 1])
        data_path = resource_filename("pypredict","data/")
        sat_u = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s1 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2.setTLE("1 44814U 19081L   20321.73053029  .00001305  00000-0  63025-4 0  9996",
                      "2 44814  97.4788  21.6285 0013387  80.2501 280.0246 15.20374749 54001")
        with open("u_TLEs.txt", 'r') as u_TLEs:
            sat_u_TLEs = u_TLEs.readlines()

        with open("s1_TLEs.txt", 'r') as s1_TLEs:
            sat_s1_TLEs = s1_TLEs.readlines()

        date0 = datetime(2020, 11, 17, 00, 13, 33, 0)
        dep_date = date0 + timedelta(minutes=10)
        initial_date = dep_date + timedelta(days=3)
        L0 = [3, 15, 76, 309, 425]     # Worst cases
        ekf = EKF()
        meas_dt = 2.0#0.2
        model_dt = 0.1
        x = zeros([6,1])
        iterations = 1
        h = model_dt/iterations
        start = datetime.utcnow()
        minutes = zeros(42)
        for i in range(21):
            minutes[2*i] = i*5.0 - meas_dt/60
            minutes[2*i + 1] = i*5.0
        rmse = zeros(int(100*60/model_dt + 2))
        L = 50#0#5000
        for l0 in L0:
            print("Despliegue: {}".format(l0))
            sat_u.setTLE(sat_u_TLEs[1 + l0*3], sat_u_TLEs[2 + l0*3])
            sat_s1.setTLE(sat_s1_TLEs[1 + l0*3], sat_s1_TLEs[2 + l0*3])
            acum_est = acum_est*0
            for l in range(L):
                studied_date = initial_date
                studied_date -= timedelta(seconds=2*meas_dt)
                d = 0
                e, GNSS_noise = self.generate_noise(std_RD, std_AOA, std_GNSS, int(6000/model_dt + 2))
                for t in range(21):
                    studied_date += timedelta(seconds=meas_dt)
                    if (l==0):
                        sat_u.updateOrbitalParameters(studied_date)
                        sat_s1.updateOrbitalParameters(studied_date)
                        sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
                        real_u[2*t] = sat_u.getXYZ().copy()
                        real_s1[2*t] = sat_s1.getXYZ().copy()
                        real_s2[2*t] = sat_s2.getXYZ().copy()
                    noisy_s1, noisy_s2 = self.loc.add_GNSS_error(real_s1[2*t], real_s2[2*t], GNSS_noise[:,d])
                    k_w_GNSS_error = self.loc.get_real_vector(real_u[2*t], noisy_s1, noisy_s2)
                    u_hat1 = self.loc.estimate(noisy_s1, noisy_s2, k_w_GNSS_error, e[:,d], Q)
                    acum_est[2*t] += u_hat1
                    MSE = self.loc.get_MSE(u_hat1, noisy_s1, noisy_s2, k_w_GNSS_error, Q)
                    e_rcrb[2*t] += MSE[0,0] + MSE[1,1] + MSE[2,2]
                    rmse[d] += (u_hat1[0,0] - real_u[2*t,0,0])**2 + (u_hat1[1,0] - real_u[2*t,1,0])**2 + (u_hat1[2,0] - real_u[2*t,2,0])**2
                    d += 1
                    studied_date += timedelta(seconds=meas_dt)
                    sat_s1.updateOrbitalParameters(studied_date)
                    if (l==0):
                        sat_u.updateOrbitalParameters(studied_date)
                        sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
                        real_u[2*t+1] = sat_u.getXYZ().copy()
                        real_s1[2*t+1] = sat_s1.getXYZ().copy()
                        real_s2[2*t+1] = sat_s2.getXYZ().copy()
                    noisy_s1, noisy_s2 = self.loc.add_GNSS_error(real_s1[2*t+1], real_s2[2*t+1], GNSS_noise[:,d])
                    k_w_GNSS_error = self.loc.get_real_vector(real_u[2*t+1], noisy_s1, noisy_s2)
                    u_hat2 = self.loc.estimate(noisy_s1, noisy_s2, k_w_GNSS_error, e[:,d], Q)
                    acum_est[2*t+1] += u_hat2
                    MSE = self.loc.get_MSE(u_hat2, noisy_s1, noisy_s2, k_w_GNSS_error, Q)
                    e_rcrb[2*t+1] += MSE[0,0] + MSE[1,1] + MSE[2,2]
                    vel = (u_hat2 - u_hat1)/meas_dt
                    self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, noisy_s1)
                    self.mainSat_v = self.Q_poQ_ip*sat_s1.getInertialVel().copy()
                    lvlh_u = self.transformPosition(u_hat2)
                    lvlh_v = self.transformVelocity(vel)
                    x[0:3] = lvlh_u.copy()
                    x[3:6] = lvlh_v.copy()
                    rmse[d] += (u_hat2[0,0] - real_u[2*t+1,0,0])**2 + (u_hat2[1,0] - real_u[2*t+1,1,0])**2 + (u_hat2[2,0] - real_u[2*t+1,2,0])**2
                    d += 1
                    rc = sat_s1.a*sqrt(1 - sat_s1.e**2)
                    sat_s1_n = sat_s1.n
                    if ((studied_date - initial_date).total_seconds() < 6000):
                        for i in range(int(300/model_dt)-2):#161):
                            studied_date += timedelta(seconds=model_dt)
                            model_pos = ekf.RK4(x, h, iterations, sat_s1_n, rc)
                            sat_u.updateOrbitalParameters(studied_date)
                            sat_s1.updateOrbitalParameters(studied_date)
                            s1 = sat_s1.getXYZ().copy()
                            noisy_s1, noisy_s2 = self.loc.add_GNSS_error(s1, real_s2[2*t+1], GNSS_noise[:,d])
                            self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, noisy_s1)
                            real_pos = self.transformPosition(sat_u.getXYZ()).copy()
                            rmse[d] += (model_pos[0,0] - real_pos[0,0])**2 + (model_pos[1,0] - real_pos[1,0])**2 + (model_pos[2,0] - real_pos[2,0])**2
                            d += 1
                            x = model_pos.copy()
                        studied_date += timedelta(seconds=(model_dt - meas_dt)*2)
            acum_est = acum_est/L
            for j, u in enumerate(real_u):
                bias[j] += (acum_est[j,0,0] - u[0,0])**2 + (acum_est[j,1,0] - u[1,0])**2 + (acum_est[j,2,0] - u[2,0])**2

        print("Time3: {}".format((studied_date - initial_date).total_seconds()))

        e_rcrb = sqrt(e_rcrb/len(L0)/L)
        rmse = sqrt(rmse/len(L0)/L)
        bias = sqrt(bias/len(L0))
        with open("e_rcrb.csv", 'w') as data_file_rcrb:
            outRCRB = csv.writer(data_file_rcrb)
            outRCRB.writerow(["RMSE", "e RCRB", "Bias"])
            for i in range(len(e_rcrb)):
                outRCRB.writerow([rmse[i], e_rcrb[i], bias[i]])
            for i in range(len(e_rcrb), len(rmse)):
                outRCRB.writerow([rmse[i]])
        finish = datetime.utcnow()
        print("Start: {}\nFinish: {}\nDelta: {} minutes".format(start, finish, (finish-start).total_seconds()/60))

    def test_case3(self, L0=[582, 2989, 963, 3714, 425], L=5000, output="rmse.csv"):
        std_RD = 10.0
        std_AOA = 1.0*pi/180.0
        std_GNSS = 10.0
        Q = self.loc.get_Q(std_RD, std_AOA)
        data_path = resource_filename("pypredict","data/")
        sat_u = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s1 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2.setTLE("1 44814U 19081L   20321.73053029  .00001305  00000-0  63025-4 0  9996",
                      "2 44814  97.4788  21.6285 0013387  80.2501 280.0246 15.20374749 54001")
        with open("u_TLEs.txt", 'r') as u_TLEs:
            sat_u_TLEs = u_TLEs.readlines()

        with open("s1_TLEs.txt", 'r') as s1_TLEs:
            sat_s1_TLEs = s1_TLEs.readlines()

        date0 = datetime(2020, 11, 17, 00, 13, 33, 0)
        dep_date = date0 + timedelta(minutes=10)
        initial_date = dep_date + timedelta(days=3)
        ekf = EKF()
        meas_dt = 2.0#0.2
        model_dt = 0.01
        x = zeros([6,1])
        iterations = 1
        h = model_dt/iterations
        start = datetime.utcnow()
        minutes = zeros(26)
        for i in range(13):
            minutes[2*i] = (i*300.0 - meas_dt + 600.0)/60.0
            minutes[2*i + 1] = (i*300.0 + 600.0)/60
        e_rcrb = zeros(len(minutes))
        bias = zeros(len(minutes))
        rmse = zeros(len(minutes))
        acum_est = zeros([len(minutes), 3, 1])
        real_u = zeros([len(minutes), 3, 1])
        real_s1 = zeros([len(minutes), 3, 1])
        real_s2 = zeros([len(minutes), 3, 1])
        model_rmse = zeros(2*int(300/model_dt) - 2)
        for l0 in L0:
            print("Despliegue: {}".format(l0))
            sat_u.setTLE(sat_u_TLEs[1 + l0*3], sat_u_TLEs[2 + l0*3])
            sat_s1.setTLE(sat_s1_TLEs[1 + l0*3], sat_s1_TLEs[2 + l0*3])
            acum_est = acum_est*0
            for l in range(L):
                studied_date = initial_date + timedelta(minutes=10)
                studied_date -= timedelta(seconds=2*meas_dt)
                d = 0
                e, GNSS_noise = self.generate_noise(std_RD, std_AOA, std_GNSS, int(6000/model_dt + 2))
                for t in range(13):
                    studied_date += timedelta(seconds=meas_dt)
                    if (l==0):
                        sat_u.updateOrbitalParameters(studied_date)
                        sat_s1.updateOrbitalParameters(studied_date)
                        sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
                        real_u[2*t] = sat_u.getXYZ().copy()
                        real_s1[2*t] = sat_s1.getXYZ().copy()
                        real_s2[2*t] = sat_s2.getXYZ().copy()
                    noisy_s1, noisy_s2 = self.loc.add_GNSS_error(real_s1[2*t], real_s2[2*t], GNSS_noise[:,d])
                    k_w_GNSS_error = self.loc.get_real_vector(real_u[2*t], noisy_s1, noisy_s2)
                    u_hat1 = self.loc.estimate(noisy_s1, noisy_s2, k_w_GNSS_error, e[:,d], Q)
                    acum_est[2*t] += u_hat1
                    MSE = self.loc.get_MSE(u_hat1, noisy_s1, noisy_s2, k_w_GNSS_error, Q)
                    e_rcrb[2*t] += MSE[0,0] + MSE[1,1] + MSE[2,2]
                    rmse[2*t] += (u_hat1[0,0] - real_u[2*t,0,0])**2 + (u_hat1[1,0] - real_u[2*t,1,0])**2 + (u_hat1[2,0] - real_u[2*t,2,0])**2
                    d += 1
                    studied_date += timedelta(seconds=meas_dt)
                    sat_s1.updateOrbitalParameters(studied_date)
                    if (l==0):
                        sat_u.updateOrbitalParameters(studied_date)
                        sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
                        real_u[2*t+1] = sat_u.getXYZ().copy()
                        real_s1[2*t+1] = sat_s1.getXYZ().copy()
                        real_s2[2*t+1] = sat_s2.getXYZ().copy()
                    noisy_s1, noisy_s2 = self.loc.add_GNSS_error(real_s1[2*t+1], real_s2[2*t+1], GNSS_noise[:,d])
                    k_w_GNSS_error = self.loc.get_real_vector(real_u[2*t+1], noisy_s1, noisy_s2)
                    u_hat2 = self.loc.estimate(noisy_s1, noisy_s2, k_w_GNSS_error, e[:,d], Q)
                    acum_est[2*t+1] += u_hat2
                    MSE = self.loc.get_MSE(u_hat2, noisy_s1, noisy_s2, k_w_GNSS_error, Q)
                    e_rcrb[2*t+1] += MSE[0,0] + MSE[1,1] + MSE[2,2]
                    vel = (u_hat2 - u_hat1)/meas_dt
                    self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, noisy_s1)
                    self.mainSat_v = self.Q_poQ_ip*sat_s1.getInertialVel().copy()
                    lvlh_u = self.transformPosition(u_hat2)
                    lvlh_v = self.transformVelocity(vel)
                    x[0:3] = lvlh_u.copy()
                    x[3:6] = lvlh_v.copy()
                    rmse[2*t+1] += (u_hat2[0,0] - real_u[2*t+1,0,0])**2 + (u_hat2[1,0] - real_u[2*t+1,1,0])**2 + (u_hat2[2,0] - real_u[2*t+1,2,0])**2
                    d += 1
                    seconds = (studied_date - initial_date).total_seconds()
                    date_before_model = studied_date
                    if (seconds == 2100):
                        rc = sat_s1.a*sqrt(1 - sat_s1.e**2)
                        sat_s1_n = sat_s1.n
                        for i in range(2*int(300/model_dt)-2):#161):
                            studied_date += timedelta(seconds=model_dt)
                            model_pos = ekf.RK4(x, h, iterations, sat_s1_n, rc)
                            sat_u.updateOrbitalParameters(studied_date)
                            sat_s1.updateOrbitalParameters(studied_date)
                            s1 = sat_s1.getXYZ().copy()
                            noisy_s1, noisy_s2 = self.loc.add_GNSS_error(s1, real_s2[2*t+1], GNSS_noise[:,d])
                            self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, noisy_s1)
                            real_pos = self.transformPosition(sat_u.getXYZ()).copy()
                            model_rmse[i] += (model_pos[0,0] - real_pos[0,0])**2 + (model_pos[1,0] - real_pos[1,0])**2 + (model_pos[2,0] - real_pos[2,0])**2
                            d += 1
                            x = model_pos.copy()
                        studied_date = date_before_model + timedelta(seconds=300-2*meas_dt)
                    else:
                        studied_date += timedelta(seconds=300-2*meas_dt)
            acum_est = acum_est/L
            for j, u in enumerate(real_u):
                bias[j] += (acum_est[j,0,0] - u[0,0])**2 + (acum_est[j,1,0] - u[1,0])**2 + (acum_est[j,2,0] - u[2,0])**2

        print("Time3: {}".format((studied_date - initial_date).total_seconds()))

        e_rcrb = sqrt(e_rcrb/len(L0)/L)
        rmse = sqrt(rmse/len(L0)/L)
        model_rmse = sqrt(model_rmse/len(L0)/L)
        bias = sqrt(bias/len(L0))
        with open(output, 'w') as data_file_rcrb:
            outRCRB = csv.writer(data_file_rcrb)
            outRCRB.writerow(["Model RMSE", "e RCRB", "Bias", "RMSE"])
            for i in range(len(e_rcrb)):
                outRCRB.writerow([model_rmse[i], e_rcrb[i], bias[i], rmse[i]])
            for i in range(len(e_rcrb), len(model_rmse)):
                outRCRB.writerow([model_rmse[i]])
        finish = datetime.utcnow()
        print("Start: {}\nFinish: {}\nDelta: {} minutes".format(start, finish, (finish-start).total_seconds()/60))

    def test_without_model(self, L0=[2824], L=5000, output="best_rmse.csv"):
        std_RD = 10.0
        std_AOA = 1.0*pi/180.0
        std_GNSS = 10.0
        Q = self.loc.get_Q(std_RD, std_AOA)
        data_path = resource_filename("pypredict","data/")
        sat_u = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s1 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2.setTLE("1 44814U 19081L   20321.73053029  .00001305  00000-0  63025-4 0  9996",
                      "2 44814  97.4788  21.6285 0013387  80.2501 280.0246 15.20374749 54001")
        with open("u_TLEs.txt", 'r') as u_TLEs:
            sat_u_TLEs = u_TLEs.readlines()

        with open("s1_TLEs.txt", 'r') as s1_TLEs:
            sat_s1_TLEs = s1_TLEs.readlines()

        date0 = datetime(2020, 11, 17, 00, 13, 33, 0)
        dep_date = date0 + timedelta(minutes=10)
        initial_date = dep_date + timedelta(days=3)
        ekf = EKF()
        meas_dt = 2.0#0.2
        x = zeros([6,1])
        iterations = 1
        start = datetime.utcnow()
        minutes = zeros(42)
        for i in range(21):
            minutes[2*i] = i*5.0 - meas_dt/60
            minutes[2*i + 1] = i*5.0
        len_minutes = len(minutes)
        e_rcrb = zeros(len_minutes)
        bias = zeros(len_minutes)
        acum_est = zeros([len_minutes, 3, 1])
        real_u = zeros([len_minutes, 3, 1])
        real_s1 = zeros([len_minutes, 3, 1])
        real_s2 = zeros([len_minutes, 3, 1])
        rmse = zeros(len_minutes)
        for l0 in L0:
            print("Despliegue: {}".format(l0))
            sat_u.setTLE(sat_u_TLEs[1 + l0*3], sat_u_TLEs[2 + l0*3])
            sat_s1.setTLE(sat_s1_TLEs[1 + l0*3], sat_s1_TLEs[2 + l0*3])
            acum_est = acum_est*0
            for l in range(L):
                studied_date = initial_date
                studied_date -= timedelta(seconds=2*meas_dt)
                d = 0
                e, GNSS_noise = self.generate_noise(std_RD, std_AOA, std_GNSS, len_minutes)
                for t in range(21):
                    studied_date += timedelta(seconds=meas_dt)
                    if (l==0):
                        sat_u.updateOrbitalParameters(studied_date)
                        sat_s1.updateOrbitalParameters(studied_date)
                        sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
                        real_u[2*t] = sat_u.getXYZ().copy()
                        real_s1[2*t] = sat_s1.getXYZ().copy()
                        real_s2[2*t] = sat_s2.getXYZ().copy()
                    noisy_s1, noisy_s2 = self.loc.add_GNSS_error(real_s1[2*t], real_s2[2*t], GNSS_noise[:,d])
                    k_w_GNSS_error = self.loc.get_real_vector(real_u[2*t], noisy_s1, noisy_s2)
                    u_hat1 = self.loc.estimate(noisy_s1, noisy_s2, k_w_GNSS_error, e[:,d], Q)
                    acum_est[2*t] += u_hat1
                    MSE = self.loc.get_MSE(u_hat1, noisy_s1, noisy_s2, k_w_GNSS_error, Q)
                    e_rcrb[2*t] += MSE[0,0] + MSE[1,1] + MSE[2,2]
                    rmse[2*t] += (u_hat1[0,0] - real_u[2*t,0,0])**2 + (u_hat1[1,0] - real_u[2*t,1,0])**2 + (u_hat1[2,0] - real_u[2*t,2,0])**2
                    d += 1
                    studied_date += timedelta(seconds=meas_dt)
                    sat_s1.updateOrbitalParameters(studied_date)
                    if (l==0):
                        sat_u.updateOrbitalParameters(studied_date)
                        sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
                        real_u[2*t+1] = sat_u.getXYZ().copy()
                        real_s1[2*t+1] = sat_s1.getXYZ().copy()
                        real_s2[2*t+1] = sat_s2.getXYZ().copy()
                    noisy_s1, noisy_s2 = self.loc.add_GNSS_error(real_s1[2*t+1], real_s2[2*t+1], GNSS_noise[:,d])
                    k_w_GNSS_error = self.loc.get_real_vector(real_u[2*t+1], noisy_s1, noisy_s2)
                    u_hat2 = self.loc.estimate(noisy_s1, noisy_s2, k_w_GNSS_error, e[:,d], Q)
                    acum_est[2*t+1] += u_hat2
                    MSE = self.loc.get_MSE(u_hat2, noisy_s1, noisy_s2, k_w_GNSS_error, Q)
                    e_rcrb[2*t+1] += MSE[0,0] + MSE[1,1] + MSE[2,2]
                    rmse[2*t+1] += (u_hat2[0,0] - real_u[2*t+1,0,0])**2 + (u_hat2[1,0] - real_u[2*t+1,1,0])**2 + (u_hat2[2,0] - real_u[2*t+1,2,0])**2
                    d += 1
                    studied_date += timedelta(seconds=300-2*meas_dt)
            acum_est = acum_est/L
            for j, u in enumerate(real_u):
                bias[j] += (acum_est[j,0,0] - u[0,0])**2 + (acum_est[j,1,0] - u[1,0])**2 + (acum_est[j,2,0] - u[2,0])**2

        print("Time3: {}".format((studied_date - initial_date).total_seconds()))

        e_rcrb = sqrt(e_rcrb/len(L0)/L)
        rmse = sqrt(rmse/len(L0)/L)
        bias = sqrt(bias/len(L0))
        with open(output, 'w') as data_file_rcrb:
            outRCRB = csv.writer(data_file_rcrb)
            outRCRB.writerow(["RMSE", "e RCRB", "Bias"])
            for i in range(len(e_rcrb)):
                outRCRB.writerow([rmse[i], e_rcrb[i], bias[i]])
            for i in range(len(e_rcrb), len(rmse)):
                outRCRB.writerow([rmse[i]])
        finish = datetime.utcnow()
        print("Start: {}\nFinish: {}\nDelta: {} minutes".format(start, finish, (finish-start).total_seconds()/60))

    def plot_data(self, rcrb_file="worst_rcrb.csv", rmse_file="rmse.csv", model=1):
        div = 32
        finer_minutes = arange(0, div*100+1, dtype="float")/div
        len_finer_minutes = len(finer_minutes)
        rcrb = zeros(len_finer_minutes)
        dist_u_s1 = zeros(len_finer_minutes)
        dist_u_s2 = zeros(len_finer_minutes)
        alpha_max = zeros(len_finer_minutes)
        alpha_mean = zeros(len_finer_minutes)
        alpha_min = zeros(len_finer_minutes)
        meas_dt = 2.0#0.13#0.2
        model_dt = 0.01
        model_rmse = zeros(2*int(300/model_dt) - 2)
        if (model):
            minutes = zeros(26)
            for i in range(13):
                minutes[2*i] = (i*300.0 - meas_dt + 600.0)/60.0
                minutes[2*i + 1] = (i*300.0 + 600.0)/60.0
            ind0 = 320
            indf = 2241
        else:
            minutes = zeros(42)
            for i in range(21):
                minutes[2*i] = i*5.0 - meas_dt/60
                minutes[2*i + 1] = i*5.0
            ind0 = 0
            indf = len_finer_minutes
        len_minutes = len(minutes)
        e_rcrb = zeros(len_minutes)
        bias = zeros(len_minutes)
        rmse = zeros(len_minutes)
        with open(rcrb_file, newline='') as f:
            reader = csv.reader(f, delimiter=",")
            next(reader)
            data = array(list(reader), dtype=object)
        for i in range(len(data)):
            rcrb[i] = float(data[i,0])
            dist_u_s1[i] = float(data[i,1])
            dist_u_s2[i] = float(data[i,2])
            alpha_min[i] = float(data[i,3])
            alpha_mean[i] = float(data[i,4])
            alpha_max[i] = float(data[i,5])
        with open(rmse_file, newline='') as f:
            reader = csv.reader(f, delimiter=",")
            next(reader)
            data = array(list(reader), dtype=object)
        if (model):
            for i in range(len(model_rmse)):
                model_rmse[i] = float(data[i][0])
            for i in range(len(e_rcrb)):
                e_rcrb[i] = float(data[i][1])
                bias[i] = float(data[i][2])
                rmse[i] = float(data[i][3])
        else:
            for i in range(len(e_rcrb)):
                rmse[i] = float(data[i][0])
                e_rcrb[i] = float(data[i][1])
                bias[i] = float(data[i][2])
        rmse_time = []
        time = 2100
        for j in range(2*int(300/model_dt) - 2):
            time += model_dt
            rmse_time.append(time/60)

        fig, ax = subplots(1, 1)
        if (model):
            ax.set_xlim(10, 70)
            ax.set_ylim(10, 1000000)
        else:
            ax.set_xlim(0,100)
            ax.set_ylim(1, 10000)
        ax2 = ax.twinx()
        ax2.set_ylim(0, 180)
        ax.set_zorder(10)
        ax.patch.set_visible(False)
        #ax2.plot(finer_minutes[ind0:indf], alpha_max[ind0:indf], '-', linewidth=2.0, markersize=12,
        #         label="{} (max)".format(r'$\alpha$'), color="k")
        ax2.plot(finer_minutes[ind0:indf], alpha_mean[ind0:indf], '-', linewidth=2.0, markersize=12,
                 label="{} (mean)".format(r'$\alpha$'), color="dimgrey")
        #ax2.plot(finer_minutes[ind0:indf], alpha_min[ind0:indf], '-', linewidth=2.0, markersize=12,
        #         label="{} (min)".format(r'$\alpha$'), color="lightgrey")
        ax.semilogy(finer_minutes[ind0:indf], rcrb[ind0:indf], '-', linewidth=2.0, markersize=12,
                    label="Root CRB", color="tab:orange")#, alpha=0.7)
        ax.semilogy(finer_minutes[ind0:indf], dist_u_s1[ind0:indf], '-', linewidth=2.0, markersize=12,
                    label=r'$||\mathbf{u} - \mathbf{s_1}||$', color="tab:purple", alpha=0.7)
        #ax.semilogy(finer_minutes, dist_u_s2, '-', linewidth=2.0, markersize=12,
        #            label=r'$||\mathbf{u} - \mathbf{s_2}||$', color="tab:red")
        #ax.semilogy(minutes, e_rcrb, '-', linewidth=2.0, markersize=12,
        #            label="Root CRB", color="tab:red", alpha=0.7)
        ax.semilogy(minutes, rmse, 'o', linewidth=2.0, markersize=12, fillstyle="none",
                    label="Measurement RMSE", color="tab:blue", alpha=0.7)
        if (model):
            ax.semilogy(rmse_time, model_rmse, '-', linewidth=2.0, markersize=12,
                        label="Model RMSE", color="tab:blue", alpha=0.7)
        ax.semilogy(minutes, bias, '+', linewidth=2.0, markersize=12,# clip_on=False,
                    label="Bias", color="tab:green")
        ax.grid()
        ax.set(xlabel="Time [min]",
               ylabel="RCRB, RMSE, bias and {} [m]".format(r'$||\mathbf{u} - \mathbf{s_1}||$'))
        ax2.set_ylabel("{} [deg]".format(r'$\alpha$'))
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax2.yaxis.label.set_size(16)
        ax.tick_params(which="both", direction="in", labelsize=14,
                       bottom=True, top=True, left=True, right=False)
        ax2.tick_params(which="both", direction="in", labelsize=14,
                        bottom=False, top=False, left=False, right=True)
        if (model):
            ax.set_xticks([10, 20, 30, 40, 50, 60, 70])
        else:
            ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        tight_layout()
        show()

    def plot_estimated_distance(self, l0=425, L=1, output="model.csv"):
        start = datetime.utcnow()
        std_RD = 10.0
        std_AOA = 1.0*pi/180.0
        std_GNSS = 10.0
        Q = self.loc.get_Q(std_RD, std_AOA)
        data_path = resource_filename("pypredict","data/")
        sat_u = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s1 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2.setTLE("1 44814U 19081L   20321.73053029  .00001305  00000-0  63025-4 0  9996",
                      "2 44814  97.4788  21.6285 0013387  80.2501 280.0246 15.20374749 54001")
        with open("u_TLEs.txt", 'r') as u_TLEs:
            sat_u_TLEs = u_TLEs.readlines()

        with open("s1_TLEs.txt", 'r') as s1_TLEs:
            sat_s1_TLEs = s1_TLEs.readlines()

        date0 = datetime(2020, 11, 17, 00, 13, 33, 0)
        dep_date = date0 + timedelta(minutes=10)
        initial_date = dep_date + timedelta(days=3)
        time = arange(0, 6000)#*2
        dist_u_s1 = zeros(len(time))
        dist_u_s2 = zeros(len(time))
        est_dist_u_s1 = zeros(len(time))
        #est_dist_u_s2 = zeros(len(time))
        dist_diff_between_meas = zeros(len(time)-1)
        mean_est_dist_u_s1 = zeros(len(time))
        #mean_est_dist_u_s2 = zeros(len(time))
        meas_error = zeros(len(time))
        mean_meas_error = zeros(len(time))
        mean_meas_and_model_error = zeros(len(time))
        u_hat = zeros([len(time), 3, 1])
        all_u_hat = zeros([len(time), 3, 1])
        all_u = zeros([len(time), 3, 1])
        all_s1 = zeros([len(time), 3, 1])
        all_s2 = zeros([len(time), 3, 1])
        sat_u.setTLE(sat_u_TLEs[1 + l0*3], sat_u_TLEs[2 + l0*3])
        sat_s1.setTLE(sat_s1_TLEs[1 + l0*3], sat_s1_TLEs[2 + l0*3])
        for i in range(len(time)):
            studied_date = initial_date + timedelta(seconds=float(time[i]))
            sat_u.updateOrbitalParameters(studied_date)
            sat_s1.updateOrbitalParameters(studied_date)
            sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
            all_u[i,:] = sat_u.getXYZ().copy()
            all_s1[i,:] = sat_s1.getXYZ().copy()
            all_s2[i,:] = sat_s2.getXYZ().copy()
            dist_u_s1[i] = linalg.norm(all_u[i,:] - all_s1[i,:])#*0.001
            dist_u_s2[i] = linalg.norm(all_u[i,:] - all_s2[i,:])#*0.001
        for l in range(L):
            current_time = datetime.utcnow()
            print("\nSimulation {}/{}, Elapsed time: {:.2f} minutes".format(l+1, L,
                                      (current_time - start).total_seconds()/60))
            e, GNSS_noise = self.generate_noise(std_RD, std_AOA, std_GNSS, len(time))
            peak1 = 0
            peak2 = 0
            peak3 = 0
            for i in range(len(time)):
                noisy_s1, noisy_s2 = self.loc.add_GNSS_error(all_s1[i,:], all_s2[i,:], GNSS_noise[:,i])
                k_w_GNSS_error = self.loc.get_real_vector(all_u[i,:], noisy_s1, noisy_s2)
                u_hat[i,:] = self.loc.estimate(noisy_s1, noisy_s2, k_w_GNSS_error, e[:,i], Q)
                est_dist_u_s1[i] = linalg.norm(u_hat[i] - noisy_s1)
                #est_dist_u_s2[i] = linalg.norm(u_hat[i] - noisy_s2)
                meas_error[i] = linalg.norm(u_hat[i,:] - all_u[i,:])
            dist_diff_between_meas = abs(est_dist_u_s1[1:] - est_dist_u_s1[:-1])
            peak1 = dist_diff_between_meas.argmax()
            if (dist_diff_between_meas[:peak1-840].size > 0):
                peak2 = dist_diff_between_meas[:peak1-840].argmax()
            if (dist_diff_between_meas[peak1+600:].size > 0):
                peak3 = dist_diff_between_meas[peak1+600:].argmax() + peak1 + 600
            if (dist_diff_between_meas[peak3] > dist_diff_between_meas[peak2]):
                peak2 = peak3
            if (peak2 < peak1):
                aux = peak1
                peak1 = peak2
                peak2 = aux
            print("Peak 1: {:.4f}, Peak 2 : {:.4f}, Peak3: {:.4f}".format(peak1/60, peak2/60, peak3/60))
            print("Peak 1: {:.2f}, Peak 2 : {:.2f}, Peak3: {:.2f}".format(dist_diff_between_meas[peak1], dist_diff_between_meas[peak2], dist_diff_between_meas[peak3]))
            #print(meas_error.argmax()/60)
            mean_est_dist_u_s1 += est_dist_u_s1
            #mean_est_dist_u_s2 += est_dist_u_s2
            mean_meas_error += meas_error
            meas_and_model_error = self.model_sim(sat_u, sat_s1, sat_s2, u_hat, peak1-420, peak1+300, initial_date, GNSS_noise, e, meas_error)
            meas_and_model_error = self.model_sim(sat_u, sat_s1, sat_s2, u_hat, peak2-420, peak2+300, initial_date, GNSS_noise, e, meas_and_model_error)
            mean_meas_and_model_error += meas_and_model_error
            all_u_hat += u_hat
        mean_est_dist_u_s1[:] = mean_est_dist_u_s1/L
        #mean_est_dist_u_s2[:] = mean_est_dist_u_s2/L
        mean_meas_error[:] = mean_meas_error/L
        mean_meas_and_model_error[:] = mean_meas_and_model_error/L
        all_u_hat[:] = all_u_hat/L
        with open(output, 'w') as data_file:
            outRCRB = csv.writer(data_file)
            outRCRB.writerow(["Time", "Mean measurement error", "Mean measurement and model error", "Mean estimated distance u-s1", "Distance u-s1", "All u_x", "All u_y", "All u_z", "All u_hat_x", "All u_hat_y", "All u_hat_z"])
            for i in range(len(time)):
                outRCRB.writerow([time[i], mean_meas_error[i], mean_meas_and_model_error[i],
                                  mean_est_dist_u_s1[i], dist_u_s1[i],
                                  all_u[i,0,0], all_u[i,1,0], all_u[i,2,0],
                                  all_u_hat[i,0,0], all_u_hat[i,1,0], all_u_hat[i,2,0]])
        end = datetime.utcnow()
        print("Duration: {}".format(end-start))
        self.plot_model(output)

    def model_sim(self, sat_u, sat_s1, sat_s2, u_hat, peak_t0, peak_tf, initial_date, GNSS_noise, e, meas_error):
        std_RD = 10.0
        std_AOA = 1.0*pi/180.0
        std_GNSS = 10.0
        t0 = float(peak_t0)
        studied_date = initial_date + timedelta(seconds=t0)
        sat_s1.updateOrbitalParameters(studied_date)
        sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
        s1 = sat_s1.getXYZ().copy()
        s2 = sat_s2.getXYZ().copy()
        noisy_s1, noisy_s2 = self.loc.add_GNSS_error(s1, s2, GNSS_noise[:,int(t0)])
        u_vel = (u_hat[peak_t0,:] - u_hat[peak_t0-1,:])/1
        model_dt = 0.01
        seconds = int(peak_tf - t0)
        N = int(seconds/model_dt)
        x = zeros([6,1])
        iterations = 1
        h = model_dt/iterations
        self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, noisy_s1)
        self.mainSat_v = self.Q_poQ_ip*sat_s1.getInertialVel().copy()
        lvlh_u = self.transformPosition(u_hat[peak_t0,:])
        lvlh_v = self.transformVelocity(u_vel)
        x[0:3] = lvlh_u.copy()
        x[3:6] = lvlh_v.copy()
        sat_s1_n = sat_s1.n
        rc = sat_s1.a*sqrt(1 - sat_s1.e**2)
        meas_and_model_error = meas_error.copy()
        studied_date = initial_date + timedelta(seconds=t0)
        one_div_model_dt = int(1/model_dt)
        for s in range(seconds):
            e, GNSS_noise = self.generate_noise(std_RD, std_AOA, std_GNSS, one_div_model_dt)
            for i in range(one_div_model_dt):#N):
                studied_date += timedelta(seconds=model_dt)
                sat_u.updateOrbitalParameters(studied_date)
                sat_s1.updateOrbitalParameters(studied_date)
                u = sat_u.getXYZ().copy()
                s1 = sat_s1.getXYZ().copy()
                noisy_s1, noisy_s2 = self.loc.add_GNSS_error(s1, s2, GNSS_noise[:,i])
                self.initTransformationData(sat_s1.incl, sat_s1.RAAN, sat_s1.theta, sat_s1.w, noisy_s1)
                model_pos = self.ekf.RK4(x, h, iterations, sat_s1_n, rc)
                x = model_pos.copy()
            new_u_hat = self.orbital2inertial(x[:3])
            meas_and_model_error[s+int(t0)+1] = linalg.norm(new_u_hat - u)
            u_hat[s+int(t0)+1,:] = new_u_hat
        return meas_and_model_error

    def plot_model(self, filename="model.csv"):
        time = arange(0, 6000)#*2
        dist_u_s1 = zeros(len(time))
        dist_u_s2 = zeros(len(time))
        est_dist_u_s1 = zeros(len(time))
        #est_dist_u_s2 = zeros(len(time))
        dist_diff_between_meas = zeros(len(time)-1)
        mean_est_dist_u_s1 = zeros(len(time))
        #mean_est_dist_u_s2 = zeros(len(time))
        meas_error = zeros(len(time))
        mean_meas_error = zeros(len(time))
        mean_meas_and_model_error = zeros(len(time))
        all_u_hat = zeros([len(time), 3, 1])
        all_u = zeros([len(time), 3, 1])
        all_s1 = zeros([len(time), 3, 1])
        all_s2 = zeros([len(time), 3, 1])
        with open(filename, newline='') as f:
            reader = csv.reader(f, delimiter=",")
            next(reader)
            data = array(list(reader), dtype=object)
        for i in range(len(data)):
            time[i] = float(data[i,0])
            mean_meas_error[i] = float(data[i,1])
            mean_meas_and_model_error[i] = float(data[i,2])
            mean_est_dist_u_s1[i] = float(data[i,3])
            dist_u_s1[i] = float(data[i,4])
            all_u[i,0,0] = float(data[i,5])
            all_u[i,1,0] = float(data[i,6])
            all_u[i,2,0] = float(data[i,7])
            all_u_hat[i,0,0] = float(data[i,8])
            all_u_hat[i,1,0] = float(data[i,9])
            all_u_hat[i,2,0] = float(data[i,10])
        time = time/60
        fig, ax = subplots(1, 1)
        ax.semilogy(time, mean_meas_error, '-', linewidth=2.0,
                    markersize=12, label="Measurement error")
        ax.semilogy(time, mean_meas_and_model_error, '-', linewidth=2.0,
                    markersize=12, label="Measurement + Model error")
        ax.set(xlabel="Time [min]",
               ylabel="Error {} [m]".format(r'$||\mathbf{u} - \mathbf{\hat{u}}||$'))
        ax.grid()
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        #ax.legend(fontsize=14)
        ax.set_ylim(10, 10**10)
        ax.tick_params(which="both", direction="in", labelsize=14,
                       bottom=True, top=True, left=True, right=False)
        tight_layout()
        savefig("vidal12", transparent=True, bbox_inches="tight", dpi=200, pad_inches=0.01)
        fig2, ax2 = subplots(1, 1)
        ax2.semilogy(time, mean_est_dist_u_s1, '-', linewidth=2.0,
                     markersize=12, label="Estimated distance u-s1")
        ax2.semilogy(time, dist_u_s1, '-', linewidth=2.0,
                     markersize=12, label="Real distance u-s1")
        ax2.set(xlabel="Time [min]",
                ylabel="Distance {} [m]".format(r'$||\mathbf{u} - \mathbf{s_1}||$'))
        ax2.grid()
        ax2.xaxis.label.set_size(16)
        ax2.yaxis.label.set_size(16)
        #ax2.legend(fontsize=14)
        ax2.tick_params(which="both", direction="in", labelsize=14,
                        bottom=True, top=True, left=True, right=False)
        tight_layout()
        savefig("vidal10_2", transparent=True, bbox_inches="tight", dpi=200, pad_inches=0.01)
        ax3 = figure().add_subplot(projection='3d')
        ax3.plot(all_u[:,0,0], all_u[:,1,0], all_u[:,2,0], '-', linewidth=2.0,
                 markersize=12, label="Real orbit")
        ax3.plot(all_u_hat[:,0,0], all_u_hat[:,1,0], all_u_hat[:,2,0], '-', linewidth=2.0,
                 markersize=12, label="Estimated orbit")
        ax3.grid()
        ax3.xaxis.label.set_size(16)
        ax3.yaxis.label.set_size(16)
        #ax3.legend(fontsize=14)
        ax3.tick_params(which="both", direction="in", labelsize=14,
                        bottom=True, top=True, left=True, right=False)
        tight_layout()
        savefig("vidal13", transparent=True, bbox_inches="tight", dpi=200, pad_inches=0.18)
        fig4, ax4 = subplots(1, 1)
        ax4.semilogy(time[1:], abs(abs(mean_est_dist_u_s1[:-1]) - abs(mean_est_dist_u_s1[1:])), '-',
                     linewidth=2.0, markersize=12, label="Distance difference between measurements")
        ax4.set(xlabel="Time [min]",
                ylabel="Distance difference [m]")
        ax4.grid()
        ax4.xaxis.label.set_size(16)
        ax4.yaxis.label.set_size(16)
        #ax4.legend(fontsize=14)
        ax4.tick_params(which="both", direction="in", labelsize=14,
                        bottom=True, top=True, left=True, right=False)
        tight_layout()
        savefig("vidal11", transparent=True, bbox_inches="tight", dpi=200, pad_inches=0.01)
        show()

    def plot_orbit(self, l0):
        std_RD = 10.0
        std_AOA = 1.0*pi/180.0
        std_GNSS = 10.0
        Q = self.loc.get_Q(std_RD, std_AOA)
        data_path = resource_filename("pypredict","data/")
        sat_u = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s1 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2.setTLE("1 44814U 19081L   20321.73053029  .00001305  00000-0  63025-4 0  9996",
                      "2 44814  97.4788  21.6285 0013387  80.2501 280.0246 15.20374749 54001")
        with open("u_TLEs.txt", 'r') as u_TLEs:
            sat_u_TLEs = u_TLEs.readlines()

        with open("s1_TLEs.txt", 'r') as s1_TLEs:
            sat_s1_TLEs = s1_TLEs.readlines()

        date0 = datetime(2020, 11, 17, 00, 13, 33, 0)
        dep_date = date0 + timedelta(minutes=10)
        initial_date = dep_date + timedelta(days=3)
        time = arange(0, 1800)
        real_u = zeros([len(time), 3])
        est_u = zeros([len(time), 3])
        sat_u.setTLE(sat_u_TLEs[1 + l0*3], sat_u_TLEs[2 + l0*3])
        sat_s1.setTLE(sat_s1_TLEs[1 + l0*3], sat_s1_TLEs[2 + l0*3])
        e, GNSS_noise = self.generate_noise(std_RD, std_AOA, std_GNSS, len(time))
        for i in range(len(time)):
            studied_date = initial_date + timedelta(seconds=float(time[i]))
            sat_u.updateOrbitalParameters(studied_date)
            sat_s1.updateOrbitalParameters(studied_date)
            sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
            u = sat_u.getXYZ().copy()
            real_u[i,:] = u.T
            s1 = sat_s1.getXYZ().copy()
            s2 = sat_s2.getXYZ().copy()
            noisy_s1, noisy_s2 = self.loc.add_GNSS_error(s1, s2, GNSS_noise[:,i])
            k_w_GNSS_error = self.loc.get_real_vector(u, noisy_s1, noisy_s2)
            u_hat = self.loc.estimate(noisy_s1, noisy_s2, k_w_GNSS_error, e[:,i], Q)
            est_u[i,:] = u_hat.T
        ax = figure().add_subplot(projection='3d')
        ax.plot(est_u[:,0], est_u[:,1], est_u[:,2], '-', linewidth=2.0, markersize=12,
                 label="Estimated orbit")
        #ax.plot(real_u[:,0], real_u[:,1], real_u[:,2], '-', linewidth=2.0, markersize=12,
        #        label="Real orbit")
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        ax.legend(fontsize=14)
        tight_layout()
        show()


    def plot_distance_over_time(self, L0=[582, 2989, 963, 3714, 425]):
        data_path = resource_filename("pypredict","data/")
        sat_u = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s1 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2 = Sat(name="FLOCK 4P-1", tlepath="{}planet.txt".format(data_path), cat="Planet Labs")
        sat_s2.setTLE("1 44814U 19081L   20321.73053029  .00001305  00000-0  63025-4 0  9996",
                      "2 44814  97.4788  21.6285 0013387  80.2501 280.0246 15.20374749 54001")
        with open("u_TLEs.txt", 'r') as u_TLEs:
            sat_u_TLEs = u_TLEs.readlines()

        with open("s1_TLEs.txt", 'r') as s1_TLEs:
            sat_s1_TLEs = s1_TLEs.readlines()

        date0 = datetime(2020, 11, 17, 00, 13, 33, 0)
        dep_date = date0 + timedelta(minutes=10)
        initial_date = dep_date + timedelta(days=3)
        time = arange(0, 18000)
        altitude = zeros(len(time))
        dist_u_s1 = zeros(len(time))
        dist_u_s2 = zeros(len(time))
        for l0 in L0:
            sat_u.setTLE(sat_u_TLEs[1 + l0*3], sat_u_TLEs[2 + l0*3])
            sat_s1.setTLE(sat_s1_TLEs[1 + l0*3], sat_s1_TLEs[2 + l0*3])
            for i in range(len(time)):
                studied_date = initial_date + timedelta(days=float(time[i]))
                sat_u.updateOrbitalParameters(studied_date)
                sat_s1.updateOrbitalParameters(studied_date)
                sat_s2.updateOrbitalParameters(studied_date + timedelta(seconds=4))
                u = sat_u.getXYZ().copy()
                s1 = sat_s1.getXYZ().copy()
                s2 = sat_s2.getXYZ().copy()
                dist_u_s1[i] += linalg.norm(u - s1)*0.001
                dist_u_s2[i] += linalg.norm(u - s2)*0.001
                altitude[i] += sat_u.getAlt()*0.001

        dist_u_s1 = dist_u_s1/len(L0)
        dist_u_s2 = dist_u_s2/len(L0)
        altitude = altitude/len(L0)
        fig, ax = subplots(1, 1)
        ax.plot(time, altitude)
        ax.plot(time, dist_u_s1)
        ax.plot(time, dist_u_s2)
        ax.xaxis.label.set_size(16)
        ax.yaxis.label.set_size(16)
        tight_layout()
        show()

