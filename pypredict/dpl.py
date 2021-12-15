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
from pypredict.calcOrbitParam import Calc
from pypredict.sat import Sat
from numpy import arange, array, asscalar, concatenate, cos, flip, matrix, pi, sin, sqrt
from sgp4.earth_gravity import wgs72, wgs84
from sgp4.io import twoline2rv

class Dpl(object):
    __slots__ = ["calc", "Q_pi", "Q_po",
                 "Q_op", "v", "r", "vel_weight"]
    def __init__(self):
        self.calc = Calc()
        self.vel_weight = 1000

    def __call__(self):
        return self

    def perifocal2inertial(self, sat, v):
        sin_RAAN = sin(sat.RAAN)
        cos_RAAN = cos(sat.RAAN)
        sin_i = sin(sat.incl)
        cos_i = cos(sat.incl)
        sin_w = sin(sat.w)
        cos_w = cos(sat.w)
        self.Q_pi = matrix([[-sin_RAAN*cos_i*sin_w + cos_RAAN*cos_w,
                             -sin_RAAN*cos_i*cos_w - cos_RAAN*sin_w,
                             sin_RAAN*sin_i],
                            [cos_RAAN*cos_i*sin_w + sin_RAAN*cos_w,
                             cos_RAAN*cos_i*cos_w - sin_RAAN*sin_w,
                             -cos_RAAN*sin_i],
                            [sin_i*sin_w, sin_i*cos_w, cos_i]])
        return self.Q_pi*v

    def perifocal2orbital(self, theta, v):
        phi = theta + pi/2
        self.Q_po = matrix([[cos(phi), sin(phi), 0],
                            [0, 0, -1],
                            [-sin(phi), cos(phi), 0]])
        return self.Q_po*v

    def orbital2perifocal(self, theta, v):
        phi = theta + pi/2
        self.Q_op = matrix([[cos(phi), 0, -sin(phi)],
                            [sin(phi), 0, cos(phi)],
                            [0, -1, 0]])
        return self.Q_op*v

    def calcPosAndVel(self, sat, vel):
        aux = sat.getXYZ()
        self.r = [aux[0,0], aux[1,0], aux[2,0]]
        v_p, v_q = sat.getPerifocalVel()
        v_peri = matrix([[v_p], [v_q], [0]])
        v_orb = self.perifocal2orbital(sat.theta, v_peri)
        v_orb[0,0] = v_orb[0,0] + vel[0]
        v_orb[1,0] = v_orb[1,0] + vel[1]
        v_orb[2,0] = v_orb[2,0] + vel[2]
        v_peri = self.orbital2perifocal(sat.theta, v_orb)
        v = self.perifocal2inertial(sat, v_peri)
        self.v = [v[0,0], v[1,0], v[2,0]]

    def updateSat(self, sat, date):
        self.calc.newCalc(r1=self.r, v=self.v)
        sat.setInclination(self.calc.i)
        sat.setRAAN0(self.calc.RAAN)
        sat.setRAAN(self.calc.RAAN)
        sat.setArgPerigee0(self.calc.w)
        sat.setArgPerigee(self.calc.w)
        sat.setEccentricity(self.calc.e_scalar)
        sat.setMeanAnomaly0(self.calc.MA)
        sat.setMeanAnomaly(self.calc.MA)
        sat.setTrueAnomaly(self.calc.theta)
        sat.setSemiMajorAxis(self.calc.a)
        sat.setMeanMotion(self.calc.n)
        sat.setSemilatusRectum(self.calc.a*(1 - self.calc.e_scalar**2))
        sat.setSpecAngMomentum(sqrt(sat.a*sat.mu*(1 - self.calc.e_scalar**2)))
        sat.updateEpoch(date=date)

    def fitTLE(self, sat, line1, line2, pos0, vel0, date, step=0.0001, den=2.0, best_err=None):
        """
        This methods searches for the second line of the TLE that better fits
        the position pos0 and velocity vel0. This is done iteratively by
        searching in the neighbourhood of the orbital elements of the line
        line2. Since iterating through six orbital elements with six nested
        for loops takes a huge amount of time, we took a different approach.
        We selected different groups of orbital elements and changed only
        them leaving the rest constant, saving the line that yields the lowest
        error in each step. The search is not linear to make it faster. The
        steps are denser in the vecinity of the current best line, but they
        can also allow a big step if we are not near the solution.
        """
        if (best_err is None):
            best_err = 1000000.0
        pos0 = [x*0.001 for x in pos0]
        vel0 = [x*0.001 for x in vel0]
        year = date.year
        month = date.month
        day = date.day
        hour = date.hour
        minute = date.minute
        sec = date.second + date.microsecond*0.000001
        str_e0 = line2[25:34]
        base_i = array([-21, -4, -1, 0, 1, 4, 21])*0.0001
        base_RAAN = array([-50, -15, -5, -1, 0, 1, 5, 15, 50])*0.0001
        base = array([-256, -128, -64, -32, -16, -8, -2,
                       0, 2, 8, 16, 32, 64, 128, 256])*step
        base_MM = array([-6000, -4500, -3000, -1000, -500, -100, -50, -10, -2,
                          0, 2, 10, 50, 100, 500, 1000, 3000, 4500, 6000])*step*0.0001
        all_i = base_i + float(line2[8:16])
        all_RAAN = base_RAAN + float(line2[17:25])
        all_e = base*0.001 + float(line2[26:33])*0.0000001
        all_MA = base + float(line2[43:52])
        all_w = base + float(line2[34:42])
        all_MM = base_MM + float(line2[52:63])
        if (all_i[0] < 0):
            all_i += abs(all_i[0])
        if (all_RAAN[0] < 0):
            all_RAAN += abs(all_RAAN[0])
        if (all_e[0] < 0):
            all_e += abs(all_e[0]) 
        if (all_w[0] < 0):
            all_w += abs(all_w[0])
        if (all_MA[0] < 0):
            all_MA += abs(all_MA[0])
        if (all_MM[0] < 0):
            all_MM += abs(all_MM[0])
        best_line2 = line2
        beginning = line2[0:17]
        ending = line2[51:]
        for RAAN in all_RAAN:
            for w in all_w:
                for MA in all_MA:
                    new_line2 = "{}{:8.4f}{}{:8.4f}{:9.4f}{}".format(beginning, RAAN, str_e0,
                                                                     w, MA, ending)
                    pos_err, vel_err = self.get_error(sat, line1, new_line2, year, month,
                                                      day, hour, minute, sec, pos0, vel0)
                    self.vel_weight = (vel_err*150000.0 - pos_err)*1000.0
                    self.vel_weight = (self.vel_weight < 20.0)*20.0 + (self.vel_weight >= 20.0)*self.vel_weight
                    self.vel_weight = (self.vel_weight > 2500.0)*2500.0 + (self.vel_weight <= 2500.0)*self.vel_weight
                    if (pos_err + self.vel_weight*vel_err < best_err):
                        best_line2 = new_line2
                        best_pos_err = pos_err
                        best_vel_err = vel_err
                        best_err = best_pos_err + self.vel_weight*best_vel_err
                        #print("Step: {:6.4f},  Pos. err: {:9.4f} m,  Vel. err: {:9.6f} m/s,  Weight: {:7.2f}".format(step, best_pos_err*1000, best_vel_err*1000, self.vel_weight))
                        #print("{}\n".format(best_line2))
        str_RAAN = " {:8.4f} ".format(float(best_line2[17:25]))
        beginning = best_line2[0:8]
        middle = best_line2[33:52]
        ending = best_line2[63:]
        for e in all_e:
            str_e = "{}".format("{:.8f}".format(e)[2:-1])
            for i in all_i:
                for MM in all_MM:
                    new_line2 = "{}{:8.4f}{}{}{}{:11.8f}{}".format(beginning, i, str_RAAN,
                                                                   str_e, middle, MM, ending)
                    pos_err, vel_err = self.get_error(sat, line1, new_line2, year, month,
                                                      day, hour, minute, sec, pos0, vel0)
                    self.vel_weight = (vel_err*150000.0 - pos_err)*1000.0
                    self.vel_weight = (self.vel_weight < 20.0)*20.0 + (self.vel_weight >= 20.0)*self.vel_weight
                    self.vel_weight = (self.vel_weight > 2500.0)*2500.0 + (self.vel_weight <= 2500.0)*self.vel_weight
                    if (pos_err + self.vel_weight*vel_err < best_err):
                        best_line2 = new_line2
                        best_pos_err = pos_err
                        best_vel_err = vel_err
                        best_err = best_pos_err + self.vel_weight*best_vel_err
                        #print("Step: {:6.4f},  Pos. err: {:9.4f} m,  Vel. err: {:9.6f} m/s,  Weight: {:7.2f}".format(step, best_pos_err*1000, best_vel_err*1000, self.vel_weight))
                        #print("{}\n".format(best_line2))
        all_e = base*0.001 + float(best_line2[26:33])*0.0000001
        all_w = base + float(best_line2[34:42])
        all_MA = base + float(best_line2[43:52])
        if (all_e[0] < 0):
            all_e += abs(all_e[0])
        if (all_w[0] < 0):
            all_w += abs(all_w[0])
        if (all_MA[0] < 0):
            all_MA += abs(all_MA[0])
        beginning = best_line2[0:26]
        ending = best_line2[51:]
        for e in all_e:
            str_e = "{} ".format("{:.8f}".format(e)[2:-1])
            for w in all_w:
                for MA in all_MA:
                    new_line2 = "{}{}{:8.4f}{:9.4f}{}".format(beginning, str_e, w, MA, ending)
                    pos_err, vel_err = self.get_error(sat, line1, new_line2, year, month,
                                                      day, hour, minute, sec, pos0, vel0)
                    self.vel_weight = (vel_err*150000.0 - pos_err)*1000.0
                    self.vel_weight = (self.vel_weight < 20.0)*20.0 + (self.vel_weight >= 20.0)*self.vel_weight
                    self.vel_weight = (self.vel_weight > 2500.0)*2500.0 + (self.vel_weight <= 2500.0)*self.vel_weight
                    if (pos_err + self.vel_weight*vel_err < best_err):
                        best_line2 = new_line2
                        best_pos_err = pos_err
                        best_vel_err = vel_err
                        best_err = best_pos_err + self.vel_weight*best_vel_err
                        #print("Step: {:6.4f},  Pos. err: {:9.4f} m,  Vel. err: {:9.6f} m/s,  Weight: {:7.2f}".format(step, best_pos_err*1000, best_vel_err*1000, self.vel_weight))
                        #print("{}\n".format(best_line2))
        all_w = base + float(best_line2[34:42])
        all_MA = base + float(best_line2[43:52])
        all_MM = base_MM + float(best_line2[52:63])
        if (all_w[0] < 0):
            all_w += abs(all_w[0])
        if (all_MA[0] < 0):
            all_MA += abs(all_MA[0])
        if (all_MM[0] < 0):
            all_MM += abs(all_MM[0])
        beginning = best_line2[0:34]
        ending = best_line2[63:]
        for w in all_w:
            for MA in all_MA:
                for MM in all_MM:
                    new_line2 = "{}{:8.4f}{:9.4f}{:12.8f}{}".format(beginning, w, MA, MM, ending)
                    pos_err, vel_err = self.get_error(sat, line1, new_line2, year, month,
                                                      day, hour, minute, sec, pos0, vel0)
                    self.vel_weight = (vel_err*150000.0 - pos_err)*1000.0
                    self.vel_weight = (self.vel_weight < 20.0)*20.0 + (self.vel_weight >= 20.0)*self.vel_weight
                    self.vel_weight = (self.vel_weight > 2500.0)*2500.0 + (self.vel_weight <= 2500.0)*self.vel_weight
                    if (pos_err + self.vel_weight*vel_err < best_err):
                        best_line2 = new_line2
                        best_pos_err = pos_err
                        best_vel_err = vel_err
                        best_err = best_pos_err + self.vel_weight*best_vel_err
                        #print("Step: {:6.4f},  Pos. err: {:9.4f} m,  Vel. err: {:9.6f} m/s,  Weight: {:7.2f}".format(step, best_pos_err*1000, best_vel_err*1000, self.vel_weight))
                        #print("{}\n".format(best_line2))
        all_i = base_i + float(best_line2[8:16])
        all_RAAN = base_RAAN + float(best_line2[17:25])
        if (all_i[0] < 0):
            all_i += abs(all_i[0])
        if (all_RAAN[0] < 0):
            all_RAAN += abs(all_RAAN[0])
        beginning = best_line2[0:8]
        ending = best_line2[25:]
        for i in all_i:
            for RAAN in all_RAAN:
                new_line2 = "{}{:8.4f}{:9.4f}{}".format(beginning, i, RAAN, ending)
                pos_err, vel_err = self.get_error(sat, line1, new_line2, year, month,
                                                  day, hour, minute, sec, pos0, vel0)
                self.vel_weight = (vel_err*150000.0 - pos_err)*1000.0
                self.vel_weight = (self.vel_weight < 20.0)*20.0 + (self.vel_weight >= 20.0)*self.vel_weight
                self.vel_weight = (self.vel_weight > 2500.0)*2500.0 + (self.vel_weight <= 2500.0)*self.vel_weight
                if (pos_err + self.vel_weight*vel_err < best_err):
                    best_line2 = new_line2
                    best_pos_err = pos_err
                    best_vel_err = vel_err
                    best_err = best_pos_err + self.vel_weight*best_vel_err
                    #print("Step: {:6.4f},  Pos. err: {:9.4f} m,  Vel. err: {:9.6f} m/s,  Weight: {:7.2f}".format(step, best_pos_err*1000, best_vel_err*1000, self.vel_weight))
                    #print("{}\n".format(best_line2))
        return best_line2, best_err

    def get_error(self, sat, line1, line2, year, month, day, hour, minute, second, pos0, vel0):
        sat.sat_model = twoline2rv(line1, line2, wgs84)#wgs72)
        pos, vel = sat.sat_model.propagate(year, month, day, hour, minute, second)
        pos_err = sqrt((pos0[0] - pos[0])**2 + (pos0[1] - pos[1])**2 + (pos0[2] - pos[2])**2)
        vel_err = sqrt((vel0[0] - vel[0])**2 + (vel0[1] - vel[1])**2 + (vel0[2] - vel[2])**2)
        return pos_err, vel_err#dist_err + self.vel_weight*vel_err

    def deploy(self, cat, dplyr, dplyr_mass, dplyd_mass, name, vel, date=None):
        """
        Simulates a deployment of a satellite from another. Calculates
        the orbital parameters of the new satellite (the deployed
        satellite) and updates the deployer satellite due to the
        change in momentum considering the mass.

        Parameters
        ----------
        cat        : str
                     The new satellite's category.
        dplyr      : pypredict.sat.Sat object
                     The satellite that is deploying the new satellite.
        dplyr_mass : float
                     The mass of the deployer.
        dplyd_mass : float
                     The mass of the deployed satellite.
        name       : str
                     The name of the deployed satellite.
        vel        : list
                     A list with the 3 components of the velocity
                     of deployment from the deployer's perspective.
        date       : datetime.datetime object, optional
                     Defaults to datetime.utcnow().

        Returns
        -------
        newsat     : pypredict.sat.Sat object
                     A satellite object with the properties of the
                     deployed satellite.
        """
        steps = [0.0002, 0.0001, 0.00006, 0.00004, 0.00002]
        newSat = Sat(name=name, line1=dplyr.line1, line2=dplyr.line2, cat=cat)
        newSat.updateOrbitalParameters(date)
        self.calcPosAndVel(newSat, vel)
        line1 = dplyr.line1
        last_line2 = dplyr.line2
        line2, err = self.fitTLE(newSat, line1, last_line2, self.r, self.v, date, step=0.0001, den=1.8)
        while (line2 != last_line2):
            last_line2 = line2
            for s in steps:
                line2, err = self.fitTLE(newSat, line1, line2, self.r, self.v, date, step=s, den=2.0, best_err=err)

        pos0 = [x*0.001 for x in self.r]
        vel0 = [x*0.001 for x in self.v]
        year = date.year
        month = date.month
        day = date.day
        hour = date.hour
        minute = date.minute
        sec = date.second + date.microsecond*0.000001
        pos_err, vel_err = self.get_error(newSat, line1, line2, year, month,
                                          day, hour, minute, sec, pos0, vel0)
        self.vel_weight = (vel_err*150000.0 - pos_err)*1000.0
        self.vel_weight = (self.vel_weight < 20.0)*20.0 + (self.vel_weight >= 20.0)*self.vel_weight
        self.vel_weight = (self.vel_weight > 2500.0)*2500.0 + (self.vel_weight <= 2500.0)*self.vel_weight
        print("Pos. err:{:10.4f} m,  Vel. err:{:10.6f} m/s,  Weight:{:8.2f}".format(pos_err*1000,
                                                                                    vel_err*1000,
                                                                                    self.vel_weight))
        newSat.setTLE(line1, line2)

        dplyr_mass = dplyr_mass - dplyd_mass
        dplyr_vel = [-vel[0]*dplyd_mass/dplyr_mass,
                     -vel[1]*dplyd_mass/dplyr_mass,
                     -vel[2]*dplyd_mass/dplyr_mass]
        self.calcPosAndVel(dplyr, dplyr_vel)
        dv = sqrt(dplyr_vel[0]**2 + dplyr_vel[1]**2 + dplyr_vel[2]**2)

        last_line2 = dplyr.line2
        line2, err = self.fitTLE(dplyr, line1, last_line2, self.r, self.v, date, step=0.0001, den=1.8)
        while (line2 != last_line2):
            last_line2 = line2
            for s in steps[1:]:
                line2, err = self.fitTLE(dplyr, line1, line2, self.r, self.v, date, step=s, den=1.8, best_err=err)

        pos0 = [x*0.001 for x in self.r]
        vel0 = [x*0.001 for x in self.v]
        pos_err, vel_err = self.get_error(dplyr, line1, line2, year, month,
                                          day, hour, minute, sec, pos0, vel0)
        self.vel_weight = (vel_err*150000.0 - pos_err)*1000.0
        self.vel_weight = (self.vel_weight < 20.0)*20.0 + (self.vel_weight >= 20.0)*self.vel_weight
        self.vel_weight = (self.vel_weight > 2500.0)*2500.0 + (self.vel_weight <= 2500.0)*self.vel_weight
        print("Pos. err:{:10.4f} m,  Vel. err:{:10.6f} m/s,  Weight:{:8.2f}".format(pos_err*1000,
                                                                                    vel_err*1000,
                                                                                    self.vel_weight))

        dplyr.setTLE(line1, line2)
        return newSat
