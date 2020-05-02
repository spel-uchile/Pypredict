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
from calcOrbitParam import Calc
from pyorbital import tlefile
from sat import Sat
from numpy import asscalar, cos, matrix, pi, sin, sqrt

class Dpl(object):
    __slots__ = ["calc", "Q_pi", "Q_po",
                 "Q_op", "v", "r"]
    def __init__(self):
        self.calc = Calc()

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

    def calcPosAndVel(self, deployer, vel):
        aux = deployer.getXYZ()
        self.r = [aux[0,0], aux[1,0], aux[2,0]]
        v_p, v_q = deployer.getPerifocalVel()
        v_peri = matrix([[v_p], [v_q], [0]])
        v_orb = self.perifocal2orbital(deployer.theta, v_peri)
        v_orb[0,0] = v_orb[0,0] + vel[0]
        v_orb[1,0] = v_orb[1,0] + vel[1]
        v_orb[2,0] = v_orb[2,0] + vel[2]
        v_peri = self.orbital2perifocal(deployer.theta, v_orb)
        v = self.perifocal2inertial(deployer, v_peri)
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

    def deploy(self, category, deployer, dplyr_mass, dplyd_mass, name, vel, date=None):
        self.calcPosAndVel(deployer, vel)
        deployer_name, line1, line2 = deployer.createTLE(date)
        tle = tlefile.read(deployer_name, line1=line1, line2=line2)
        newSat = Sat(name=name, tle=tle, cat=category)
        self.updateSat(newSat, date)
        B = (2*(0.034*0.084 + 0.034*0.028 + 0.084*0.028))/6/dplyd_mass
        newSat.setBallisticCoeff(B)
        newSat.createTLE(date)
        dplyr_mass = dplyr_mass - dplyd_mass
        dplyr_vel = [-vel[0]*dplyd_mass/dplyr_mass,
                     -vel[1]*dplyd_mass/dplyr_mass,
                     -vel[2]*dplyd_mass/dplyr_mass]
        self.calcPosAndVel(deployer, dplyr_vel)
        self.updateSat(deployer, date)
        B = (0.1*0.1*2 + 4*0.3*0.1)/6/dplyr_mass
        deployer.setBallisticCoeff(B)
        deployer.createTLE(date)
        return newSat
