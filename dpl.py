from calcOrbitParam import Calc
from pyorbital import tlefile
from sat import Sat
from numpy import asscalar, cos, matrix, pi, sin

class Dpl(object):
    __slots__ = ["calc", "Q_pi", "Q_po", "Q_op", "v", "r"]
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
        deployer.updateWithDragEffect()
        aux = deployer.getXYZ()
        self.r = [aux[0,0], aux[1,0], aux[2,0]]
        v_p, v_q = deployer.getPerifocalVel()
        v_peri = matrix([[v_p], [v_q], [0]])
        v_orb = self.perifocal2orbital(deployer.theta, v_peri)
        v_orb[0][0] = asscalar(v_orb[0][0]) + vel[0]
        v_orb[1][0] = asscalar(v_orb[1][0]) + vel[1]
        v_orb[2][0] = asscalar(v_orb[2][0]) + vel[2]
        v_peri = self.orbital2perifocal(deployer.theta, v_orb)
        v_x, v_y, v_z = self.perifocal2inertial(deployer, v_peri)
        self.v = [asscalar(v_x), asscalar(v_y), asscalar(v_z)]

    def updateSat(self, sat):
        self.calc.newCalc(r1=self.r, v=self.v)
        sat.updateEpoch()
        sat.updateGST0()
        sat.setInclination(self.calc.i)
        sat.setRAAN0(self.calc.RAAN)
        sat.setArgPerigee0(self.calc.w)
        sat.setEccentricity(self.calc.e_scalar)
        sat.setMeanAnomaly0(self.calc.MA)
        sat.setTrueAnomaly(self.calc.theta)
        sat.setSemiMajorAxis(self.calc.a)
        sat.setSemilatusRectum(self.calc.a*(1 - self.calc.e_scalar**2))
        sat.setMeanVelocity(self.calc.n)
        sat.p = self.calc.a*(1 - self.calc.e_scalar**2)
        sat.updateWithDragEffect()

    def deploy(self, category, deployer, dplyr_mass, dplyd_mass, name, vel):
        self.calcPosAndVel(deployer, vel)
        newSat = Sat(name=name,
                tle=tlefile.read(deployer.name, "TLE/cubesat.txt"),
                cat=category)
        self.updateSat(newSat)
        B = (2*(0.034*0.084 + 0.034*0.028 + 0.084*0.028)/dplyd_mass)/3
        newSat.setBallisticCoeff(B)
        dplyr_mass = dplyr_mass - dplyd_mass
        dplyr_vel = [-vel[0]*dplyd_mass/dplyr_mass,
                     -vel[1]*dplyd_mass/dplyr_mass,
                     -vel[2]*dplyd_mass/dplyr_mass]
        self.calcPosAndVel(deployer, dplyr_vel)
        self.updateSat(deployer)
        deployer.setBallisticCoeff((2*0.1*(2*0.3 + 0.1)/dplyr_mass)/3)
        return newSat
