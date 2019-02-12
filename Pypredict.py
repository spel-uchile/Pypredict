from pyorbital import tlefile
from calcOrbitParam import Calc
from sat import Sat
from tdoa import TDOA
from gui import GUI

SUCHAI = Sat(name="SUCHAI", tle=tlefile.read("SUCHAI", "TLE/cubesat.txt"), cat="CubeSat")
HODOYOSHI3 = Sat(name="HODOYOSHI-3", tle=tlefile.read("HODOYOSHI-3", "TLE/resource.txt"), cat="Earth Resources")
HODOYOSHI4 = Sat(name="HODOYOSHI-4", tle=tlefile.read("HODOYOSHI-4", "TLE/resource.txt"), cat="Earth Resources")
#tdoa = TDOA()
[x, y, z] = SUCHAI.getXYZ()
r = [x, y, z]
[v_x, v_y, v_z] = SUCHAI.getVelocityVector()
v = [v_x, v_y, v_z]
calc = Calc(r, v=v)
FE1 = Sat(name="FE1", tle=tlefile.read("SUCHAI", "TLE/cubesat.txt"), cat="Femto-satellite")
FE1.updateEpoch()
FE1.updateGST0()
FE1.setInclination(calc.i)
FE1.setRAAN0(calc.RAAN)
FE1.setArgPerigee0(calc.w)
FE1.setEccentricity(calc.e_scalar)
FE1.setMeanAnomaly0(calc.MA)
FE1.setTrueAnomaly(calc.theta)
FE1.setSemiMajorAxis(calc.a)
FE1.setSemilatusRectum(calc.a*(1 - calc.e_scalar**2))
FE1.setMeanVelocity(calc.n)
FE1.p = calc.a*(1 - calc.e_scalar**2)
FE1.updateOrbitalParameters()

Sats = [SUCHAI, HODOYOSHI3, HODOYOSHI4, FE1]
#x, y, z = tdoa.calculateLocation(SUCHAI, IRIDIUM83, IRIDIUM90, IRIDIUM91, IRIDIUM95)
#print("NOI real x: " + str(IRIDIUM83.x) + "    NOI calculated x: " + str(x))
#print("NOI real y: " + str(IRIDIUM83.y) + "    NOI calcualted y: " + str(y))
GUI = GUI(Sats=Sats)
