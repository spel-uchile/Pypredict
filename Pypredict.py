from pyorbital import tlefile
from sat import Sat
from tdoa import TDOA
from gui import GUI

SUCHAI = Sat(name="SUCHAI", tle=tlefile.read("SUCHAI", "TLE/cubesat.txt"), cat="CubeSat")
HODOYOSHI3 = Sat(name="HODOYOSHI-3", tle=tlefile.read("HODOYOSHI-3", "TLE/resource.txt"), cat="Earth Resources")
HODOYOSHI4 = Sat(name="HODOYOSHI-4", tle=tlefile.read("HODOYOSHI-4", "TLE/resource.txt"), cat="Earth Resources")
CUBESATXI_IV = Sat(name="CUBESAT XI-IV (CO-57)", tle=tlefile.read("CUBESAT XI-IV (CO-57)", "TLE/cubesat.txt"), cat="CubeSat")
#tdoa = TDOA()

Sats = [SUCHAI, HODOYOSHI3, HODOYOSHI4, CUBESATXI_IV]
#x, y, z = tdoa.calculateLocation(SUCHAI, IRIDIUM83, IRIDIUM90, IRIDIUM91, IRIDIUM95)
#print("NOI real x: " + str(IRIDIUM83.x) + "    NOI calculated x: " + str(x))
#print("NOI real y: " + str(IRIDIUM83.y) + "    NOI calcualted y: " + str(y))
GUI = GUI(Sats=Sats)
