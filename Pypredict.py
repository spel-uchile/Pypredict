from pyorbital import tlefile
from calcOrbitParam import Calc
from sat import Sat
from tdoa import TDOA
from gui import GUI

#tlefile.TLE_URLS = ("https://celestrak.com/NORAD/elements/cubesat.txt", )
#tlefile.fetch("TLE/cubesat.txt")
#tlefile.TLE_URLS = ('https://celestrak.com/NORAD/elements/iridium.txt', )
#tlefile.fetch("TLE/iridium.txt")
tlefile.TLE_URLS =  ("https://celestrak.com/NORAD/elements/argos.txt",
                     "https://celestrak.com/NORAD/elements/cubesat.txt",
                     "https://celestrak.com/NORAD/elements/dmc.txt",
                     "https://celestrak.com/NORAD/elements/dmc.txt",
                     "https://celestrak.com/NORAD/elements/intelsat.txt",
                     "https://celestrak.com/NORAD/elements/iridium.txt",
                     "https://celestrak.com/NORAD/elements/iridium-NEXT.txt",
                     "https://celestrak.com/NORAD/elements/noaa.txt",
                     "https://celestrak.com/NORAD/elements/planet.txt",
                     "https://celestrak.com/NORAD/elements/resource.txt",
                     "https://celestrak.com/NORAD/elements/sarsat.txt",
                     "https://celestrak.com/NORAD/elements/spire.txt",
                     "https://celestrak.com/NORAD/elements/tdrss.txt",
                     "https://celestrak.com/NORAD/elements/visual.txt",
                     "https://celestrak.com/NORAD/elements/weather.txt", )
SUCHAI = Sat(name="SUCHAI", tle=tlefile.read("SUCHAI", "TLE/cubesat.txt"), cat="CubeSat")
IRIDIUM90 = Sat(name="IRIDIUM 90", tle=tlefile.read("IRIDIUM 90 [-]", "TLE/iridium.txt"), cat="Iridium")
IRIDIUM91 = Sat(name="IRIDIUM 91", tle=tlefile.read("IRIDIUM 91 [+]", "TLE/iridium.txt"), cat="Iridium")
IRIDIUM95 = Sat(name="IRIDIUM 95", tle=tlefile.read("IRIDIUM 95 [+]", "TLE/iridium.txt"), cat="Iridium")
#tdoa = TDOA()
Sats = [SUCHAI, IRIDIUM90, IRIDIUM91, IRIDIUM95]
#x, y, z = tdoa.calculateLocation(SUCHAI, IRIDIUM83, IRIDIUM90, IRIDIUM91, IRIDIUM95)
#print("NOI real x: " + str(IRIDIUM83.x) + "    NOI calculated x: " + str(x))
#print("NOI real y: " + str(IRIDIUM83.y) + "    NOI calcualted y: " + str(y))
GUI = GUI(Sats=Sats)
