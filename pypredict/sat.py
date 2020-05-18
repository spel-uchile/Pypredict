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
from datetime import datetime, timedelta
from numpy import abs, arccos, arcsin, arctan, arctan2, array, cos, matrix, pi, sin, sqrt, tan
from pypredict.node import Node
from sgp4.earth_gravity import wgs72
from sgp4.ext import rv2coe
from sgp4.io import twoline2rv

class Sat(Node):
    __slots__ = ["cat", "incl", "RAAN0", "RAAN", "e", "w0", "w", "MA0", "MA",
                 "n", "epoch_year", "epoch_day", "theta", "GST0", "a", "mu",
                 "t0", "p",  "x", "y", "z", "Eq_r", "Po_r", "P_r", "Er_Pr2",
                 "J2", "P_w", "r", "h", "v", "tray_lat", "tray_lng", "B",
                 "bstar", "v_peri", "v_iner", "r_iner", "sat_model", "ndot",
                 "id_launch_data", "element_number", "satnumber", "tray_alt",
                 "line1", "line2", "tlepath"]
    def __init__(self, name="", line1=None, line2=None, tlepath=None, cat=""):
        """
        Parameters
        ----------
        name    : str, optional
                  The satellite's name (default is ).
        line1   : str, optional
                  First line of the two line elements.
        line2   : str, optional
                  Second line of the two line elements.
        tlepath : str, optional
                  Path to the satellite's TLE (default is None).
        cat     : str, optional
                  The satellite's category (default is ).
        """
        self.name = name
        self.cat = cat
        self.tray_lat = []
        self.tray_lng = []
        self.tray_alt = []
        if (tlepath is not None):
            self.readTLE(tlepath)
            self.tlepath = tlepath
        else:
            self.line1 = line1
            self.line2 = line2
        self.sat_model = twoline2rv(self.line1, self.line2, wgs72)
        self.incl = self.sat_model.inclo
        self.RAAN0 = self.sat_model.nodeo
        self.e = self.sat_model.ecco
        self.w0 = self.sat_model.argpo
        self.MA0 = self.sat_model.mo
        self.n = self.sat_model.no_kozai/60
        self.epoch_year = self.sat_model.epochyr
        self.epoch_day = self.sat_model.epochdays
        self.bstar = self.sat_model.bstar
        self.B = 2*self.bstar/(2.461*10**(-5)*6378.135)
        self.ndot = self.sat_model.ndot*1440**2/(2*pi)
        self.id_launch_data = self.sat_model.intldesg
        self.element_number = self.sat_model.elnum
        self.satnumber = self.sat_model.satnum
        print(self.name + " TLE found!")
        self.RAAN = self.RAAN0
        self.w = self.w0
        self.MA = self.MA0
        self.theta = self.MA
        G = 6.67408*10**(-11)                          # Gravitational constant
        Mt = 5.9722*10**24                             # Earth mass
        self.t0 = self.epoch_day - int(self.epoch_day)
        self.t0 = self.t0*24*3600                      # Initial time in seconds
        self.updateGST0()                              # Get Greenwich sideral time
        self.a = (G*Mt/self.n**2)**(1/3)               # Semi-major axis
        self.v_peri = matrix([[0.0], [0.0]])           # Velocity in the perifocal frame
        self.v_iner = matrix([[0.0], [0.0], [0.0]])    # Velocity in the inertial frame
        self.r_iner = matrix([[0.0], [0.0], [0.0]])    # Position in the inertial frame
        self.changePlanet()                            # Sets all the planet's parameters
        self.updateOrbitalParameters()
        self.getLat()
        self.getLng()
        self.getAlt()

    def __call__(self):
        return self

    def readTLE(self, file_name):
        """
        Search for the satellite name inside the file_name. If found,
        sets the next two lines as the two line elements of this
        satellite.
        """
        with open(file_name, 'r') as f:
            for count, line in enumerate(f):
                if (line.strip() == self.name):
                    self.line1 = next(f).strip()
                    self.line2 = next(f).strip()
                    break

    def setMu(self, mu):
        """
        Sets the standard gravitational parameter.

        Parameters
        ----------
        mu : float
             The new standard gravitational parameter.
        """
        self.mu = mu

    def setInclination(self, incl):
        """
        Sets the inclination of the orbit.

        Parameters
        ----------
        incl : float
               The new satellite's orbit inclination in radians.
        """
        self.incl = incl

    def setRAAN(self, RAAN):
        """
        Sets the orbit's right ascension of the ascending node.

        Parameters
        ----------
        RAAN : float
               The new Right Ascension of the Ascending Node in radians.
        """
        self.RAAN = RAAN

    def setRAAN0(self, RAAN0):
        """
        Sets the orbit's initial RAAN.

        Parameters:
        -----------
        RAAN0 : float
                The initial Right Ascension of the Ascending Node in
                radians.
        """
        self.RAAN0 = RAAN0

    def setArgPerigee(self, w):
        """
        Sets the orbit's argument of the perigee.

        Parameters
        ----------
        w : float
            The new orbit's argument of the perigee in radians.
        """
        self.w = w

    def setArgPerigee0(self, w0):
        """
        Sets the orbit's initial argument of the perigee.

        Parameters
        ----------
        w0 : float
             The initial orbit's argument of the perigee in radians.
        """
        self.w0 = w0

    def setEccentricity(self, e):
        """
        Sets the orbit's eccentricity.

        Parameters
        ----------
        e : float
            The new orbit's eccentricity.
        """
        self.e = e

    def setMeanAnomaly(self, MA):
        """
        Sets the orbit's mean anomaly.

        Parameters
        ----------
        MA : float
             The new orbit's mean anomaly in radians.
        """
        self.MA = MA

    def setMeanAnomaly0(self, MA0):
        """
        Sets the orbit's initial mean anomaly.

        Parameters
        ----------
        MA0 : float
              The initial orbit's mean aomaly in radians.
        """
        self.MA0 = MA0

    def setMeanMotion(self, n):
        """
        Sets the satellite's mean motion.

        Parameters
        ----------
        n = float
            The satellite's mean motion.
        """
        self.n = n

    def setTrueAnomaly(self, theta):
        """
        Sets the orbit's true anomaly.

        Parameters
        ----------
        theta : float
                The new orbit's true anomaly in radians.
        """
        self.theta = theta

    def setSemiMajorAxis(self, a):
        """
        Sets the orbit's semi-major axis.

        Parameters
        ----------
        a = float
            The new orbit's semi-major axis in meters.
        """
        self.a = a

    def setSemilatusRectum(self, p):
        """
        Sets the orbit's semilatus rectum.

        Parameters
        ----------
        p : float
            The new orbit's semilatus rectum in meters.
        """
        self.p = p

    def setBallisticCoeff(self, B):
        """
        Sets the ballistic coefficient of the satellite.

        Parameters
        ----------
        B : float
            Ballistic coefficient.
        """
        self.B = B

    def setSpecAngMomentum(self, h):
        """
        Sets the satellite's specific relative angular momentum.

        Parameters
        ----------
        h : float
            The new specific relative angular momentum in squared meters
            per second.
        """
        self.h = h

    def setName(self, name):
        """
        Sets the satellite's name.

        Parameters
        ----------
        name : str
               The satellite's new name.
        """
        self.name = name
    
    def setCategory(self, cat):
        """
        Sets the satellite's category.

        Parameters
        ----------
        cat : str
              The satellite's new category.
        """
        self.cat = cat

    def getCategory(self):
        """Returns the satellite's category."""
        return self.cat

    def getInclination(self):
        """Returns the inclination of the orbit."""
        return self.incl

    def getRAAN(self):
        """Returns the orbit's RAAN."""
        return self.RAAN

    def getArgPerigee(self):
        """Returns the argument of the perigee of the orbit."""
        return self.w
    
    def getEccentricity(self):
        """Returns the eccentricity of the orbit."""
        return self.e
    
    def getSemiMajorAxis(self):
        """Returns the orbit's semi-major axis."""
        return self.a

    def getAnomaly(self):
        """Returns the orbit's true anomaly."""
        return self.theta

    def getMeanAnomaly(self):
        """Returns the orbit's mean anomaly."""
        return self.MA

    def getSpecAngMomentum(self):
        """
        Returns the satellite's specific relative angular momentum.
        """
        return self.h

    def getSpeed(self):
        """Returns the speed of the satellite."""
        mu_div_h = self.mu/self.h
        v_peri_p = -sin(self.theta)*mu_div_h
        v_peri_q = (self.e + cos(self.theta))*mu_div_h
        self.v = sqrt(v_peri_p**2 + v_peri_q**2)
        return self.v

    def getPerifocalVel(self):
        """
        Returns the satellite's velocity with respect to the perifocal
        frame.
        """
        mu_div_h = self.mu/self.h
        self.v_peri[0,0] = -sin(self.theta)*mu_div_h
        self.v_peri[1,0] = (self.e + cos(self.theta))*mu_div_h
        return self.v_peri[0,0], self.v_peri[1,0]

    def getInertialVel(self):
        """
        Calculates the velocity relative to its Perifocal Frame and
        transforms it to the inertial frame (Geocentric Equatorial
        Frame).

        Returns
        -------
        self.v_iner
            The velocity with respect to the inertial frame.
        """
        return self.v_iner

    def getXYZ(self):
        """
        Returns the position of the satellite with respect to the
        inertial frame.
        """
        return self.r_iner

    def getLat(self, rad2deg=180/pi):
        """
        Returns the last latitude of the satellite.
        """
        self.lat = arcsin(self.z/self.r)*rad2deg
        return self.lat

    def getLng(self, rad2deg=180/pi, twopi=2*pi, date=None):
        """
        Returns the last longitude of the satellite.
        """
        tnow = self.getTnow(date)
        RAAN_s = arctan2(self.y, self.x)
        lng = (RAAN_s - self.GST0 - self.P_w*tnow)  % twopi
        self.lng = ((lng > pi)*(lng - twopi) + (lng <= pi)*lng)*rad2deg
        return self.lng

    def getAlt(self):
        """
        Calculates the satellite's altitude considering the planet's
        radius at the current coordinates.

        Returns
        -------
        float
            Altitude in meters
        """
        self.alt = self.r - self.getPlanetRadius()
        return self.alt

    def getDMY(self):
        """
        Calculates the day, month and year from the TLE's epoch data and
        returns the result.
        """
        year = self.epoch_year
        ep_day = self.epoch_day
        if (year % 4 == 0):
            if (year % 100 == 0):
                if (year % 400 == 0):
                    day, month = self.leapYearDM(ep_day)
                else:
                    day, month = self.notLeapYearDM(ep_day)
            else:
                day, month = self.leapYearDM(ep_day)
        else:
            day, month = self.notLeapYearDM(ep_day)
        return day, month, year

    def leapYearDM(self, ep_day):
        """
        Calculates the days and months from an epoch day parameter,
        assuming it is a leap year.

        Parameters
        ----------
        ep_day : float
                 Satellite's epoch day from TLE.
        """
        month = 1 + (ep_day > 31) + (ep_day > 60) + (ep_day > 91)
        month += (ep_day > 121) + (ep_day > 152) + (ep_day > 182)
        month += (ep_day > 213) + (ep_day > 244) + (ep_day > 274)
        month += (ep_day > 305) + (ep_day > 335)

        day = ep_day - (ep_day > 31)*31 - (ep_day > 59)*29
        day = day - (ep_day > 90)*31 - (ep_day > 120)*30
        day = day - (ep_day > 151)*31 - (ep_day > 181)*30
        day = day - (ep_day > 212)*31 - (ep_day > 243)*31
        day = day - (ep_day > 273)*30 - (ep_day > 304)*31
        day = day - (ep_day > 334)*30
        return day, month

    def notLeapYearDM(self, ep_day):
        """
        Calculates the days and months from an epoch day parameter,
        assuming it is not a leap year.

        Parameters
        ----------
        ep_day : float
                 Satellite's epoch day from TLE.
        """
        month = 1 + (ep_day > 31) + (ep_day > 59) + (ep_day > 90)
        month += (ep_day > 120) + (ep_day > 151) + (ep_day > 181)
        month += (ep_day > 212) + (ep_day > 243) + (ep_day > 273)
        month += (ep_day > 304) + (ep_day > 334)

        day = ep_day - (ep_day > 31)*31 - (ep_day > 59)*28
        day = day- (ep_day > 90)*31 - (ep_day > 120)*30
        day = day - (ep_day > 151)*31 - (ep_day > 181)*30
        day = day - (ep_day > 212)*31 - (ep_day > 243)*31
        day = day - (ep_day > 273)*30 - (ep_day > 304)*31
        day = day - (ep_day > 334)*30
        return day, month

    def getGST(self, D, M, Y):
        """
        Returns the Greenwich Sidereal Time.

        Parameters
        ----------
        D : int
            Days since the beginning of the month.
        M : int
            Current month.
        Y : int
            Current year.

        Returns
        -------
        GST0
            The Greenwich Sidereal Time at t0.
        """
        JD = 367*Y - int(7/4*(Y + int((M + 9)/12))) + int(275*M/9) + D + 1721013.5 # Julian day
        T0 = (JD - 2451545)/36525
        GST0 = (100.4606184 + 36000.77004*T0 + 0.000387933*T0**2 - 2.583*10**(-8)*T0**3)*pi/180
        return GST0

    def getCurrentTimeInSeconds(self, date):
        """
        Returns the number of seconds since the beginning of the day.

        Parameters
        ----------
        date : datetime.utcnow()
               Current date.
        """
        seconds = date.hour*3600 + date.minute*60 + date.second
        return seconds + date.microsecond*0.000001

    def month2days(self, date):
        """
        Returns the number of months in days.

        Parameters
        ----------
        date : datetime
               Date from the datetime library.
        """
        days = date.day
        month = date.month
        year = date.year
        days = days + (month > 1)*31 + (month > 2)*28 + (month > 3)*31 + (month > 4)*30
        days = days + (month > 5)*31 + (month > 6)*30 + (month > 7)*31 + (month > 8)*31
        days = days + (month > 9)*30 + (month > 10)*31 + (month > 11)*30
        if (year % 4 == 0):
            if (year % 100 == 0):
                if (year % 400 == 0):
                    days = days + (month > 2)
            else:
                days = days + (month > 2)
        return days

    def M2E(self, M):
        """
        Transforms the mean anomaly to the eccentric anomaly.

        Parameters
        ----------
        M : float
            The mean anomaly in radians.
        """
        E0 = 0
        E = M
        while (abs(E - E0) > 0.0000001):
            E0 = E
            E = M + self.e*sin(E0)
        return E

    def E2theta(self, E):
        """
        Transforms the eccentric anomaly to the true anomaly.

        Parameters
        ----------
        E : float
            The eccentric anomaly in radians.
        """
        theta = 2*arctan(sqrt((1 + self.e)/(1 - self.e))*tan(E/2))
        return theta

    def getPlanetRadius(self):
        """
        Returns the planet radius considering the polar radius and the
        equatorial radius, given the current latitude and longitude.
        """
        cos_lat = cos(self.lat)
        sin_lat = sin(self.lat)
        aux = ((self.Eq_r**2*cos_lat)**2 + (self.Po_r**2*sin_lat)**2)
        radius = sqrt(aux/((self.Eq_r*cos_lat)**2 + (self.Po_r*sin_lat)**2))
        return radius

    def getPeriod(self):
        """
        Returns
        -------
        float
            The period of the satellite's orbit.
        """
        return 2*pi/self.n

    def getCoverage(self):
        """
        Returns
        -------
        float
            The angle of the coverage considering the planet's radius at
            the satellite's coordinates and the satellite's altitude. To
            return the angle in degrees, it is multiplied by
            57.29577951308232 which is 180/pi.
            self.r is the Planet's radius plus the satellite's altitude.
        """
        radius = self.getPlanetRadius()
        ang = arccos(radius/self.r)*57.29577951308232
        return ang

    def getComCoverage(self, E_r, c=299792458):
        """
        Calculates the angle of the coverage of the comunication link
        considering the transmit and receive power, the antenna gain,
        the sensibility and the losses.

        Parameters
        ----------
        E_r : float
              The Earth radius in meters.
        c   : int
              The speed of light in meters per second.

        Returns
        -------
        float
            The angle of the coverage in degrees.
        """
        # FSPL = Ptx + Ant_tx + Ant_rx - Cable + Sens - margin
        #SUCHAI
        alt = self.getAlt()
        FSPL = 30 + 2 + 18.9 - 1 + 117 - 12
        d = 10**(FSPL/20)*(c/self.freq)/(4*pi)
        ang = arccos((E_r**2 + (E_r + alt)**2 - d**2)/(2*E_r*(E_r + alt)))*180/pi
        return ang

    def getTnow(self, date=None):
        """
        Calculates the number of seconds since the last TLE data.

        Parameters
        ----------
        date : datetime
               The date used to calculate the number of seconds.

        Returns
        -------
        float
            Number of seconds.
        """
        if (date is None):
            date = datetime.utcnow()              # Use current time in UTC.
        dayinsec = 86400                          # Day in seconds
        tnow = self.getCurrentTimeInSeconds(date) # Current time in seconds
        days = self.month2days(date)              # Days in current time
        daysdiff = days - int(self.epoch_day)     # Difference between TLE date and current date
        return tnow + daysdiff*dayinsec           # Time in seconds from TLE to present time

    def getTrayectory(self, T, dt, date=None):
        """
        Calculates the future trayectory of the satellite using the
        SGP4 library.

        Parameters
        ----------
        T    : int
               Period of time in seconds to be calculated.
        dt   : int
               Time step in seconds.
        date : datetime
               Date from which the trayectory is calculated.
        """
        self.tray_lat[:] = []
        self.tray_lng[:] = []
        self.tray_alt[:] = []
        if (date is None):
            date = datetime.utcnow()
        current_date = date
        for t in range(0, T, dt):
            self.updateOrbitalParameters(date)
            date += timedelta(seconds=dt)
            self.tray_lat.append(self.getLat())
            self.tray_lng.append(self.getLng(date=date))
            self.getAlt()
            self.tray_alt.append(self.alt)
        self.updateOrbitalParameters(current_date)
        return self.tray_lat, self.tray_lng

    def changePlanet(self, M=5.9722*10**24, P_r=6371000, Eq_r=6378000, Po_r=6356000, J2=0.00108263, P_w=7.29211505*10**(-5)):
        """
        Changes the planet's parameters.

        Parameters
        ----------
        M    : float
               The planet's mass in kilograms.
        P_r  : int
               The planet's mean radius in meters.
        Eq_r : int
               The planet's equatorial radius in meters.
        Po_r : int
               The planet's polar radius in meters.
        J2   : float
               The planet's second degree harmonic model.
        P_w  : float
               The planet's angular velocity in radians per second.
        """
        G = 6.67408*10**(-11)              # Gravitational constant
        self.mu = G*M                      # Gravitational parameter
        self.n = sqrt(self.mu/(self.a**3)) # Mean motion
        self.p = self.a*(1 - self.e**2)
        self.h = sqrt(self.p*self.mu)
        self.P_r = P_r
        self.Eq_r = Eq_r
        self.Po_r = Po_r
        self.Er_Pr2 = (Eq_r/Po_r)**2
        self.J2 = J2
        self.P_w = P_w

    def updateOrbitalParameters(self, date=None):
        """
        Updates all the orbital parameters with the SGP4 propagation
        model.

        date : datetime
               The orbital elements are updated for this date.
        """
        if (date is None):
            date = datetime.utcnow()
        year = date.year
        month = date.month
        day = date.day
        hour = date.hour
        minute = date.minute
        second = date.second + date.microsecond*0.000001
        pos, vel = self.sat_model.propagate(year, month, day,
                                            hour, minute, second)
        p, a, e, i, raan, w, theta, m, argl, tlon, lonp = rv2coe(pos, vel,
                                                           self.mu*0.000000001)
        self.x = pos[0]*1000
        self.y = pos[1]*1000
        self.z = pos[2]*1000
        self.r = sqrt(self.x**2 + self.y**2 + self.z**2)
        self.r_iner[0,0] = self.x
        self.r_iner[1,0] = self.y
        self.r_iner[2,0] = self.z
        self.v_iner[0,0] = vel[0]*1000
        self.v_iner[1,0] = vel[1]*1000
        self.v_iner[2,0] = vel[2]*1000
        self.p = p*1000
        self.a = a*1000
        self.e = e
        self.incl = i
        self.RAAN = raan
        self.w = w
        self.theta = theta
        self.MA = m
        self.n = sqrt(self.mu/(self.a**3))
        self.h = sqrt(self.a*self.mu*(1 - self.e**2))

    def updateGST0(self):
        """
        Returns the Greenwich Sidereal Time of the TLE data.
        """
        D, Month, Y = self.getDMY()
        self.GST0 = self.getGST(int(D), Month, Y)

    def updateEpoch(self, date=None):
        """
        Updates the epoch day, year and the GST0 to now or to the
        received date.

        Parameters
        ----------
        date : datetime
               Date received to update the epoch.
        """
        if (date is None):
            date = datetime.utcnow()
        tnow = self.getCurrentTimeInSeconds(date)
        days = self.month2days(date)
        self.epoch_day = days + tnow/86400
        self.t0 = tnow
        self.epoch_year = date.year
        self.updateGST0()

    def getTLE(self):
        """
        Returns
        -------
        str
            The satellite's TLE as a string.
        """
        return "{}\n{}\n{}".format(self.name, self.line1, self.line2)

    def checksum(self, line):
        """
        Calculates the line's checksum.

        line : str
               String used to calculate its checksum.
        """
        check = 0
        for char in line[:-1]:
            if char.isdigit():
                check += int(char)
            if char == "-":
                check += 1
        return check % 10

    def createTLE(self, date):
        """
        Creates a new TLE for the received date.

        Parameters
        ----------
        date : datetime
               Date used to obtain the epoch of the TLE.
        """
        rad2deg = 180/pi
        self.updateEpoch(date)
        aux="{:+.9f}".format(self.ndot)
        ndot = "{}{}".format(aux[0],aux[2:-1])
        #aux = "{:+.6f}".format(self.bstar*10000)
        #BSTAR = "{}{}-4".format(aux[0],aux[3:-1])
        aux="{:+.6f}".format(self.B*2.461*10**(-5)*6378.135*0.5*100)
        BSTAR = "{}{}-2".format(aux[0],aux[3:-1])
        aux="{:.8f}".format(self.e)
        e = "{}".format(aux[2:-1])
        tle_num = "{:4d}".format(self.element_number)
        epoch_year = self.epoch_year - 2000
        line1 = "1 {}U {:9}{}{:012.8f} {} +{} {} 0 {}7".format(self.satnumber,
                                                           self.id_launch_data,
                                                           epoch_year,
                                                           self.epoch_day,
                                                           ndot,
                                                           "00000-0",
                                                           BSTAR,
                                                           tle_num)
        line2 = "2 {} {:8.4f} {:8.4f} {} {:8.4f} {:08.4f} {:02.8f}{}7".format(self.satnumber,
                                                             self.incl*rad2deg,
                                                             self.RAAN*rad2deg,
                                                             e,
                                                             self.w*rad2deg,
                                                             self.MA*rad2deg,
                                                             86400/self.getPeriod(),
                                                             "    0")
        checksum1 = self.checksum(line1)
        checksum2 = self.checksum(line2)
        line1 = "{}{}".format(line1[0:-1], checksum1)
        line2 = "{}{}".format(line2[0:-1], checksum2)
        self.sat_model = twoline2rv(line1, line2, wgs72)
        return self.name, line1, line2
