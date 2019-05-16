"""
                                Pypredict
    Orbit prediction software. Displays the satellites' position and
    orbital parameters in real time. Simulates satellite localization
    and deployment.
    
    Copyright (C) 2018-2019, Matías Vidal Valladares, matvidal.
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
from numpy import abs, arange, arccos, arcsin, arctan, arctan2, array, cos, matrix, pi, sin, sqrt, tan
from node import Node
from datetime import datetime
from pyorbital import tlefile
from Atmospheric_data import denUSSA76
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from sgp4.ext import rv2coe

class Sat(Node):
    __slots__ = ["cat", "incl", "RAAN0", "RAAN", "e", "w0", "w", "MA0", "MA", "n",
                 "epoch_year", "epoch_day", "theta", "GST0", "a", "mu", "t0", "p", 
                 "x", "y", "z", "Eq_r", "Po_r", "P_r", "Er_Pr2", "J2", "P_w", "r",
                 "h", "v", "tray_lat", "tray_lng", "tlast", "B", "bstar", "v_atm",
                 "v_peri", "v_iner", "r_iner", "Rot_Mat", "r_vect0", "sat"]
    def __init__(self, name="", lat=0, lng=0, alt=0, freq=437225000, tle=None, cat=""):
        """
        Parameters
        ----------
        name : str, optional
            The satellite's name (default is )
        lat : float, optional
            The satellite's, latitude (default is 0)
        lng : float, optional
            The satellite's longitude (default is 0)
        alt : float, optional
            The satellite's altitude in meters (default is 0)
        freq : int, optional
            The satellite's frequency in Hz (default is 437225000)
        tle : Tle, optional
            The satellite's TLE (default is None)
        cat : str
            The satellite's category (default is )
        """
        deg2rad = pi/180
        self.name = name
        self.lat = lat
        self.lng = lng
        self.alt = alt
        self.freq = freq
        self.cat = cat
        self.tray_lat = []
        self.tray_lng = []
        dayinsec = 24*3600                                        # Day in seconds
        if (tle is not None):
            self.incl = tle.inclination*deg2rad
            self.RAAN0 = tle.right_ascension*deg2rad
            self.e = tle.excentricity
            self.w0 = tle.arg_perigee*deg2rad
            self.MA0 = tle.mean_anomaly*deg2rad
            rpd2radps = (2*pi/dayinsec)                           # Revolutions per day to radians per second
            self.n = tle.mean_motion*rpd2radps                    # Radians per second
            self.epoch_year = tle.epoch_year
            self.epoch_day = tle.epoch_day
            self.bstar = tle.bstar
            self.B = 2*self.bstar/(2.461*10**(-5)*6378.135)
            self.sat = twoline2rv(tle.line1, tle.line2, wgs72)
            print(self.name + " TLE found!")
        else:
            self.incl = 97.39*deg2rad
            self.RAAN0 = 281.63*deg2rad
            self.e = 0.001
            self.w0 = 248*deg2rad
            self.MA0 = 112*deg2rad
            self.n = 0.0011
            self.B = 2.2*0.0025*0.0085/0.08
            self.updateEpoch()
            print("No TLE found, used default parameters instead")
        self.RAAN = self.RAAN0
        self.w = self.w0
        self.MA = self.MA0
        self.theta = self.MA
        G = 6.67408*10**(-11)                          # Gravitational constant
        Mt = 5.9722*10**24                             # Earth mass
        self.t0 = self.epoch_day - int(self.epoch_day)
        self.t0 = self.t0*dayinsec                     # Initial time in seconds
        self.tlast = self.t0                           # Time since last update, since t0
        self.updateGST0()                              # Get Greenwich sideral time
        self.a = (G*Mt/self.n**2)**(1/3)               # Semi-major axis
        self.v_atm = matrix([[0.0], [0.0], [0.0]])     # Velocity of the atmosphere
        self.v_peri = matrix([[0.0], [0.0]])           # Velocity in the perifocal frame
        self.v_iner = matrix([[0.0], [0.0], [0.0]])    # Velocity in the inertial frame
        self.r_iner = matrix([[0.0], [0.0], [0.0]])    # Position in the inertial frame
        self.r_vect0 = array([[0], [0]])               # Position in the plane
        self.Rot_Mat = matrix([[0.0, 0.0],             # Rotation matrix to move from
                               [0.0, 0.0],             # the perifocal frame to the
                               [0.0, 0.0]])            # inertial frame
        self.changePlanet()                            # Sets all the planet's parameters

    def __call__(self):
        return self

    def setMu(self, mu):
        """Sets the standard gravitational parameter."""
        self.mu = mu

    def setInclination(self, incl):
        """Sets the inclination of the orbit."""
        self.incl = incl

    def setRAAN(self, RAAN):
        """Sets the orbit's right ascension of the ascending node."""
        self.RAAN = RAAN

    def setRAAN0(self, RAAN0):
        """Sets the orbit's initial RAAN."""
        self.RAAN0 = RAAN0

    def setArgPerigee(self, w):
        """Sets the orbit's argument of the perigee."""
        self.w = w

    def setArgPerigee0(self, w0):
        """Sets the orbit's initial argument of the perigee."""
        self.w0 = w0

    def setEccentricity(self, e):
        """Sets the orbit's eccentricity."""
        self.e = e

    def setMeanAnomaly(self, MA):
        """Sets the orbit's mean anomaly."""
        self.MA = MA

    def setMeanAnomaly0(self, MA0):
        """Sets the orbit's initial mean anomaly."""
        self.MA0 = MA0

    def setMeanMotion(self, n):
        """Sets the satellite's mean motion."""
        self.n = n

    def setTrueAnomaly(self, theta):
        """Sets the orbit's true anomaly."""
        self.theta = theta

    def setSemiMajorAxis(self, a):
        """Sets the orbit's semi-major axis."""
        self.a = a

    def setSemilatusRectum(self, p):
        """Sets the orbit's semilatus rectum."""
        self.p = p

    def setMeanVelocity(self, n):
        """Sets the orbit's mean motion."""
        self.n = n

    def setBallisticCoeff(self, B):
        """Sets the ballistic coefficient of the satellite."""
        self.B = B

    def setSpecAngMomentum(self, h):
        """Sets the satellite's specific relative angular momentum."""
        self.h = h

    def setName(self, name):
        """Sets the satellite's name."""
        self.name = name
    
    def setCategory(self, cat):
        """Sets the satellite's category."""
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
        """Returns the satellite's specific relative angular momentum."""
        return self.h

    def getSpeed(self):
        """Returns the speed of the satellite."""
        self.v = sqrt(self.v_peri[0,0]**2 + self.v_peri[1,0]**2)
        return self.v

    def getPerifocalVel(self):
        """
        Returns the satellite's velocity with respect to the
        perifocal frame.
        """
        return self.v_peri[0,0], self.v_peri[1,0]

    def getInertialVel(self):
        """
        Calculates the velocity relative to its Perifocal Frame and 
        transforms it to the inertial frame (Geocentric Equatorial Frame).
        Returns the velocity with respect to the inertial frame.
        """
        return self.v_iner

    def getXYZ(self):
        """
        Returns the position of the satellite with respect to
        the inertial frame.
        """
        return self.r_iner

    def getLat(self, rad2deg=180/pi):
        """
        Returns the last latitude of the satellite.
        """
        self.lat = arcsin(self.z/self.r)*rad2deg
        return self.lat

    def getLng(self, rad2deg=180/pi, twopi=2*pi):
        """
        Returns the last longitude of the satellite.
        """
        tnow = self.getTnow()
        RAAN_s = arctan2(self.y, self.x)
        lng = (RAAN_s - self.GST0 - self.P_w*tnow)  % twopi
        self.lng = ((lng > pi)*(lng - twopi) + (lng <= pi)*lng)*rad2deg
        return self.lng

    def getDMY(self):
        """
        Calculates the day, month and year from the TLE's epoch data and
        returns the result.
        """
        Y = int(self.epoch_year) + 2000                                             # Year
        ep_day = self.epoch_day
        Month = 1 + (ep_day > 31) + (ep_day > 59) + (ep_day > 90) + (ep_day > 120)
        Month = Month + (ep_day > 151) + (ep_day > 181) + (ep_day > 212) + (ep_day > 243)
        Month = Month + (ep_day > 273) + (ep_day > 304) + (ep_day > 334)            # Month

        D = ep_day - (ep_day > 31)*31 - (ep_day > 59)*28 - (ep_day > 90)*31 - (ep_day > 120)*30
        D = D - (ep_day > 151)*31 - (ep_day > 181)*30 - (ep_day > 212)*31 - (ep_day > 243)*31
        D = D - (ep_day > 273)*30 - (ep_day > 304)*31 - (ep_day > 334)*30           # Day
        return D, Month, Y

    def getGST(self, D, M, Y):
        """
        Returns the Greenwich Sidereal Time.

        Parameters
        ----------
        D : int
            Days since the beginning of the month
        M : int
            Current month
        Y : int
            Current year
        """
        JD = 367*Y - int(7/4*(Y + int((M + 9)/12))) + int(275*M/9) + D + 1721013.5  # Julian day
        T0 = (JD - 2451545)/36525
        GST0 = (100.4606184 + 36000.77004*T0 + 0.000387933*T0**2 - 2.583*10**(-8)*T0**3)*pi/180
        return GST0

    def getCurrentTimeInSeconds(self, date):
        """
        Returns the number of seconds since the beginning of the day.

        Parameters
        ----------
        date : datetime.utcnow()
            Current date
        """
        return date.hour*3600 + date.minute*60 + date.second

    def month2days(self, month):
        """
        Returns the number of months in days.

        Parameters
        ----------
        month : int
            Number of months to be transformed to days
        """
        days = (month > 1)*31 + (month > 2)*28 + (month > 3)*31 + (month > 4)*30
        days = days + (month > 5)*31 + (month > 6)*30 + (month > 7)*31 + (month > 8)*31
        days = days + (month > 9)*30 + (month > 10)*31 + (month > 11)*30
        return days

    def M2E(self, M):
        """
        Transforms the mean anomaly to the eccentric anomaly.

        Parameters
        ----------
        M : float
            The mean anomaly in radians
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
            The eccentric anomaly in radians
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
        radius = sqrt(((self.Eq_r**2*cos_lat)**2 + (self.Po_r**2*sin_lat)**2)/((self.Eq_r*cos_lat)**2 + (self.Po_r*sin_lat)**2))
        return radius

    def getPeriod(self):
        """
        Returns the period of the satellite's orbit.
        """
        return 2*pi/self.n

    def getCoverage(self):
        """
        Returns the angle of the coverage considering planet's radius
        at the satellite's coordinates and the satellite's altitude.
        """
        radius = self.getPlanetRadius()
        ang = arccos(radius/(radius + self.alt))*180/pi
        return ang

    def getComCoverage(self, E_r):
        """
        Returns the angle of the coverage of the comunication link
        considering the transmit and receive power, the antenna
        gain, the sensibility and the losses.
        """
        # FSPL = Ptx + Ant_tx + Ant_rx - Cable + Sens - margin
        #SUCHAI
        FSPL = 30 + 2 + 18.9 - 1 + 117 - 12
        c = 299792458
        d = 10**(FSPL/20)*(c/self.freq)/(4*pi)
        ang = arccos((E_r**2 + (E_r + self.alt)**2 - d**2)/(2*E_r*(E_r + self.alt)))*180/pi
        return ang

    def getTnow(self):
        """
        Returns the number of seconds since the last TLE data.
        """
        dayinsec = 86400                                # Day in seconds
        date = datetime.utcnow()                        # Use current time in UTC.
        tnow = self.getCurrentTimeInSeconds(date)       # Current time in seconds
        days = self.month2days(date.month) + date.day   # Days in current time
        daysdiff = days - int(self.epoch_day)           # Difference between TLE date and current date in days
        return tnow + daysdiff*dayinsec                 # Time in seconds from TLE to present time

    def getDrag(self):
        """
        Returns the drag at the satellite's altitude.
        """
        alt = self.alt
        if (alt > 30000):
            rho = denUSSA76(alt*0.001)
        elif (alt > 25000):
            T = -131.21 + 0.00299*alt
            P = 2.488*((T + 273.1)/216.6)**-11.388
            rho = P/(0.2869*(T - 273.1))
        elif (alt > 11000):
            T = -56.46
            P = 22.65*exp(1.73 - 0.000157*alt)
            rho = P/(0.2869*(T - 273.1))
        elif (alt > 0):
            T = 15.04 - 0.00649*alt
            P = 101.29*((T + 273.1)/288.08)**5.256
            rho = P/(0.2869*(T - 273.1))
        else:
            print("{}{:0.2f}".format("Invalid altitude: ", alt))
            rho = 0
        return rho

    def updateWithDragEffect(self, tnow=None):
        if (tnow is None):
            tnow = self.getTnow()
        dt = tnow - self.tlast
        while (dt > self.getPeriod()):
            dt = dt - self.getPeriod()
            drag = self.getDrag()
            da = -2*pi*self.B*drag*self.a**2
            self.a = self.a + da
            self.n = sqrt(self.mu/(self.a**3))
            self.p = self.a*(1 - self.e**2)
            self.h = sqrt(self.p*self.mu)
            self.e = sqrt(1 - self.h**2/(self.a*self.mu))
            self.updateOrbitalParameters(tnow - dt)
        drag = self.getDrag()
        da = -2*pi*self.B*drag*self.a**2*dt/self.getPeriod()
        self.a = self.a + da
        self.n = sqrt(self.mu/(self.a**3))
        self.p = self.a*(1 - self.e**2)
        self.h = sqrt(self.p*self.mu)
        self.e = sqrt(1 - self.h**2/(self.a*self.mu))
        self.updateOrbitalParameters(tnow)
        self.tlast = tnow

    def getLocation(self, T, dt):
        tnow = self.getTnow()
        self.tray_lat[:] = []
        self.tray_lng[:] = []
        rad2deg = 180/pi                             # Radian to degrees
        twopi = 2*pi                                 # Two times pi
        aux = self.n*(self.P_r**2/self.p**2)*self.J2 # Calculate only one time
        cos_incl = cos(self.incl)
        sin_incl = sin(self.incl)
        for t in range(tnow, T + tnow, dt):
            RAAN = self.RAAN0 - 1.5*aux*cos_incl*(t - self.t0)
            w = self.w0 + 0.75*aux*(5*cos_incl**2 - 1)*(t - self.t0)
            MA = self.MA0 + (self.n + 0.75*aux*(sqrt(1 - self.e**2))*(2 - 3*sin_incl**2))*(t - self.t0)
            E = self.M2E(MA)
            theta = self.E2theta(E)

            cos_RAAN = cos(RAAN)
            sin_RAAN = sin(RAAN)
            cos_w = cos(w)
            sin_w = sin(w)
            cos_theta = cos(theta)
            sin_theta = sin(theta)

            self.Rot_Mat[0,0] = cos_RAAN*cos_w - sin_RAAN*cos_incl*sin_w
            self.Rot_Mat[0,1] = -sin_RAAN*cos_incl*cos_w - cos_RAAN*sin_w
            self.Rot_Mat[1,0] = cos_RAAN*cos_incl*sin_w + sin_RAAN*cos_w
            self.Rot_Mat[1,1] = cos_RAAN*cos_incl*cos_w - sin_RAAN*sin_w
            self.Rot_Mat[2,0] = sin_incl*sin_w
            self.Rot_Mat[2,1] = sin_incl*cos_w
            r0 = self.p/(1 + self.e*cos_theta)
            self.r_vect0[0,0] = r0*cos_theta
            self.r_vect0[1,0] = r0*sin_theta
            r_vectf = self.Rot_Mat*self.r_vect0
            r_s = sqrt(r_vectf.item(0)**2 + r_vectf.item(1)**2 + r_vectf.item(2)**2)
            
            #decl_s = arcsin(r_vectf.item(2)/r_s)
            self.tray_lat.append(arcsin(r_vectf.item(2)/r_s)*rad2deg)
            RAAN_s = arctan2(r_vectf.item(1),r_vectf.item(0))
            #lat = arctan(self.Er_Pr2*tan(decl_s)) % pi
            lng = (RAAN_s - self.GST0 - self.P_w*t)  % twopi
            #self.tray_lat.append(((lat > pi/2)*(lat - pi) + (lat <= pi/2)*lat)*rad2deg)
            self.tray_lng.append(((lng > pi)*(lng - twopi) + (lng <= pi)*lng)*rad2deg)
        return self.tray_lat, self.tray_lng

    def changePlanet(self, M=5.9722*10**24, P_r=6371000, Eq_r=6378000, Po_r=6356000, J2=0.00108263, P_w=7.29211505*10**(-5)):
        G = 6.67408*10**(-11)                       # Gravitational constant
        self.mu = G*M                               # Gravitational parameter
        self.n = sqrt(self.mu/(self.a**3))          # Mean motion
        self.p = self.a*(1 - self.e**2)
        self.h = sqrt(self.p*self.mu)
        self.P_r = P_r                              # Planet radius
        self.Eq_r = Eq_r                            # Equatorial radius
        self.Po_r = Po_r                            # Polar radius
        self.Er_Pr2 = (Eq_r/Po_r)**2                # (Equatorial radius/Polar radius)^2
        self.J2 = J2                                # The planet's second degree harmonic model
        self.P_w = P_w                              # The planet's angular velocity
        if (self.tlast == self.t0):
            self.updateOrbitalParameters(self.t0)
        #self.updateOrbitalParameters3()

        #tnow = self.getTnow()
        #for t in arange(int(self.t0)+1, tnow, 0.5):
        #    self.updateOrbitalParameters2(t)
        
        #self.updateWithDragEffect()

    def updateOrbitalParameters(self, tnow=None):
        if (tnow is None):
            tnow = self.getTnow()
        rad2deg = 180/pi                                # Radian to degrees
        twopi = 2*pi                                    # Two times pi
        aux = self.n*(self.P_r**2/self.p**2)*self.J2    # Calculate only one time
        cos_incl = cos(self.incl)
        sin_incl = sin(self.incl)

        self.RAAN = self.RAAN0 - 1.5*aux*cos_incl*(tnow - self.t0)
        self.w = self.w0 + 0.75*aux*(5*cos_incl**2 - 1)*(tnow - self.t0)
        self.MA = (self.MA0 + (self.n + 0.75*aux*(sqrt(1 - self.e**2))*(2 - 3*sin_incl**2))*(tnow - self.t0)) % twopi
        E = self.M2E(self.MA)
        self.theta = self.E2theta(E) % twopi
        cos_RAAN = cos(self.RAAN)
        sin_RAAN = sin(self.RAAN)
        cos_w = cos(self.w)
        sin_w = sin(self.w)
        cos_theta = cos(self.theta)
        sin_theta = sin(self.theta)

        #Rot_Mat = matrix((((cos_RAAN*cos_w - sin_RAAN*cos_incl*sin_w), (-sin_RAAN*cos_incl*cos_w - cos_RAAN*sin_w), (sin_RAAN*sin_incl)),
        #                  ((cos_RAAN*cos_incl*sin_w + sin_RAAN*cos_w), (cos_RAAN*cos_incl*cos_w - sin_RAAN*sin_w), (-cos_RAAN*sin_incl)),
        #                  ((sin_incl*sin_w), (sin_incl*cos_w), (cos_incl))))
        self.Rot_Mat[0,0] = cos_RAAN*cos_w - sin_RAAN*cos_incl*sin_w
        self.Rot_Mat[0,1] = -sin_RAAN*cos_incl*cos_w - cos_RAAN*sin_w
        self.Rot_Mat[1,0] = cos_RAAN*cos_incl*sin_w + sin_RAAN*cos_w
        self.Rot_Mat[1,1] = cos_RAAN*cos_incl*cos_w - sin_RAAN*sin_w
        self.Rot_Mat[2,0] = sin_incl*sin_w
        self.Rot_Mat[2,1] = sin_incl*cos_w
        r0 = self.p/(1 + self.e*cos_theta)
        self.r_vect0[0,0] = r0*cos_theta
        self.r_vect0[1,0] = r0*sin_theta
        r_vectf = self.Rot_Mat*self.r_vect0
        self.x = r_vectf.item(0)
        self.y = r_vectf.item(1)
        self.z = r_vectf.item(2)
        self.r = sqrt(self.x**2 + self.y**2 + self.z**2)
        #decl_s = arcsin(self.z/self.r)
        self.lat = arcsin(self.z/self.r)*rad2deg
        RAAN_s = arctan2(self.y, self.x)

        #lat = arctan(self.Er_Pr2*tan(decl_s)) % pi
        lng = (RAAN_s - self.GST0 - self.P_w*tnow)  % twopi

        #self.lat = ((lat > pi/2)*(lat - pi) + (lat <= pi/2)*lat)*rad2deg
        self.lng = ((lng > pi)*(lng - twopi) + (lng <= pi)*lng)*rad2deg
        self.alt = self.r - self.getPlanetRadius()
        self.v_peri[0,0] = -self.mu*sin_theta/self.h
        self.v_peri[1,0] = self.mu*(self.e + cos_theta)/self.h
        self.v_iner = self.Rot_Mat*self.v_peri
        #self.vt = self.h/self.r
        #self.vr = self.mu*self.e*sin(self.theta)/self.h
        #self.v = sqrt(self.v_p**2 + self.v_q**2)

    def updateOrbitalParameters2(self, tnow=None, rad2deg=180/pi, twopi=2*pi):
        if (tnow is None):
            tnow = self.getTnow()
        dt = tnow - self.tlast

        dh, de, de2, dtheta, dRAAN, di, dw, da = self.getDeltas()
        self.h = self.h + dh*dt
        self.e = self.e + de*dt
        self.theta = (self.theta + dtheta*dt) % twopi
        self.RAAN = (self.RAAN + dRAAN*dt) % twopi
        self.incl = self.incl + di*dt
        self.w = (self.w + dw*dt) % twopi
        self.a = self.h**2/(self.mu*(1 - self.e**2))
        self.a = self.a + da*dt
        self.e = abs(self.e + de2*dt)
        self.p = self.a*(1 - self.e**2)
        self.h = sqrt(self.a*self.mu*(1 - self.e**2))
        self.n = sqrt(self.mu/(self.a**3))

        cos_incl = cos(self.incl)
        sin_incl = sin(self.incl)
        cos_RAAN = cos(self.RAAN)
        sin_RAAN = sin(self.RAAN)
        cos_w = cos(self.w)
        sin_w = sin(self.w)
        cos_theta = cos(self.theta)
        sin_theta = sin(self.theta)
        
        self.Rot_Mat[0,0] = cos_RAAN*cos_w - sin_RAAN*cos_incl*sin_w
        self.Rot_Mat[0,1] = -sin_RAAN*cos_incl*cos_w - cos_RAAN*sin_w
        self.Rot_Mat[1,0] = cos_RAAN*cos_incl*sin_w + sin_RAAN*cos_w
        self.Rot_Mat[1,1] = cos_RAAN*cos_incl*cos_w - sin_RAAN*sin_w
        self.Rot_Mat[2,0] = sin_incl*sin_w
        self.Rot_Mat[2,1] = sin_incl*cos_w
        r0 = self.p/(1 + self.e*cos_theta)
        r_vec0 = r0*cos_theta
        r_vec1 = r0*sin_theta
        self.r_iner[0,0] = self.Rot_Mat[0,0]*r_vec0 + self.Rot_Mat[0,1]*r_vec1
        self.r_iner[1,0] = self.Rot_Mat[1,0]*r_vec0 + self.Rot_Mat[1,1]*r_vec1
        self.r_iner[2,0] = self.Rot_Mat[2,0]*r_vec0 + self.Rot_Mat[2,1]*r_vec1
        self.x = self.r_iner[0,0]
        self.y = self.r_iner[1,0]
        self.z = self.r_iner[2,0]
        self.r = sqrt(self.x**2 + self.y**2 + self.z**2)

        self.alt = self.r - self.getPlanetRadius()
        self.v_peri[0,0] = -self.mu*sin_theta/self.h
        self.v_peri[1,0] = self.mu*(self.e + cos_theta)/self.h
        #self.v_iner = self.Rot_Mat*self.v_peri
        self.v_iner[0,0] = self.Rot_Mat[0,0]*self.v_peri[0,0] + self.Rot_Mat[0,1]*self.v_peri[1,0]
        self.v_iner[1,0] = self.Rot_Mat[1,0]*self.v_peri[0,0] + self.Rot_Mat[1,1]*self.v_peri[1,0]
        self.v_iner[2,0] = self.Rot_Mat[2,0]*self.v_peri[0,0] + self.Rot_Mat[2,1]*self.v_peri[1,0]
        self.tlast = tnow 

    def updateOrbitalParameters3(self, date=None):
        if (date is None):
            date = datetime.utcnow()
        year = date.year
        month = date.month
        day = date.day
        hour = date.hour
        minute = date.minute
        second = date.second + date.microsecond*0.000001
        pos, vel = self.sat.propagate(year, month, day,
                hour, minute, second)
        p, a, e, i, raan, w, theta, m, argl, tlon, lonp = rv2coe(pos, vel,
                self.mu*0.000000001)
        self.x = pos[0]*1000
        self.y = pos[1]*1000
        self.z = pos[2]*1000
        self.r = sqrt(self.x**2 + self.y**2 + self.z**2)
        self.alt = self.r - self.getPlanetRadius()
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
        self.v_peri[0,0] = -self.mu*sin(theta)/self.h
        self.v_peri[1,0] = self.mu*(e + cos(theta))/self.h

    def updateGST0(self):
        D, Month, Y = self.getDMY()
        self.GST0 = self.getGST(int(D), Month, Y)

    def updateEpoch(self):
        date = datetime.utcnow()                        # Use current time in UTC.
        tnow = self.getCurrentTimeInSeconds(date)       # Current time in seconds
        days = self.month2days(date.month) + date.day   # Days in current time
        self.epoch_day = days + tnow/86400
        self.t0 = tnow
        self.tlast = self.t0
        self.epoch_year = date.year - 2000

    def J2_effect(self):
        u = self.w + self.theta
        sin_u = sin(u)
        sin_2u = sin(2*u)
        cos_incl = cos(self.incl)
        sin_incl = sin(self.incl)
        cos_theta = cos(self.theta)
        sin_theta = sin(self.theta)
        aux1 = 1.5*self.J2*self.mu*self.P_r**2/self.r**3 # Calculate only one time
        aux2 = aux1/self.h
        aux3 = aux2/self.e
        aux4 = self.h**2*(3*sin_incl**2*sin_u**2 - 1)/(self.mu*self.r)
        aux5 = sin_2u*sin_incl**2
        aux6 = 2 + self.e*cos_theta
        aux7 = 2*self.e*cos_incl**2*sin_u**2
        aux8 = self.h/self.r**2

        dh = -aux1*aux5
        de = aux4*sin_theta - aux5*(aux6*cos_theta + self.e)
        de = aux2*de
        dtheta = aux4*cos_theta + aux6*aux5*sin_theta
        dtheta = aux8 + aux3*dtheta
        dRAAN = -2*aux2*sin_u**2*cos_incl
        di = -0.5*aux2*sin_2u*sin(2*self.incl)
        dw = -aux4*cos_theta - aux6*aux5*sin_theta + aux7
        dw = aux3*dw
        return dh, de, dtheta, dRAAN, di, dw

    def drag_effect(self):
        rho = self.getDrag()
        self.v_atm[0,0] = -self.P_w*self.y
        self.v_atm[1,0] = self.P_w*self.x
        v_rel = self.v_iner - self.v_atm
        v_rel_s = v_rel[0,0]**2 + v_rel[1,0]**2 + v_rel[2,0]**2
        v_rel_s = sqrt(v_rel_s)
        aux = self.B*rho*v_rel_s

        da = -aux*v_rel_s**2/(self.a*self.n**2)
        de = -aux*(self.e + cos(self.theta))
        dw = -aux*sin(self.theta)/self.e
        return da, de, dw

    def getDeltas(self):
        dh, de, dtheta, dRAAN, di, dw = self.J2_effect()
        da, de2, dw2 = self.drag_effect()
        #de = de + de2
        dw = dw + dw2
        return dh, de, de2, dtheta, dRAAN, di, dw, da

    def createTLE(self, tle):
        aux="{:+.9f}".format(tle.mean_motion_derivative)
        mean_motion_derivative = "{}{}".format(aux[0],aux[2:-1])
        #aux = "{:+.6f}".format(tle.bstar*10000)
        #BSTAR = "{}{}-4".format(aux[0],aux[3:-1])
        aux="{:+.6f}".format(self.B*2.461*10**(-5)*6378.135*0.5*100)
        BSTAR = "{}{}-2".format(aux[0],aux[3:-1])
        aux="{:.8f}".format(self.e)
        e = "{}".format(aux[2:-1])
        aux = tle.element_number+1
        tle_num = "{:4d}".format(aux)
        TLE1 = "{} {}{} {}{}{} {}{:08.8f} {} +{} {} {} {}{}".format("1",
                                tle.satnumber,
                                "U",
                                tle.id_launch_year,
                                tle.id_launch_number,
                                tle.id_launch_piece,
                                self.epoch_year,
                                self.epoch_day,
                                mean_motion_derivative,
                                #tle.mean_motion_sec_derivative,
                                "00000-0",
                                BSTAR,
                                "0",
                                tle_num,
                                "7")
        TLE2 = "{} {} {:8.4f} {:08.4f} {} {:08.4f} {:08.4f} {:02.8f}{}{}".format("2",
                                                             tle.satnumber,
                                                             self.incl*180/pi,
                                                             self.RAAN*180/pi,
                                                             e,
                                                             self.w*180/pi,
                                                             self.MA*180/pi,
                                                             86400/self.getPeriod(),
                                                             "    0",
                                                             "7")
        print(self.name)
        print(TLE1)
        print(TLE2)
        self.sat = twoline2rv(TLE1, TLE2, wgs72)

