'''
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
'''
from numpy import abs, arccos, arcsin, arctan, arctan2, array, cos, matrix, pi, sin, sqrt, tan
from node import Node
from datetime import datetime
from pyorbital import tlefile

class Sat(Node):
    __slots__ = ["cat", "incl", "RAAN0", "RAAN", "e", "w0", "w", "MA0", "MA", "n",
                 "epoch_year", "epoch_day", "theta", "GST0", "a", "mu", "t0", "p", 
                 "x", "y", "z", "Eq_r", "Po_r", "P_r", "Er_Pr2", "J2", "P_w", "r",
                 "h", "vt", "vr", "v", "tray_lat", "tray_lng", "tlast", "B", "bstar"]
    def __init__(self, name="", lat=0, lng=0, alt=0, freq=437225000, tle=None, cat=""):
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
            self.B = 2*self.bstar/(2.461*10**(-5)*6371)
            print(self.name + " TLE found!")
        else:
            self.incl = 97.39*deg2rad
            self.RAAN0 = 281.63*deg2rad
            self.e = 0.001
            self.w0 = 248*deg2rad
            self.MA0 = 112*deg2rad
            self.n = 0.0011
            self.epoch_year = 18
            self.epoch_day = 222.64
            self.B = 2.2*0.0025*0.0085/0.08
            print("No TLE found, used default parameters instead")
        self.RAAN = self.RAAN0
        self.w = self.w0
        self.MA = self.MA0
        self.theta = self.MA
        G = 6.67408*10**(-11)                                     # Gravitational constant
        Mt = 5.9722*10**24                                        # Earth mass
        self.t0 = (self.epoch_day - int(self.epoch_day))*dayinsec # Initial time in seconds
        self.tlast = self.t0
        D, Month, Y = self.getDMY()                               # Get day, month and year from TLE data
        self.GST0 = self.getGST(int(D), Month, Y)                 # Get Greenwich sideral time
        self.a = (G*Mt/self.n**2)**(1/3)                          # Semi-major axis
        self.changePlanet()                                       # Sets all the planet's parameters

    def __call__(self):
        return self

    def setMu(self, mu):
        self.mu = mu

    def setInclination(self, incl):
        self.incl = incl

    def setRAAN(self, RAAN):
        self.RAAN = RAAN

    def setRAAN0(self, RAAN0):
        self.RAAN0 = RAAN0

    def setArgPerigee(self, w):
        self.w = w

    def setArgPerigee0(self, w0):
        self.w0 = w0

    def setEccentricity(self, e):
        self.e = e

    def setMeanAnomaly(self, MA):
        self.MA = MA

    def setMeanAnomaly0(self, MA0):
        self.MA0 = MA0

    def setTrueAnomaly(self, theta):
        self.theta = theta

    def setSemiMajorAxis(self, a):
        self.a = a

    def setSemilatusRectum(self, p):
        self.p = p

    def setMeanVelocity(self, n):
        self.n = n

    def setBallisticCoeff(self, B):
        self.B = B

    def setName(self, name):
        self.name = name
    
    def setCategory(self, cat):
        self.cat = cat

    def getCategory(self):
        return self.cat

    def getInclination(self):
        return self.incl

    def getRAAN(self):
        return self.RAAN

    def getArgPerigee(self):
        return self.w
    
    def getEccentricity(self):
        return self.e
    
    def getSemiMajorAxis(self):
        return self.a

    def getAnomaly(self):
        return self.theta

    def getMeanAnomaly(self):
        return self.MA

    def getSpecAngMomentum(self):
        return self.h

    def getSpeed(self):
        return self.v

    def getPerifocalVel(self):
        v_p = -self.mu*sin(self.theta)/self.h
        v_q = self.mu*(self.e + cos(self.theta))/self.h
        return v_p, v_q

    # Calculates the velocity relative to its Perifocal Frame and
    # transforms it to the Geocentric Equatorial Frame
    def getInertialVel(self):
        v_p, v_q = self.getPerifocalVel()
        sin_RAAN = sin(self.RAAN)
        cos_RAAN = cos(self.RAAN)
        sin_i = sin(self.incl)
        cos_i = cos(self.incl)
        sin_w = sin(self.w)
        cos_w = cos(self.w)
        v_x = v_p*(-sin_RAAN*cos_i*sin_w + cos_RAAN*cos_w)
        v_x -= v_q*(sin_RAAN*cos_i*cos_w + cos_RAAN*sin_w)
        v_y = v_p*(cos_RAAN*cos_i*sin_w + sin_RAAN*cos_w)
        v_y += v_q*(cos_RAAN*cos_i*cos_w - sin_RAAN*sin_w)
        v_z = v_p*sin_i*sin_w + v_q*sin_i*cos_w
        return v_x, v_y, v_z


    def getXYZ(self):
        return matrix([[self.x], [self.y], [self.z]])

    def getDMY(self):
        Y = int(self.epoch_year) + 2000                                                         # Year
        ep_day = self.epoch_day
        Month = 1 + (ep_day > 31) + (ep_day > 59) + (ep_day > 90) + (ep_day > 120)
        Month = Month + (ep_day > 151) + (ep_day > 181) + (ep_day > 212) + (ep_day > 243)
        Month = Month + (ep_day > 273) + (ep_day > 304) + (ep_day > 334)                        # Month

        D = ep_day - (ep_day > 31)*31 - (ep_day > 59)*28 - (ep_day > 90)*31 - (ep_day > 120)*30
        D = D - (ep_day > 151)*31 - (ep_day > 181)*30 - (ep_day > 212)*31 - (ep_day > 243)*31
        D = D - (ep_day > 273)*30 - (ep_day > 304)*31 - (ep_day > 334)*30                       # Day
        return D, Month, Y

    def getGST(self, D, M, Y):
        JD = 367*Y - int(7/4*(Y + int((M + 9)/12))) + int(275*M/9) + D + 1721013.5  # Julian day
        T0 = (JD - 2451545)/36525
        GST0 = (100.4606184 + 36000.77004*T0 + 0.000387933*T0**2 - 2.583*10**(-8)*T0**3)*pi/180
        return GST0

    def getCurrentTimeInSeconds(self, date):
        return date.hour*3600 + date.minute*60 + date.second

    def month2days(self, month):
        days = (month > 1)*31 + (month > 2)*28 + (month > 3)*31 + (month > 4)*30
        days = days + (month > 5)*31 + (month > 6)*30 + (month > 7)*31 + (month > 8)*31
        days = days + (month > 9)*30 + (month > 10)*31 + (month > 11)*30
        return days

    def M2E(self, M):
        E0 = 0
        E = M
        while (abs(E - E0) > 0.0000001):
            E0 = E
            E = M + self.e*sin(E0)
        return E

    def E2theta(self, E):
        theta = 2*arctan(sqrt((1 + self.e)/(1 - self.e))*tan(E/2))
        return theta

    def getPlanetRadius(self):
        cos_lat = cos(self.lat)
        sin_lat = sin(self.lat)
        radius = sqrt(((self.Eq_r**2*cos_lat)**2 + (self.Po_r**2*sin_lat)**2)/((self.Eq_r*cos_lat)**2 + (self.Po_r*sin_lat)**2))
        return radius

    def getPeriod(self):
        return 2*pi/self.n

    def getCoverage(self):
        radius = self.getPlanetRadius()
        ang = arccos(radius/(radius + self.alt))*180/pi
        return ang

    def getComCoverage(self, E_r):
        # FSPL = Ptx + Ant_tx + Ant_rx - Cable + Sens - margin
        #SUCHAI
        FSPL = 30 + 2 + 18.9 - 1 + 117 - 12
        c = 299792458
        d = 10**(FSPL/20)*(c/self.freq)/(4*pi)
        ang = arccos((E_r**2 + (E_r + self.alt)**2 - d**2)/(2*E_r*(E_r + self.alt)))*180/pi
        return ang

    def getTnow(self):
        dayinsec = 24*3600                              # Day in seconds
        date = datetime.utcnow()                        # Use current time in UTC.
        tnow = self.getCurrentTimeInSeconds(date)       # Current time in seconds
        days = self.month2days(date.month) + date.day   # Days in current time
        daysdiff = days - int(self.epoch_day)           # Difference between  TLE date and current date in days
        return tnow + daysdiff*dayinsec                 # Time in seconds from TLE to present time

    def getDrag(self):
        if (self.alt < 11000):
            T = 15.04 - 0.00649*self.alt
            P = 101.29*((T + 273.1)/288.08)**5.256
        elif (self.alt < 25000):
            T = -56.46
            P = 22.65*exp(1.73 - 0.000157*self.alt)
        elif (self.alt >= 25000):
            T = -131.21 + 0.00299*self.alt
            P = 2.488*((T + 273.1)/216.6)**-11.388
        else:
            print("{}{:0.2f}".format("Invalid altitude: ", self.alt))
            P = 0
            T = 0
        return P/(0.2869*(T - 273.1))

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

    def getLocation(self, tf, dt):
        tnow = self.getTnow()
        self.tray_lat[:] = []
        self.tray_lng[:] = []
        rad2deg = 180/pi                             # Radian to degrees
        twopi = 2*pi                                 # Two times pi
        aux = self.n*(self.P_r**2/self.p**2)*self.J2 # Calculate only one time
        cos_incl = cos(self.incl)
        sin_incl = sin(self.incl)
        for i in range(0, tf, dt):
            t = tnow + i
            j = int(i/dt)

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

            Rot_Mat = matrix((((cos_RAAN*cos_w - sin_RAAN*cos_incl*sin_w), (-sin_RAAN*cos_incl*cos_w - cos_RAAN*sin_w), (sin_RAAN*sin_incl)),
                              ((cos_RAAN*cos_incl*sin_w + sin_RAAN*cos_w), (cos_RAAN*cos_incl*cos_w - sin_RAAN*sin_w), (-cos_RAAN*sin_incl)),
                              ((sin_incl*sin_w), (sin_incl*cos_w), (cos_incl))))
            r0 = self.p/(1 + self.e*cos_theta)
            r_vect0 = array([[r0*cos_theta], [r0*sin_theta], [0]])
            r_vectf = Rot_Mat*r_vect0
            r_s = sqrt(r_vectf.item(0)**2 + r_vectf.item(1)**2 + r_vectf.item(2)**2)
            decl_s = arcsin(r_vectf.item(2)/r_s)
            RAAN_s = arctan2(r_vectf.item(1),r_vectf.item(0))

            lat = arctan(self.Er_Pr2*tan(decl_s)) % pi
            lng = (RAAN_s - self.GST0 - self.P_w*t)  % twopi

            self.tray_lat.append(((lat > pi/2)*(lat - pi) + (lat <= pi/2)*lat)*rad2deg)
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
        self.updateWithDragEffect()

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

        Rot_Mat = matrix((((cos_RAAN*cos_w - sin_RAAN*cos_incl*sin_w), (-sin_RAAN*cos_incl*cos_w - cos_RAAN*sin_w), (sin_RAAN*sin_incl)),
                          ((cos_RAAN*cos_incl*sin_w + sin_RAAN*cos_w), (cos_RAAN*cos_incl*cos_w - sin_RAAN*sin_w), (-cos_RAAN*sin_incl)),
                          ((sin_incl*sin_w), (sin_incl*cos_w), (cos_incl))))
        r0 = self.p/(1 + self.e*cos_theta)
        r_vect0 = array([[r0*cos_theta], [r0*sin_theta], [0]])
        r_vectf = Rot_Mat*r_vect0
        self.x = r_vectf.item(0)
        self.y = r_vectf.item(1)
        self.z = r_vectf.item(2)
        self.r = sqrt(self.x**2 + self.y**2 + self.z**2)
        decl_s = arcsin(self.z/self.r)
        RAAN_s = arctan2(self.y, self.x)

        lat = arctan(self.Er_Pr2*tan(decl_s)) % pi
        lng = (RAAN_s - self.GST0 - self.P_w*tnow)  % twopi

        self.lat = ((lat > pi/2)*(lat - pi) + (lat <= pi/2)*lat)*rad2deg
        self.lng = ((lng > pi)*(lng - twopi) + (lng <= pi)*lng)*rad2deg
        self.alt = self.r - self.getPlanetRadius()
        self.vt = self.h/self.r
        self.vr = self.mu*self.e*sin(self.theta)/self.h
        self.v = sqrt(self.vt**2 + self.vr**2)

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

