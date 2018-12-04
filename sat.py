from numpy import abs, arccos, arcsin, arctan, arctan2, array, cos, matrix, pi, sin, sqrt, tan
from node import Node
from datetime import datetime
from pyorbital import tlefile

class Sat(Node):

    def __init__(self, name="", lat=0, lng=0, alt=0, freq=437225000, tle=None):
        deg2rad = pi/180
        self.name = name
        self.lat = lat
        self.lng = lng
        self.alt = alt
        self.freq = freq
        dayinsec = 24*3600                                        # Day in seconds
        if (tle is not None):
            self.incl = tle.inclination*deg2rad
            self.RAN0 = tle.right_ascension*deg2rad
            self.e = tle.excentricity
            self.w0 = tle.arg_perigee*deg2rad
            self.MA0 = tle.mean_anomaly*deg2rad
            rpd2radps = (2*pi/dayinsec)                           # Revolutions per day to radians per second
            self.n = tle.mean_motion*rpd2radps                    # Radians per second
            self.epoch_year = tle.epoch_year
            self.epoch_day = tle.epoch_day
            print(self.name + " TLE found!")
        else:
            self.incl = 97.39*deg2rad
            self.RAN0 = 281.63*deg2rad
            self.e = 0.001
            self.w0 = 248*deg2rad
            self.MA0 = 112*deg2rad
            self.n = 0.0011
            self.epoch_year = 18
            self.epoch_day = 222.64
            print("No TLE found, used default parameters instead")
        self.RAN = self.RAN0
        self.w = self.w0
        self.MA = self.MA0
        self.theta = self.MA
        G = 6.67408*10**(-11)                                     # Gravitational constant
        Mt = 5.9722*10**24                                        # Earth mass
        ut = G*Mt                                                 # Gravitational parameter for Earth
        self.a = (ut/self.n**2)**(1/3)                            # Semi-major axis
        self.p = self.a*(1 - self.e**2)                           # Semilatus-rectum
        self.t0 = (self.epoch_day - int(self.epoch_day))*dayinsec # Initial time in seconds
        D, Month, Y = self.getDMY()                               # Get day, month and year from TLE data
        self.GST0 = self.getGST(int(D), Month, Y)                 # Get Greenwich sideral time
        self.updateOrbitalParameters()                            # Update orbital parameters

    def __call__(self):
        return self

    def setInclination(self, incl):
        self.incl = incl

    def setRAN(self, RAN):
        self.RAN = RAN

    def setEccentricity(self, e):
        self.e = e

    def setMeanAnomaly(self, MA):
        self.MA = MA

    def getInclination(self):
        return self.incl

    def getRAN(self):
        return self.RAN

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

    def getXYZ(self):
        return self.x, self.y, self.z

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

    def getEarthRadius(self):
        Eq_r = 6378000   # Equatorial radius
        Po_r = 6356000   # Polar radius
        cos_lat = cos(self.lat)
        sin_lat = sin(self.lat)
        Earth_r = sqrt(((Eq_r**2*cos_lat)**2 + (Po_r**2*sin_lat)**2)/((Eq_r*cos_lat)**2 + (Po_r*sin_lat)**2))
        return Earth_r

    def getPeriod(self):
        return 2*pi/self.n

    def getCoverage(self):
        E_r = self.getEarthRadius()
        ang = arccos(E_r/(E_r + self.alt))*180/pi
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

    def getLocation(self, tf, dt):
        tnow = self.getTnow()
        lat = []
        lng = []
        J2 = 1.083*10**(-3)                # Second degree harmonic model of Earth
        rad2deg = 180/pi                   # Radian to degrees
        re_rp2 = (6378/6356)**2            # (equatorial radius/polar radius)^2
        twopi = 2*pi                       # Two times pi
        E_r = 6371000                      # Earth radius
        wt = twopi/(23*3600 + 56*60 + 4.1) # Earth's angular velocity
        aux = self.n*(E_r**2/self.p**2)*J2 # Calculate only one time
        cos_incl = cos(self.incl)
        sin_incl = sin(self.incl)
        for i in range(0, tf, dt):
            t = tnow + i
            j = int(i/dt)

            RAN = self.RAN0 - 1.5*aux*cos_incl*(t - self.t0)
            w = self.w0 + 0.75*aux*(5*cos_incl**2 - 1)*(t - self.t0)
            MA = self.MA0 + (self.n + 0.75*aux*(sqrt(1 - self.e**2))*(2 - 3*sin_incl**2))*(t - self.t0)
            E = self.M2E(MA)
            theta = self.E2theta(E)

            cos_RAN = cos(RAN)
            sin_RAN = sin(RAN)
            cos_w = cos(w)
            sin_w = sin(w)
            cos_theta = cos(theta)
            sin_theta = sin(theta)

            Rot_Mat = matrix((((cos_RAN*cos_w - sin_RAN*cos_incl*sin_w), (-sin_RAN*cos_incl*cos_w - cos_RAN*sin_w), (sin_RAN*sin_incl)),
                              ((cos_RAN*cos_incl*sin_w + sin_RAN*cos_w), (cos_RAN*cos_incl*cos_w - sin_RAN*sin_w), (-cos_RAN*sin_incl)),
                              ((sin_incl*sin_w), (sin_incl*cos_w), (cos_incl))))
            r0 = self.p/(1 + self.e*cos_theta)
            r_vect0 = array([[r0*cos_theta], [r0*sin_theta], [0]])
            r_vectf = Rot_Mat*r_vect0
            r_s = sqrt(r_vectf.item(0)**2 + r_vectf.item(1)**2 + r_vectf.item(2)**2)
            decl_s = arcsin(r_vectf.item(2)/r_s)
            RAN_s = arctan2(r_vectf.item(1),r_vectf.item(0))

            lat.append(arctan(re_rp2*tan(decl_s)) % pi)
            lng.append((RAN_s - self.GST0 - wt*t)  % twopi)

            lat[j] = ((lat[j] > pi/2)*(lat[j] - pi) + (lat[j] <= pi/2)*(lat[j]))*rad2deg
            lng[j] = ((lng[j] > pi)*(lng[j] - twopi) + (lng[j] <= pi)*(lng[j]))*rad2deg
        return lat, lng

    def updateOrbitalParameters(self):
        tnow = self.getTnow()
        J2 = 1.083*10**(-3)                             # Second degree harmonic model of Earth
        rad2deg = 180/pi                                # Radian to degrees
        re_rp2 = (6378/6356)**2                         # (equatorial radius/polar radius)^2
        twopi = 2*pi                                    # Two times pi
        E_r = 6371000                                   # Earth radius
        wt = twopi/(23*3600 + 56*60 + 4.1)              # Earth's angular velocity
        aux = self.n*(E_r**2/self.p**2)*J2              # Calculate only one time
        cos_incl = cos(self.incl)
        sin_incl = sin(self.incl)

        self.RAN = self.RAN0 - 1.5*aux*cos_incl*(tnow - self.t0)
        self.w = self.w0 + 0.75*aux*(5*cos_incl**2 - 1)*(tnow - self.t0)
        self.MA = self.MA0 + (self.n + 0.75*aux*(sqrt(1 - self.e**2))*(2 - 3*sin_incl**2))*(tnow - self.t0)
        E = self.M2E(self.MA)
        self.theta = self.E2theta(E)
        cos_RAN = cos(self.RAN)
        sin_RAN = sin(self.RAN)
        cos_w = cos(self.w)
        sin_w = sin(self.w)
        cos_theta = cos(self.theta)
        sin_theta = sin(self.theta)

        Rot_Mat = matrix((((cos_RAN*cos_w - sin_RAN*cos_incl*sin_w), (-sin_RAN*cos_incl*cos_w - cos_RAN*sin_w), (sin_RAN*sin_incl)),
                          ((cos_RAN*cos_incl*sin_w + sin_RAN*cos_w), (cos_RAN*cos_incl*cos_w - sin_RAN*sin_w), (-cos_RAN*sin_incl)),
                          ((sin_incl*sin_w), (sin_incl*cos_w), (cos_incl))))
        r0 = self.p/(1 + self.e*cos_theta)
        r_vect0 = array([[r0*cos_theta], [r0*sin_theta], [0]])
        r_vectf = Rot_Mat*r_vect0
        self.x = r_vectf.item(0)
        self.y = r_vectf.item(1)
        self.z = r_vectf.item(2)
        r_s = sqrt(self.x**2 + self.y**2 + self.z**2)
        decl_s = arcsin(self.z/r_s)
        RAN_s = arctan2(self.y, self.x)

        lat = arctan(re_rp2*tan(decl_s)) % pi
        lng = (RAN_s - self.GST0 - wt*tnow)  % twopi

        self.lat = ((lat > pi/2)*(lat - pi) + (lat <= pi/2)*lat)*rad2deg
        self.lng = ((lng > pi)*(lng - twopi) + (lng <= pi)*lng)*rad2deg
        self.alt = r_s - self.getEarthRadius()
