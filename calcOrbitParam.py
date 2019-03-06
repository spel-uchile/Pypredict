from numpy import abs, arccos, arctan, argmin, cos, cross, dot, pi, sin, sqrt, tan, zeros 

class Calc():

    def __init__(self, r1=None, dt=None, r2=None, v=None, Mt=5.9722*10**24):
        self.G = 6.67408*10**(-11)
        if (r1 is not None):
            self.newCalc(r1, dt=dt, r2=r2, v=v, Mt=Mt)
    
    def __call__(self):
        return self
    
    def newCalc(self, r1, dt=None, r2=None, v=None, Mt=5.9722*10**24):
        self.r1 = r1
        self.mu = self.G*Mt
        if (dt is not None and r2 is not None):
            self.dt = dt
            self.r2 = r2
            self.i = 0 # Default value
            self.getVelocities()
            self.getOrbitalParameters(self.r1, self.v1)
            if (self.i >= pi/2):
                self.getVelocities()
                self.getOrbitalParameters(self.r1, self.v1)
        elif (v is not None):
            self.v1 = v
            self.getOrbitalParameters(self.r1, self.v1)
        else:
            print("Two positions and time or a position and a velocity are needed.")
        
    def getdtheta(self):
        r1xr2 = cross(self.r1, self.r2)
        r1dotr2 = dot(self.r1, self.r2)
        if (r1xr2[2] >= 0):
            dtheta = arccos(r1dotr2/(self.r1_scalar*self.r2_scalar))
        elif (r1xr2[2] < 0):
            dtheta = 2*pi - arccos(r1dotr2/(self.r1_scalar*self.r2_scalar))
        if (self.i >= pi/2 or (self.dt < 2000 and dtheta > pi)):
            dtheta = 2*pi - dtheta
        return dtheta
        
    def getA(self, dtheta):
        A = sin(dtheta)*sqrt(self.r1_scalar*self.r2_scalar/(1 - cos(dtheta)))
        return A
        
    def gety(self, A, C, S, z):
        r1_km = self.r1_scalar/1000
        r2_km = self.r1_scalar/1000
        A_km = A/1000
        y = r1_km + r2_km + A_km*(z*S - 1)/sqrt(C)
        return y*1000
        
    def getyp(self, A, c):
        yp = A/4*sqrt(C)
        return yp
        
    def getf(self, y):   
        f = 1 - y/self.r1_scalar
        return f
        
    def getg(self, A, y):
        g = A*sqrt(y/self.mu)
        return g
        
    def getgp(self, y):
        gp = 1 - y/self.r2_scalar
        return gp
        
    def getF(self, A, C, S, y):
        F = (y/C)**(3/2)*S + A*sqrt(y) - sqrt(self.mu)*self.dt
        return F
        
    def getFp(self, A, C, S, y, z):
        if (z == 0):
            Fp = sqrt(2)/40*y**(3/2) + A/8*(sqrt(y) + A*sqrt(1/(2*y)))
        else:
            Fp = (y/C)**(3/2)*(1/(2*z)*(C - 3*S/(2*C)) + 3/4*S**2/C) + A/8*(3*S/C*sqrt(y) + A*sqrt(C/y))
        return Fp
        
    def getz(self, z0, F, Fp):
        z = z0 - F/Fp
        return z
        
    def getzNewton(self, A, z0):
        z = z0
        z0 = 0
        while (abs(z0 - z) > 0.00001):
            z0 = z
            C = 1/2 - z0/24
            S = 1/6 - z0/120
            y = self.gety(A, C, S, z0)
            F = self.getF(A, C, S, y)
            Fp = self.getFp(A, C, S, y, z0)
            z = self.getz(z0, F, Fp)
        return z
        
    def getz0(self, A):
        F_array = zeros(3000)
        z_array = zeros(3000)
        for i in range(0, 3000, 1):
            z_array[i] = i/1000
            C = 1/2 - z_array[i]/24
            S = 1/6 - z_array[i]/120
            y = self.gety(A, C, S, z_array[i])
            F_array[i] = self.getF(A, C, S, y)
        z0 = z_array[argmin(abs(F_array))]
        return z0
        
    def getMagnitude(self, vect):
        return sqrt(dot(vect, vect))
        
    def getVelocities(self):
        self.r1_scalar = self.getMagnitude(self.r1)
        self.r2_scalar = self.getMagnitude(self.r2)
        dtheta = self.getdtheta()
        A = self.getA(dtheta)
        z0 = self.getz0(A)
        z = self.getzNewton(A, z0)
        C = 1/2 - z/24
        S = 1/6 - z/120
        y = self.gety(A, C, S, z)
        f = self.getf(y)
        g = self.getg(A, y)
        gp = self.getgp(y)
        self.v1 = [0, 0, 0]
        self.v2 = [0, 0, 0]
        for i in range(3):
            self.v1[i] = 1/g*(self.r2[i] - f*self.r1[i])
            self.v2[i] = 1/g*(gp*self.r2[i] - self.r1[i])
        
    def getRadialVel(self, r, r_scalar, v):
        self.vr = dot(r, v)/r_scalar
        
    def getInclination(self, h):
        self.i = arccos(h[2]/self.h_scalar)
        
    def getRAAN(self, N, N_scalar):
        if (N[1] >= 0):
            self.RAAN = arccos(N[0]/N_scalar) % (2*pi)
        else:
            self.RAAN = (2*pi - arccos(N[0]/N_scalar)) % (2*pi)
        
    def getExcentricity(self, r, v):
        r_scalar = self.getMagnitude(r)
        v_scalar = self.getMagnitude(v)
        arg1 = dot(v, v) - self.mu/r_scalar
        arg2 = r_scalar*self.vr
        self.e = [0, 0, 0]
        for i in range(3):
            self.e[i] = (arg1*r[i] - arg2*v[i])/self.mu
        
    def getArgOfPerigee(self, N, N_scalar):
        arg = dot(N, self.e)/(N_scalar*self.e_scalar)
        if (self.e[2] >= 0):
            self.w = arccos(arg)
        else:
            self.w = 2*pi - arccos(arg)
        self.w = self.w % (2*pi)
        
    def getTrueAnomaly(self, r, r_scalar):
        arg = 1/(self.e_scalar*r_scalar)*dot(self.e, r)
        if (self.vr >= 0):
            self.theta = arccos(arg)
        else:
            self.theta = 2*pi - arccos(arg)

    def theta2E(self):
        return 2*arctan(sqrt((1 - self.e_scalar)/(1 + self.e_scalar))*tan(self.theta/2))

    def getMeanAnomaly(self):
        E = self.theta2E()
        self.MA = (E - self.e_scalar*sin(E)) % (2*pi)
        
    def getSemimajorAxis(self):
        rp = self.h_scalar**2/self.mu*(1/(1 + self.e_scalar*cos(0)))
        ra = self.h_scalar**2/self.mu*(1/(1 + self.e_scalar*cos(pi)))
        self.a = 0.5*(rp + ra)

    def getMeanVelocity(self):
        self.n = sqrt(self.mu/self.a**3)

    def getPeriod(self):
        self.T = 2*pi/self.n

    def getMeanMotion(self):
        self.MM = self.n*12*3600/pi
        
    def getOrbitalParameters(self, r, v):
        r_scalar = self.getMagnitude(r)
        h = cross(r, v)
        self.h_scalar = self.getMagnitude(h/1000000)*1000000
        self.getRadialVel(r, r_scalar, v)
        self.getExcentricity(r, v)
        self.e_scalar = self.getMagnitude(self.e)
        self.getSemimajorAxis()
        N = cross([0, 0, 1], h)
        N_scalar = self.getMagnitude(N/1000000)*1000000
        self.getRAAN(N, N_scalar)
        self.getInclination(h)
        self.getArgOfPerigee(N, N_scalar)
        self.getTrueAnomaly(r, r_scalar)
        self.getMeanAnomaly()
        self.getMeanVelocity()
        self.getMeanMotion()
        self.getPeriod()

    def getTLE(self):
        self.TLE1 = "{} {}{} {}{}{} {}{} {} {} {} {} {}{}".format("1",
                                                                  "     ",
                                                                  "U",
                                                                  "  ",
                                                                  "   ",
                                                                  "   ",
                                                                  "  ",
                                                                  "            ",
                                                                  "          ",
                                                                  "        ",
                                                                  "        ",
                                                                  "        ",
                                                                  "0",
                                                                  "    ",
                                                                  " ")
        self.TLE2 = "{} {} {:08.4f} {:08.4f} {:07d} {:08.4f} {:08.4f} {:011d}{}{}".format("2",
                                                         "     ",
                                                         self.i*180/pi,
                                                         self.RAAN*180/pi,
                                                         self.e*10000000,
                                                         self.w*180/pi,
                                                         self.MA*180/pi,
                                                         self.MM,
                                                         "     ",
                                                         " ")
        return "{}\n{}".format(self.TLE1, self.TLE2)
