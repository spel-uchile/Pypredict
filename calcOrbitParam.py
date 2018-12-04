import numpy as np

class Calc():

    def __init__(self, dt, r1, r2, Mt=5.9722*10**24):
        self.newCalc(dt, r1, r2)
        self.getVelocities()
        self.getOrbitalParameters(self.r1, self.v1)
    
    def __call__(self):
        return self
    
    def newCalc(self, dt, r1, r2, Mt=5.9722*10**24):
        G = 6.67408*10**(-11)
        self.dt = dt
        self.r1 = r1
        self.r2 = r2
        self.mu = G*Mt
        
    def getdtheta(self):
        r1xr2 = np.cross(self.r1, self.r2)
        if (r1xr2[0, 2] >= 0):
            dtheta = np.arccos(np.dot(self.r1, self.r2.transpose())/(self.r1_scalar*self.r2_scalar))
        elif (r1xr2[0, 2] < 0):
            dtheta = 2*np.pi - np.arccos(np.dot(self.r1, self.r2.transpose())/(self.r1_scalar*self.r2_scalar))
        return np.asscalar(dtheta)
        
    def getA(self, dtheta):
        A = np.sin(dtheta)*np.sqrt(self.r1_scalar*self.r2_scalar/(1 - np.cos(dtheta)))
        return A
        
    def gety(self, A, C, S, z):
        y = self.r1_scalar + self.r2_scalar + A*(z*S - 1)/np.sqrt(C)
        return y
        
    def getyp(self, A, c):
        yp = A/4*np.sqrt(C)
        return yp
        
    def getf(self, y):   
        f = 1 - y/self.r1_scalar
        return f
        
    def getg(self, A, y):
        g = A*np.sqrt(y/self.mu)
        return g
        
    def getgp(self, y):
        gp = 1 - y/self.r2_scalar
        return gp
        
    def getF(self, A, C, S, y):
        F = (y/C)**(3/2)*S + A*np.sqrt(y) - np.sqrt(self.mu)*self.dt
        return F
        
    def getFp(self, A, C, S, y, z):
        if (z == 0):
            Fp = np.sqrt(2)/40*y**(3/2) + A/8*(np.sqrt(y) + A*np.sqrt(1/(2*y)))
        else:
            Fp = (y/C)**(3/2)*(1/(2*z)*(C - 3*S/(2*C)) + 3/4*S**2/C) + A/8*(3*S/C*np.sqrt(y) + A*np.sqrt(C/y))
        return Fp
        
    def getz(self, z0, F, Fp):
        z = z0 - F/Fp
        return z
        
    def getzNewton(self, A, z0):
        z = z0
        z0 = 0
        while (np.abs(z0 - z) > 0.00001):
            z0 = z
            C = 1/2 - z0/24
            S = 1/6 - z0/120
            y = self.gety(A, C, S, z0)
            F = self.getF(A, C, S, y)
            Fp = self.getFp(A, C, S, y, z0)
            z = self.getz(z0, F, Fp)
        return z
        
    def getz0(self, A):
        F_array = np.zeros(30)
        z_array = np.zeros(30)
        for i in range(0, 30, 1):
            z_array[i] = i/10
            C = 1/2 - z_array[i]/24
            S = 1/6 - z_array[i]/120
            y = self.gety(A, C, S, z_array[i])
            F_array[i] = self.getF(A, C, S, y)
        z0 = z_array[np.argmin(np.abs(F_array))]
        return z0
        
    def getMagnitude(self, vect):
        scalar = np.asscalar(np.sqrt(np.dot(vect, vect.transpose())))
        return scalar
        
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
        self.v1 = 1/g*(self.r2 - f*self.r1)
        self.v2 = 1/g*(gp*self.r2 - self.r1)
        
    def getRadialVel(self, r, r_scalar, v):
        vr = np.asscalar(np.dot(r, v.transpose()))/r_scalar
        return vr
        
    def getInclination(self, h):
        self.i = np.arccos(h[0, 2]/self.h_scalar)
        
    def getRAN(self, N, N_scalar):
        if (N[0, 1] >= 0):
            self.RAN = np.arccos(N[0, 0]/N_scalar)
        else:
            self.RAN = 2*np.pi - np.arccos(N[0, 0]/N_scalar)
        
    def getExcentricity(self, r, v, vr):
        r_scalar = self.getMagnitude(r)
        v_scalar = self.getMagnitude(v)
        e = 1/self.mu*((v_scalar**2 - self.mu/r_scalar)*r - r_scalar*vr*v)
        return e
        
    def getArgOfPerigee(self, N, N_scalar, e):
        arg = 1/(N_scalar*self.e_scalar)*np.asscalar(np.dot(N, e.transpose()))
        if (e[0,2] >= 0):
            self.w = np.arccos(arg)
        else:
            self.w = 2*np.pi - np.arccos(arg)
        
    def getTrueAnomaly(self, e, r, r_scalar, vr):
        arg = 1/(self.e_scalar*r_scalar)*np.asscalar(np.dot(e, r.transpose()))
        if (vr >= 0):
            self.theta = np.arccos(arg)
        else:
            self.theta = 2*np.pi - np.arccos(arg)
        
    def getSemimajorAxis(self):
        rp = self.h_scalar**2/self.mu*(1/(1 + self.e_scalar*np.cos(0)))
        ra = self.h_scalar**2/self.mu*(1/(1 + self.e_scalar*np.cos(np.pi)))
        self.a = 0.5*(rp + ra)
        
    def getPeriod(self):
        self.T = 2*np.pi/np.sqrt(self.mu)*self.a**(3/2)
        
    def getOrbitalParameters(self, r, v):
        r_scalar = self.getMagnitude(r)
        h = np.cross(r, v)
        self.h_scalar = self.getMagnitude(h)
        vr = self.getRadialVel(r, r_scalar, v)
        e = self.getExcentricity(r, v, vr)
        self.e_scalar = self.getMagnitude(e)
        self.getSemimajorAxis()
        N = np.cross([[0, 0, 1]], h)
        N_scalar = self.getMagnitude(N)
        self.getRAN(N, N_scalar)
        self.getInclination(h)
        self.getArgOfPerigee(N, N_scalar, e)
        self.getTrueAnomaly(e, r, r_scalar, vr)
        self.getPeriod()
