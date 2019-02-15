from numpy import matrix, sqrt
from sat import Sat

class TDOA(object):
    __slots__ = ["c"]
    def __init__(self, c=299792458):
        self.c = c

    def __call__(self):
        return self

    def calculateLocation(self, Sat0, Sat1, Sat2, Sat3, Sat4=None):
        x1 = Sat1.x
        y1 = Sat1.y
        x21 = Sat2.x - x1
        y21 = Sat2.y - y1
        x31 = Sat3.x - x1
        y31 = Sat3.y - y1
        r0 = Sat0.getXYZ()
        r1 = Sat1.getXYZ()
        r2 = Sat2.getXYZ()
        r3 = Sat3.getXYZ()
        d21 = self.c*self.getdt(r0, r1, r2)
        d31 = self.c*self.getdt(r0, r1, r3)
        if (Sat4 is None):
            print("This result is not valid for 3D coordinates, one more satellite is needed")
            C1 = x1**2 + y1**2
            C2 = r2[0]**2 + r2[1]**2
            C3 = r3[0]**2 + r3[1]**2
            a = ((y21*d31 - y31*d21)**2 + (x31*d21 - x21*d31)**2 - (x31*y21 - x21*y31)**2)
            b = 2*(y21*d31 - y31*d21)*(x1*(x31*y21 - x21*y31)**2 + y21/2*(d31**2 - C3 + C1) - y31/2*(d21**2 - C2 + C1))
            b = b + 2*(x31*d21 - x21*d31)*(y1*(x31*y21 - x21*y31)**2 + x31/2*(d21**2 - C2 + C1) - x21/2*(d31**2 - C3 + C1))
            c = (x31/2*(d21**2 - C2 + C1) - x21/2*(d31**2 - C3 + C1))**2 + (y21/2*(d31**2 - C3 + C1) - y31/2*(d21**2 - C2 + C1))**2
            c = c + C1*(x31*y21 - x21*y31)**4
            c = c + 2*x1*(x31*y21 - x21*y31)**2*(y21/2*(d31**2 - C3 + C1) - y31/2*(d21**2 - C2 + C1))
            c = c + 2*y1*(x31*y21 - x21*y31)**2*(x31/2*(d21**2 - C2 + C1) - x21/2*(d31**2 - C3 + C1))
            print("a: " + str(a))
            print("b: " + str(b))
            print("c: " + str(c))
            d1 = (-b + sqrt(b**2 - 4*a*c))/(2*a)
            print("d1 real: " + str(self.getDistance(r0, r1)) + "    d1 calculated: " + str(d1))
            #M = matrix(((x21, y21),
            #               (x31, y31)))
            #H = self.getH(r0, r1, r2, r3)
            #x, y = -M.I*(matrix(((d21), (d31))).transpose()*d1 + H)
            x = -y21/(x31*y21 - x21*y31)*(d31*d1 + 1/2*(d31**2 - C3 + C1))
            x = x + y31/(x31*y21 - x21*y31)*(d21*d1 + 1/2*(d21**2 - C2 + C1))
            y = -x31/(x31*y21 - x21*y31)*(d21*d1 + 1/2*(d21**2 - C2 + C1))
            y = y + x21/(x31*y21 - x21*y31)*(d31*d1 + 1/2*(d31**2 - C3 + C1))
            z = r1[2]
        else:
            print("This result is valid for 3D coordinates")
            z1 = Sat1.z
            C1 = self.getC(Sat1.getXYZ())
            C2 = self.getC(Sat2.getXYZ())
            C3 = self.getC(Sat3.getXYZ())
            C4 = self.getC(Sat4.getXYZ())
            z21 = Sat2.z - z1
            z31 = Sat3.z - z1
            x41 = Sat4.x - x1
            y41 = Sat4.y - y1
            z41 = Sat4.z - z1
            r4 = Sat4.getXYZ()
            d41 = self.c*self.getdt(r0, r1, r4)
            D = (x31*y41 - x41*y31)*z21 + (x41*y21 - x21*y41)*z31 + (x21*y31 - x31*y21)*z41
            a1 = (y31*z41 - y41*z31)*d21 - (y21*z41 - y41*z21)*d31 - (y31*z21 - y21*z31)*d41
            a2 = (x31*z41 - x41*x31)*d21 - (x21*z41 - x41*z21)*d31 + (x21*z31 - x31*z21)*d41
            a3 = (x31*y41 - x41*y31)*d21 - (x21*y41 - x41*y21)*d31 - (x31*y21 - x21*y31)*d41
            a = a1**2 + a2**2 + a3**2 - D
            b = a1*(2*x1*D + (y31*z41 - y41*z31)*(d21**2 - C2 + C1) - (y21*z41 - y41*z21)*(d31**2 - C3 + C1) - (y31*z21 - y21*z31)*(d41**2 - C4 + C1))
            b += a2*(2*y1*D - (x31*z41 - x41*z31)*(d21**2 - C2 + C1) + (x21*z41 - x41*z21)*(d31**2 - C3 + C1) - (x21*z31 - x31*z21)*(d41**2 - C4 + C1))
            b += a3*(2*z1*D + (x31*y41 - x41*y31)*(d21**2 - C2 + C1) - (x21*y41 - x41*y21)*(d31**2 - C3 + C1) - (x31*y21 - x21*y31)*(d41**2 - C4 + C1))
            c = x1**2*D**2 + 1/4*(y31*z41 - y41*z31)**2*(d21**2 - C2 + C1)**2 + 1/4*(y21*z41 - y41*z21)**2*(d31**2 - C3 + C1)**2
            c += 1/4*(y31*z21 - y21*z31)**2*(d41**2 - C4 + C1)**2 + x1*D*(y31*z41 - y41*z31)*(d21**2 - C2 + C1)
            c += -x1*D*(y21*z41 - y41*z21)*(d31**2 - C3 + C1) - x1*D*(y31*z21 - y21*z31)*(d41**2 - C4 + C1)
            c += y1**2*D**2 + 1/4*(x31*z41 - x41*x31)**2*(d21**2 - C2 + C1)**2 + 1/4*(x21*z41 - x41*z21)**2*(d31**2 - C3 + C1)**2
            c += -1/4*(x21*z31 - x31*z21)**2*(d41**2 - C4 + C1)**2 - y1*D*(x31*z41 - x41*z31)*(d21**2 - C2 + C1)
            c += y1*D*(x21*z41 - x41*z21)*(d31**2 - C3 + C1) - y1*D*(x21*z31 - x31*z21)*(d41**2 - C4 + C1)
            c += z1**2*D**2 + 1/4*(x31*y41 - x41*y31)**2*(d21**2 - C2 + C1)**2 + 1/4*(x21*y41 - x41*y21)**2*(d31**2 - C3 + C1)**2
            c += 1/4*(x31*y21 - x21*y31)**2*(d41**2 - C4 + C1)**2 + z1*D*(x31*y41 - x41*y31)*(d21**2 - C2 + C1)
            c += -z1*D*(x21*y41 - x41*y21)*(d31**2 - C3 + C1) - z1*D*(x31*y21 - x21*y31)*(d41**2 - C4 + C1)
            c += -1/2*(y21*z41 - y41*z21)*(y31*z41 - y41*z31)*(d21**2 - C2 + C1)*(d31**2 - C3 + C1)
            c += -1/2*(y31*z21 - y21*z31)*(y31*z41 - y41*z31)*(d21**2 - C2 + C1)*(d41**2 - C4 + C1)
            c += 1/2*(y21*z41 - y41*z21)*(y31*z21 - y21*z31)*(d31**2 - C3 + C1)*(d41**2 - C4 + C1)
            c += -1/2*(x21*z41 - x41*z21)*(x31*z41 - x41*z31)*(d21**2 - C2 + C1)*(d31**2 - C3 + C1)
            c += 1/2*(x21*z31 - x31*z21)*(x31*z41 - x41*z31)*(d21**2 - C2 + C1)*(d41**2 - C4 + C1)
            c += -1/2*(x21*z31 - x31*z21)*(x21*z41 - x41*z21)*(d31**2 - C3 + C1)*(d41**2 - C4 + C1)
            c += -1/2*(x21*y41 - x41*y21)*(x31*y41 - x41*y31)*(d21**2 - C2 +C1)*(d31**2 - C3 + C1)
            c += -1/2*(x31*y21 - x21*y31)*(x31*y41 - x41*y31)*(d21**2 - C2 + C1)*(d41**2 - C4 + C1)
            c += 1/2*(x21*y41 - x41*y21)*(x31*y21 - x21*y31)*(d31**2 - C3 + C1)*(d41**2 - C4 + C1)
            print("a: " + str(a))
            print("b: " + str(b))
            print("c: " + str(c))
            d1 = (-b + sqrt(b**2 - 4*a*c))/(2*a)
            M = matrix(((x21, y21, z21),
                           (x31, y31, z31),
                           (x41, y41, z41)))
            H = self.getH(r0, r1, r2, r3, r4)
            x, y, z = -M.I*(matrix(((d21), (d31), (d41))).transpose()*d1 + H)
        return x, y, z

    def getdt(self, r0, r1, r2):
        d1 = self.getDistance(r0, r1)
        d2 = self.getDistance(r0, r2)
        return (d2 - d1)/self.c

    def getDistance(self, ra, rb):
        xa = ra[0]
        ya = ra[1]
        za = ra[2]
        xb = rb[0]
        yb = rb[1]
        zb = rb[2]
        return sqrt((xa - xb)**2 + (ya - yb)**2 + (za -  zb)**2)

    def getC(self, r):
        return sqrt(r[0]**2 + r[1]**2 + r[2]**2)

    def getH(self, r0, r1, r2, r3, r4=None):
        d21 = self.c*self.getdt(r0, r1, r2)
        d31 = self.c*self.getdt(r0, r1, r3)
        C1 = self.getC(r1)
        C2 = self.getC(r2)
        C3 = self.getC(r3)
        if (r4 is None):
            H = 1/2*matrix(((d21**2 - (C2 - C1)),
                               (d31**2 - (C3 - C1))))
        else:
            d41 = self.c*self.getdt(r0, r1, r4)
            C4 = self.getC(r4)
            H = 1/2*matrix(((d21**2 - (C2 - C1)),
                               (d31**2 - (C3 - C1)),
                               (d41**2 - (C4 - C1))))
        return H.transpose()

    def getG(self, r0, r1, r2, r3):
        d21 = self.c*self.getdt(r0, r1, r2)
        d31 = self.c*self.getdt(r0, r1, r3)
        x21 = r2[0] - r1[0]
        x31 = r3[0] - r1[0]
        y21 = r2[1] - r1[1]
        y31 = r3[1] - r1[1]
        G = matrix(((x21, y21, d21),
                       (x31, y31, d31)))
        return G

    def get_d1(self):
        a = ((y21*d31 - y31*d21)**2 + (x31*d21 - x21*d31)**2 - (x31*y21 - x21y31)**2)
        b = 2*(y21*d31 - y31*d21)*(x1*(x31*y21 - x21*y31)**2 + y21/2*(d31**2 - C3 + C1) - y31/2*(d21**2 - C2 + C1))
        b = b + 2*(x31*d21 - x21*d31)*(y1*(x31*y21 - x21*y31) + x31/2*(d21**2 - C2 + C1) - x21/2*(d31**2 - C2 + C1))
        c = (x31/2*(d21**2 - C2 + C1) - x21/2*(d31**2 - C3 + C1))**2 + (y21/2*(d31**2 - C3 + C1) - y31/2*(d21**2 - C2 + C1))**2
        c = c + C1*(x31*y21 - x21*y31)**4
        c = c + 2*x1*(x31*y21 - x21*y31)**2*(y21/2*(d31**2 - C3 + C1) - y31/2*(d21**2 - C2 + C1))
        c = c + 2*y1*(x31*y21 - x21*y31)**2*(x31/2*(d21**2 - C2 + C1) - x21/2*(d31**2 - C3 + C1))
        return (-b + sqrt(b**2 - 4*a*c))/(2*a)
