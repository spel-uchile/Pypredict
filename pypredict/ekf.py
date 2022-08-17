"""
                                Pypredict
    Orbit prediction software. Displays the satellites' position and
    orbital parameters in real time. Simulates satellite localization
    and deployment.
    
    Copyright (C) 2018-2022, Matías Vidal Valladares, matvidal.
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
from numpy import linalg, matrix, pi, sqrt, transpose
from numpy.linalg import inv

class EKF(object):
    __slots__ = ["A", "H", "H_tr", "I", "K",
                 "P_pred", "P_up", "R", "Q",
                 "S", "x_pred", "x_up", "y"]
    def __init__(self):
        self.initMatrices()

    def __call__(self):
        return self

    def __str__(self):
        return "EKF applied to estimate a satellite's position and velocity"

    def newIteration(self, dt, n, rc, x, s1, s2, z):
        """
        Runs one iteration of the Kalman filter's algorithm.

        Parameters
        ----------
        dt : float
             Time step since the last iteration
        n  : float
             Mean motion of the main satellite
        rc : float
             Current orbit radius of the chief satellite
        x  : matrix
             State vector with the last estimated position and velocity
        s1 : matrix
             Position vector of the s1 satellite
        s2 : matrix
             Position vector of the s2 satellite
        z  : matrix
             Measurement of the satellite's new position
        """
        self.Predict(dt, n, rc, x)
        self.Update(s1, s2, z)

    def initMatrices(self):
        """Initialize all the matrices needed to run the kalman filter."""
        self.initCovarianceMatrix()
        self.initH()
        self.initIdentityMatrix()
        self.initR()
        self.initQ()
        self.initPropagationMatrix()

    def HCW(self, x, n, rc):
        #G = 6.67408*10**(-11)
        #M = 5.9722*10**24
        mu = 398589405759999.94 #G*M
        mu_div_norm_rd3 = mu/linalg.norm(matrix([[0], [0], [-rc]]) + x[0:3,0])**3
        n2_minus_mu_div_norm_rd3 = n**2 - mu_div_norm_rd3
        two_n = 2*n
        #return matrix([[x[3,0]],
        #               [x[4,0]],
        #               [x[5,0]],
        #               [n2_minus_mu_div_norm_rd3*x[0,0] - two_n*x[5,0]],
        #               [-mu_div_norm_rd3*x[1,0]],
        #               [n2_minus_mu_div_norm_rd3*(x[2,0] - rc) + two_n*x[3,0]]])
        return matrix([[x[3,0]],
                       [x[4,0]],
                       [x[5,0]],
                       [n2_minus_mu_div_norm_rd3*x[0,0] + two_n*x[5,0]],
                       [-mu_div_norm_rd3*x[1,0]],
                       [n2_minus_mu_div_norm_rd3*(x[2,0] - rc) - two_n*x[3,0]]])

    def RK4(self, x, h, iterations, n, rc):
        half_h = h*0.5
        h_div_6 = h/6
        t = 0
        for k in range(iterations):
            #k1 = HCW(t, x, n, rc)
            #k2 = HCW(t + half_h, x + half_h*k1, n, rc)
            #k3 = HCW(t + half_h, x + half_h*k2, n, rc)
            #k4 = HCW(t + h, x + h*k3, n, rc)
            t_plus_half_h = t + half_h
            t_plus_h = t + h
            k1 = self.HCW(x, n, rc)*t
            k2 = self.HCW(x + half_h*k1, n, rc)*t_plus_half_h
            k3 = self.HCW(x + half_h*k2, n, rc)*t_plus_half_h
            k4 = self.HCW(x + h*k3, n, rc)*t_plus_h
            x = x + h_div_6*(k1 + 2*k2 + 2*k3 + k4)
            t = t_plus_h
        return x

    def initPropagationMatrix(self):
        """
        Initialize the propagation matrix A used to propagate the
        state vector.
        """
        self.A = matrix([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

    def initH(self, model=1):
        """
        Initialize the matrix that describes the observation model.
        """
        #self.h = matrix([[linalg.norm(u - s2) - linalg.norm(u - s1)],
        #                 [arctan2(u[1,0] - s1[1,0], u[0,0] - s1[0,0])],
        #                 [arctan2(u[2,0] - s1[2,0], sqrt((u[0,0] - s1[0,0])**2 + (u[1,0] - s1[1,0])**2))],
        #                 [arctan2(u[1,0] - s2[1,0], u[0,0] - s2[0,0])],
        #                 [arctan2(u[2,0] - s2[2,0], sqrt((u[0,0] - s2[0,0])**2 + (u[1,0] - s2[1,0])**2))]])
        if (model == 0):
            self.H = matrix([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
        else:
            self.H = matrix([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]])
        self.H_tr = transpose(self.H)

    def initCovarianceMatrix(self):
        """Initialize the kalman filter's covariance matrix."""
        self.P_up = matrix([[200.0**2, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 200.0**2, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 200.0**2, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 200.0**2, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 200.0**2, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 200.0**2]])

    def initIdentityMatrix(self):
        """Initialize the 6x6 identity matrix."""
        self.I = matrix([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])

    def initR(self, model=1):
        """Initialize the covariance matrix of the observation noise."""
        if (model == 0):
            std_RD = 10.0
            std_AOA = 1.0*pi/180.0
            self.R = matrix([[std_RD**2, 0.0, 0.0, 0.0, 0.0],
                             [0.0, std_AOA**2, 0.0, 0.0, 0.0],
                             [0.0, 0.0, std_AOA**2, 0.0, 0.0],
                             [0.0, 0.0, 0.0, std_AOA**2, 0.0],
                             [0.0, 0.0, 0.0, 0.0, std_AOA**2]])
        else:
            std_pos = 2#30/sqrt(3)
            self.R = matrix([[std_pos**2, 0.0, 0.0],
                             [0.0, std_pos**2, 0.0],
                             [0.0, 0.0, std_pos**2]])

    def initQ(self):
        """Initialize the covariance matrix of the process noise."""
        self.Q = matrix([[5.0**2, 0.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 5.0**2, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 5.0**2, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 2.0**2, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 2.0**2, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 2.0**2]])

    def updateH(self, x, s1, s2, model=1):
        """
        Updates the matrix that describes the observation model.

        x  :    matrix
                State vector
        s1 :    matrix
                Position vector of the s1 satellite
        s2 :    matrix
                Position vector of the s2 satellite
        model : int
                Selects the measurement model
        """
        if (model == 0):
            norm_u_s1 = linalg.norm(x[:3,0] - s1)
            norm_u_s2 = linalg.norm(x[:3,0] - s2)
            self.H[0,0] = (x[0,0] - s2[0,0])/norm_u_s2 - (x[0,0] - s1[0,0])/norm_u_s1
            self.H[0,1] = (x[1,0] - s2[1,0])/norm_u_s2 - (x[1,0] - s1[1,0])/norm_u_s1
            self.H[0,2] = (x[2,0] - s2[2,0])/norm_u_s2 - (x[2,0] - s1[2,0])/norm_u_s1
            self.H[1,0] = -(x[1,0] - s1[1,0])/((x[0,0] - s1_x[0,0])**2 + (x[1,0] - s1[1,0])**2)
            self.H[1,1] = (x[0,0] - s1[0,0])/((x[0,0] - s1[0,0])**2 + (x[1,0] - s1[1,0])**2)
            self.H[2,0] = -(x[0,0] - s1[0,0])*(x[2,0] - s1[2,0])/(sqrt((-s1[0,0] + x[0,0])**2 + (-s1[1,0] + x[1,0])**2)*norm_u_s1**2)
            self.H[2,1] = -(x[1,0] - s1[1,0])*(x[2,0] - s1[2,0])/(sqrt((-s1[0,0] + x[0,0])**2 + (-s1[1,0] + x[1,0])**2)*norm_u_s1**2)
            self.H[2,2] = sqrt((-s1[0,0] + x[0,0])**2 + (-s1[1,0] + x[1,0])**2)/norm_u_s1**2
            self.H[3,0] = -(x[1,0] - s2[1,0])/((x[0,0] - s2[0,0])**2 + (x[1,0] - s2[1,0])**2)
            self.H[3,1] = (x[0,0] - s2[0,0])/((x[0,0] - s2[0,0])**2 + (x[1,0] - s2[1,0])**2)       
            self.H[4,0] = -(x[0,0] - s2[0,0])*(x[2,0] - s2[2,0])/(sqrt((-s2[0,0] + x[0,0])**2 + (-s2[1,0] + x[1,0])**2)*norm_u_s2**2)
            self.H[4,1] = -(x[1,0] - s2[1,0])*(x[2,0] - s2[2,0])/(sqrt((-s2[0,0] + x[0,0])**2 + (-s2[1,0] + x[1,0])**2)*norm_u_s2**2)
            self.H[4,2] = sqrt((x[0,0] - s2[0,0])**2 + (x[1,0] - s2[1,0])**2)/norm_u_s2**2
        else:
            self.H[0,0] = x[0,0]
            self.H[1,1] = x[1,0]
            self.H[2,2] = x[2,0]
        self.H_tr = transpose(self.H)

    def updatePropagationMatrix(self, n, rc, x):
        """
        Initialize the propagation matrix A used to propagate the
        state vector.

        Parameters
        ----------
        n  : float
             Mean motion of the main satellite

        rc : float
             Distance between the chief satellite and the center of the Earth

        x  : matrix
             State vector
        """
        G = 6.67408*10**(-11)
        M = 5.9722*10**24
        mu = G*M
        z_minus_rc = x[2,0] - rc
        norm_rd = sqrt(x[0,0]**2 + x[1,0]**2 + z_minus_rc**2)
        mu_div_norm_rd3 = mu/norm_rd**3
        mu_div_norm_rd5 = mu/norm_rd**5
        xy = x[0,0]*x[1,0]
        self.A[3,0] = 3*x[0,0]**2*mu_div_norm_rd5 - mu_div_norm_rd3 + n**2
        self.A[3,1] = 3*xy*mu_div_norm_rd5
        self.A[3,2] = 3*x[0,0]*z_minus_rc*mu_div_norm_rd5
        self.A[3,5] = 2*n
        self.A[4,0] = 3*xy*mu_div_norm_rd5
        self.A[4,1] = 3*x[1,0]*mu_div_norm_rd5 - mu_div_norm_rd3
        self.A[4,2] = 3*z_minus_rc*mu_div_norm_rd5
        self.A[5,0] = 3*x[0,0]*z_minus_rc*mu_div_norm_rd5
        self.A[5,1] = 3*x[1,0]*z_minus_rc*mu_div_norm_rd5
        self.A[5,2] = 3*z_minus_rc**2*mu_div_norm_rd5 - mu_div_norm_rd3 + n**2
        self.A[5,3] = 2*n

    def Predict(self, dt, n, rc, x):
        """
        Predicts the state vector and the covariance matrix.

        Parameters
        ----------
        dt : float
             Time step since the last iteration
        n  : float
             Mean motion of the main satellite
        rc : float
             Current orbit radius of the chief satellite
        x  : matrix
             The state vector, which contains the satellite's
             position and velocity
        """
        #self.x_pred = self.A*x
        iterations = 1#2#100
        h = dt/iterations
        self.x_pred = self.RK4(x, h, iterations, n, rc)
        self.P_pred = self.A*self.P_up*transpose(self.A) + self.Q

    def Update(self, s1, s2, z):
        """
        The update function of the kalman filter. This function updates
        the state vector combining the model with the measurements.

        Parameters
        ----------
        s1 : matrix
             Position vector of the s1 satellite
        s2 : matrix
             Position vector of the s2 satellite
        z  : matrix
             Measurement vector with the measured position of the satellite
        """
        self.updateH(self.x_pred, s1, s2)
        #print("H: {}\nP_pred: {}\nH_tr: {}\nR: {}".format(self.H, self.P_pred, self.H_tr, self.R))
        self.S = self.H*self.P_pred*self.H_tr + self.R
        self.K = self.P_pred*self.H_tr*inv(self.S)
        self.y = z - self.H*self.x_pred
        self.x_up = self.x_pred + self.K*self.y
        self.P_up = (self.I - self.K*self.H)*self.P_pred
        #print("x_up: {}".format(self.x_up))

    def getPosition(self):
        """
        Returns the estimated position of the satellite.
        """
        return self.x_up[0:3]

    def getVelocity(self):
        """
        Returns the estimated velocity of the satellite.
        """
        return self.x_up[3:6]

    def getScalarK(self):
        """Returns a scalar version of the kalman gain."""
        s_K = self.K[0,0] + self.K[1,1] + self.K[2,2]
        return s_K/3

    def getScalarP(self):
        """Returns a scalar version of the covariance matrix."""
        s_P = self.P_up[0,0] + self.P_up[1,1] + self.P_up[2,2]
        s_P += self.P_up[3,3] + self.P_up[4,4] + self.P_up[5,5]
        return s_P/6

    def getPdiag(self):
        """Returns the trace of the covariance matrix."""
        return [self.P_up[0,0], self.P_up[1,1], self.P_up[2,2], self.P_up[3,3], self.P_up[4,4], self.P_up[5,5]]
