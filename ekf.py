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
from numpy import matrix, sqrt, transpose
from numpy.linalg import inv

class EKF(object):
    __slots__ = ["A", "H", "H_tr", "I", "K",
                 "P_pred", "P_up", "Qv", "Qw",
                 "S", "x_pred", "x_up", "y"]
    def __init__(self):
        self.initMatrices()

    def __call__(self):
        return self

    def __str__(self):
        return "EKF applied to estimate a satellite's position and velocity"

    def newIteration(self, dt, n, x, z):
        """
        Runs one iteration of the Kalman filter's algorithm.

        Parameters
        ----------
        dt : float
             Time step since the last iteration
        n : float
            Mean motion of the main satellite
        x : matrix
            State vector with the last estimated position and velocity
        z : matrix
            Measurement of the satellite's new position
        """
        self.initPropagationMatrix(dt, n)
        self.Predict(x)
        self.Update(z)

    def initMatrices(self):
        """Initialize all the matrices needed to run the kalman filter."""
        self.initCovarianceMatrix()
        self.initH()
        self.initIdentityMatrix()
        self.initQv()
        self.initQw()

    def initPropagationMatrix(self, n, dt):
        """
        Initialize the propagation matrix A used to propagate the
        state vector.

        Parameters
        ----------
        n : float
            Mean motion of the main satellite
        dt : float
             Time step in seconds since the last iteration
        """
        self.A = matrix([[1, 0, 0, dt, 0, n*dt**2],
                         [0, 1-n**2*dt**2/2, 0, 0, dt, 0],
                         [0, 0, 1 + 3*n**2*dt**2/2, -n*dt**2, 0, dt],
                         [0, 0, 0, 1, 0, 2*n*dt],
                         [0, -n**2*dt, 0, 0, 1, 0],
                         [0, 0, 3*n**2*dt, -2*n*dt, 0, 1]])

    def initH(self):
        """
        Initialize the matrix that describes the observation model.
        """
        self.H = matrix([[1, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0]])
        self.H_tr = transpose(self.H)

    def initCovarianceMatrix(self):
        """Initialize the kalman filter's covariance matrix."""
        self.P_up = matrix([[200**2, 0, 0, 0, 0, 0],
                            [0, 200**2, 0, 0, 0, 0],
                            [0, 0, 200**2, 0, 0, 0],
                            [0, 0, 0, 200**2, 0, 0],
                            [0, 0, 0, 0, 200**2, 0],
                            [0, 0, 0, 0, 0, 200**2]])

    def initIdentityMatrix(self):
        """Initialize the 6x6 identity matrix."""
        self.I = matrix([[1, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 0, 1]])

    def initQv(self):
        """Initialize the covariance matrix of the observation noise."""
        self.Qv = matrix([[3.2**2, 0, 0],
                          [0, 3.2**2, 0],
                          [0, 0, 3.2**2]])

    def initQw(self):
        """Initialize the covariance matrix of the process noise."""
        self.Qw = matrix([[0**2, 0, 0, 0, 0, 0],
                          [0, 0**2, 0, 0, 0, 0],
                          [0, 0, 0**2, 0, 0, 0],
                          [0, 0, 0, 2**2, 0, 0],
                          [0, 0, 0, 0, 2**2, 0],
                          [0, 0, 0, 0, 0, 5**2]])

    def Predict(self, x):
        """
        Predicts the state vector and the covariance matrix.

        Parameters
        ----------
        x : matrix
            The state vector, which contains the satellite's
            position and velocity
        """
        self.x_pred = self.A*x
        self.P_pred = self.A*self.P_up*transpose(self.A) + self.Qw

    def Update(self, z):
        """
        The update function of the kalman filter. This function updates
        the state vector combining the model with the measurements.

        Parameters
        ----------
        z : matrix
            Measurement vector with the measured position of the satellite
        """
        self.S = self.H*self.P_pred*self.H_tr + self.Qv
        self.K = self.P_pred*self.H_tr*inv(self.S)
        self.y = z - self.H*self.x_pred
        self.x_up = self.x_pred + self.K*self.y
        self.P_up = (self.I - self.K*self.H)*self.P_pred

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
