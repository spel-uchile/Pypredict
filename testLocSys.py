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
from pyorbital import tlefile
from sat import Sat
from localizationSystem import Loc
from dpl import Dpl
from numpy import arctan, cos, transpose, matrix, sin, sqrt, tan
from matplotlib import rcParams
from matplotlib.pyplot import plot, show, subplots, subplots_adjust
from datetime import datetime, timedelta

SUCHAI = Sat(name="SUCHAI", tle=tlefile.read("SUCHAI", "TLE/cubesat.txt"), cat="CubeSat")
HODOYOSHI3 = Sat(name="HODOYOSHI-3", tle=tlefile.read("HODOYOSHI-3", "TLE/resource.txt"), cat="Earth Resources")
HODOYOSHI4 = Sat(name="HODOYOSHI-4", tle=tlefile.read("HODOYOSHI-4", "TLE/resource.txt"), cat="Earth Resources")
CUBESATXI_IV = Sat(name="CUBESAT XI-IV (CO-57)", tle=tlefile.read("CUBESAT XI-IV (CO-57)", "TLE/cubesat.txt"), cat="CubeSat")

date = datetime.utcnow()
SUCHAI.updateOrbitalParameters3(date)
HODOYOSHI3.updateOrbitalParameters3(date)
HODOYOSHI4.updateOrbitalParameters3(date)
CUBESATXI_IV.updateOrbitalParameters3(date)

dpl = Dpl()
SUCHAI_mass = 3
FE1_mass = 0.08
v = [2, 0, 0]
FE1 = dpl.deploy("Femto", SUCHAI, SUCHAI_mass, FE1_mass, "FE1", v, date)
FE1.updateOrbitalParameters3(date)
up_dt = 2
pred_dt = 0.005
r0_i = FE1.getXYZ()
r1_i = SUCHAI.getXYZ()
r2_i = HODOYOSHI3.getXYZ()
r3_i = HODOYOSHI4.getXYZ()
r4_i = CUBESATXI_IV.getXYZ()
v = FE1.getInertialVel()
x = matrix([[0], [0], [0], [0], [0], [0]])
x[0:3] = r0_i[0:3]
x[3:6] = v[0:3]
loc = Loc()
test = matrix([[0], [0], [0]])
real_x = []
real_y = []
real_z = []
real_vx = []
real_vy = []
real_vz = []
est_x = []
est_y = []
est_z = []
est_vx = []
est_vy = []
est_vz = []
me_x = []
me_y = []
me_z = []
x_diff = []
y_diff = []
z_diff = []
vx_diff = []
vy_diff = []
vz_diff = []
P_x = []
P_y = []
P_z = []
P_vx = []
P_vy = []
P_vz = []
K = []
P = []
t = []
t_m = []

def predict(incl, RAAN, theta, w, r, v):
    t.append(i*up_dt + j*pred_dt)
    loc.initTransformationData(incl, RAAN, theta, w, r)
    loc.mainSat_v = loc.Q_po*loc.Q_ip*v
    loc.ekf.Predict(loc.ekf.x_up)
    loc.ekf.x_up = loc.ekf.x_pred
    loc.ekf.P_up = loc.ekf.P_pred
    FE1.updateOrbitalParameters3(date)
    real_pos = FE1.getXYZ()
    real_x.append(real_pos[0,0])
    real_y.append(real_pos[1,0])
    real_z.append(real_pos[2,0])
    real_vel = FE1.getInertialVel()
    real_vx.append(real_vel[0,0])
    real_vy.append(real_vel[1,0])
    real_vz.append(real_vel[2,0])
    est_pos = loc.getEstimatedPos()
    est_x.append(est_pos[0,0])
    est_y.append(est_pos[1,0])
    est_z.append(est_pos[2,0])
    est_vel = loc.getEstimatedVel()
    est_vx.append(est_vel[0,0])
    est_vy.append(est_vel[1,0])
    est_vz.append(est_vel[2,0])
    x_diff.append(real_pos[0,0] - est_pos[0,0])
    y_diff.append(real_pos[1,0] - est_pos[1,0])
    z_diff.append(real_pos[2,0] - est_pos[2,0])
    vx_diff.append(real_vel[0,0] - est_vel[0,0])
    vy_diff.append(real_vel[1,0] - est_vel[1,0])
    vz_diff.append(real_vel[2,0] - est_vel[2,0])
    K.append(loc.ekf.getScalarK())
    P.append(loc.ekf.getScalarP())
    Pdiag = loc.ekf.getPdiag()
    P_x.append(sqrt(Pdiag[0]))
    P_y.append(sqrt(Pdiag[1]))
    P_z.append(sqrt(Pdiag[2]))
    P_vx.append(sqrt(Pdiag[3]))
    P_vy.append(sqrt(Pdiag[4]))
    P_vz.append(sqrt(Pdiag[5]))

t0 = SUCHAI.getTnow()
#loc.initTransformationData(SUCHAI.incl, SUCHAI.RAAN, SUCHAI.theta,
#        SUCHAI.w, r1_i)
#loc.mainSat_v = loc.Q_po*loc.Q_ip*SUCHAI.getInertialVel()
#loc.ekf.initPropagationMatrix(SUCHAI.n, pred_dt)
#loc.ekf.x_up = matrix([[0], [0], [0], [0], [0], [0]])
#loc.ekf.K = matrix([[1, 1, 1],
#                    [1, 1, 1],
#                    [1, 1, 1],
#                    [1, 1, 1],
#                    [1, 1, 1],
#                    [1, 1, 1]])
#loc.ekf.x_up[0:3] = loc.transformPosition(x[0:3])
#loc.ekf.x_up[3:6] = loc.transformVelocity(x[3:6])
for i in range(0, 600):
    #for j in range(1, int(up_dt/pred_dt)):
    #    SUCHAI.updateWithDragEffect(t0 + (i-1)*up_dt + j*pred_dt)
    #    r1_i = SUCHAI.getXYZ()
    #    SUCHAI_v = SUCHAI.getInertialVel()
    #    predict(SUCHAI.incl, SUCHAI.RAAN, SUCHAI.theta,
    #            SUCHAI.w, r1_i, SUCHAI_v)
        #print((i-1)*up_dt + j*pred_dt)
    #print(i*up_dt)
    if (i > 0):
        x[0:3] = loc.getEstimatedPos()
        x[3:6] = loc.getEstimatedVel()
    date = date + timedelta(seconds=pred_dt)
    t.append(i*up_dt)
    t_m.append(i*up_dt)
    FE1.updateOrbitalParameters3(date)
    SUCHAI.updateOrbitalParameters3(date)
    HODOYOSHI3.updateOrbitalParameters3(date)
    HODOYOSHI4.updateOrbitalParameters3(date)
    CUBESATXI_IV.updateOrbitalParameters3(date)
    r1_i = SUCHAI.getXYZ()
    r2_i = HODOYOSHI3.getXYZ()
    r3_i = HODOYOSHI4.getXYZ()
    r4_i = CUBESATXI_IV.getXYZ()
    real_pos = FE1.getXYZ()
    r0_i = real_pos
    SUCHAI_v = SUCHAI.getInertialVel()
    test = None#real_pos
    loc.estimateLocation(SUCHAI.incl, SUCHAI.RAAN, SUCHAI.theta, SUCHAI.w,
                         SUCHAI_v, pred_dt, SUCHAI.n, r0_i, r1_i, r2_i,
                         r3_i, r4_i, x, test)
    real_x.append(real_pos[0,0])
    real_y.append(real_pos[1,0])
    real_z.append(real_pos[2,0])
    real_vel = FE1.getInertialVel()
    real_vx.append(real_vel[0,0])
    real_vy.append(real_vel[1,0])
    real_vz.append(real_vel[2,0])
    est_pos = loc.getEstimatedPos()
    est_x.append(est_pos[0,0])
    est_y.append(est_pos[1,0])
    est_z.append(est_pos[2,0])
    est_vel = loc.getEstimatedVel()
    est_vx.append(est_vel[0,0])
    est_vy.append(est_vel[1,0])
    est_vz.append(est_vel[2,0])
    x[0:3] = est_pos[0:3]
    x[3:6] = est_vel[0:3]
    x_diff.append(real_pos[0,0] - est_pos[0,0])
    y_diff.append(real_pos[1,0] - est_pos[1,0])
    z_diff.append(real_pos[2,0] - est_pos[2,0])
    vx_diff.append(real_vel[0,0] - est_vel[0,0])
    vy_diff.append(real_vel[1,0] - est_vel[1,0])
    vz_diff.append(real_vel[2,0] - est_vel[2,0])
    z_iner = loc.orbital2inertial(loc.z[0:3])
    me_x.append(z_iner[0,0])
    me_y.append(z_iner[1,0])
    me_z.append(z_iner[2,0])
    K.append(loc.ekf.getScalarK())
    P.append(loc.ekf.getScalarP())
    Pdiag = loc.ekf.getPdiag()
    P_x.append(sqrt(Pdiag[0]))
    P_y.append(sqrt(Pdiag[1]))
    P_z.append(sqrt(Pdiag[2]))
    P_vx.append(sqrt(Pdiag[3]))
    P_vy.append(sqrt(Pdiag[4]))
    P_vz.append(sqrt(Pdiag[5]))
    print("Distance real-estimated: {} m".format(sqrt((real_pos[0,0]-est_pos[0,0])**2 + (real_pos[1,0] - est_pos[1,0])**2 + (real_pos[2,0] - est_pos[2,0])**2)))
    for j in range(1, int(up_dt/pred_dt)):
        date = date + timedelta(seconds=pred_dt)
        SUCHAI.updateOrbitalParameters3(date)
        r1_i = SUCHAI.getXYZ()
        SUCHAI_v = SUCHAI.getInertialVel()
        predict(SUCHAI.incl, SUCHAI.RAAN, SUCHAI.theta,
                SUCHAI.w, r1_i, SUCHAI_v)

print("{}{}{}{:0.2f}".format("Qw:\n", loc.ekf.Qw, "   P min: ", min(P)))
rcParams["axes.formatter.useoffset"] = False
fig_pos, (ax_pv_x, ax_pv_y, ax_pv_z) = subplots(3, 2, sharey=False)
ax_pv_x[0].plot(t, est_x, label="Est.")
ax_pv_x[0].plot(t, real_x, label="Real", color="tab:red")
ax_pv_x[0].plot(t_m, me_x, label="Me.", marker='.',
        linestyle = 'None', color="tab:green")
ax_pv_x[0].legend()
ax_pv_x[0].set(ylabel="Position in x axis [m]",
        title="Estimated and real positions over time")
ax_pv_x[0].grid()
ax_pv_x[1].plot(t, est_vx, label="Est.")
ax_pv_x[1].plot(t, real_vx, label="Real", color="tab:red")
ax_pv_x[1].legend()
ax_pv_x[1].set(ylabel="Velocity in x [m/s]",
        title="Estimated and real velocities over time")
ax_pv_x[1].grid()

ax_pv_y[0].plot(t, est_y, label="Est.")
ax_pv_y[0].plot(t, real_y, label="Real", color="tab:red")
ax_pv_y[0].plot(t_m, me_y, label="Me.", marker='.',
        linestyle = 'None', color="tab:green")
ax_pv_y[0].legend()
ax_pv_y[0].set(ylabel="Position in y axis [m]")
ax_pv_y[0].grid()
ax_pv_y[1].plot(t, est_vy, label="Est.")
ax_pv_y[1].plot(t, real_vy, label="Real", color="tab:red")
ax_pv_y[1].legend()
ax_pv_y[1].set(ylabel="Velocity in y [m/s]")
ax_pv_y[1].grid()

ax_pv_z[0].plot(t, est_z, label="Est.")
ax_pv_z[0].plot(t, real_z, label="Real", color="tab:red")
ax_pv_z[0].plot(t_m, me_z, label="Me.", marker='.',
        linestyle = 'None', color="tab:green")
ax_pv_z[0].legend()
ax_pv_z[0].set(ylabel="Position in z axis [m]",
        xlabel="Time [s]")
ax_pv_z[0].grid()
ax_pv_z[1].plot(t, est_vz, label="Est.")
ax_pv_z[1].plot(t, real_vz, label="Real", color="tab:red")
ax_pv_z[1].legend()
ax_pv_z[1].set(xlabel="Time [s]", ylabel="Velocity in z [m/s]")
ax_pv_z[1].grid()
subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.96, wspace=0.2, hspace=0.2)

fig_diff, (ax_diff_x, ax_diff_y, ax_diff_z) = subplots(3, 2, sharey=False)
ax_diff_x[0].plot(t, x_diff, label="Real-Est.")
ax_diff_x[0].plot(t, P_x, label="±Sqrt(Px)")
ax_diff_x[0].plot(t, [ -x for x in P_x], color="tab:orange")
ax_diff_x[0].legend()
ax_diff_x[0].set(ylabel="Difference x [m]",
        title="Difference between real and estimated position")
ax_diff_x[0].grid()
ax_diff_x[1].plot(t, vx_diff, label="Real-Est.")
ax_diff_x[1].plot(t, P_vx, label="±Sqrt(Pvx)")
ax_diff_x[1].plot(t, [ -x for x in P_vx], color="tab:orange")
ax_diff_x[1].legend()
ax_diff_x[1].set(ylabel="Difference vx [m/s]",
        title="Difference between real and estimated velocity")
ax_diff_x[1].grid()

ax_diff_y[0].plot(t, y_diff, label="Real-Est.")
ax_diff_y[0].plot(t, P_y, label="±Sqrt(Py)")
ax_diff_y[0].plot(t, [ -x for x in P_y], color="tab:orange")
ax_diff_y[0].legend()
ax_diff_y[0].set(ylabel="Difference y [m]")
ax_diff_y[0].grid()
ax_diff_y[1].plot(t, vy_diff, label="Real-Est.")
ax_diff_y[1].plot(t, P_vy, label="±Sqrt(Pvy)")
ax_diff_y[1].plot(t, [ -x for x in P_vy], color="tab:orange")
ax_diff_y[1].legend()
ax_diff_y[1].set(ylabel="Difference vy [m/s]")
ax_diff_y[1].grid()

ax_diff_z[0].plot(t, z_diff, label="Real-Est.")
ax_diff_z[0].plot(t, P_z, label="±Sqrt(Pz)")
ax_diff_z[0].plot(t, [ -x for x in P_z], color="tab:orange")
ax_diff_z[0].legend()
ax_diff_z[0].set(xlabel="Time [s]", ylabel="Difference z [m]")
ax_diff_z[0].grid()
ax_diff_z[1].plot(t, vz_diff, label="Real-Est.")
ax_diff_z[1].plot(t, P_vz, label="±Sqrt(Pvz)")
ax_diff_z[1].plot(t, [ -x for x in P_vz], color="tab:orange")
ax_diff_z[1].legend()
ax_diff_z[1].set(xlabel="Time [s]", ylabel="Difference vz [m/s]")
ax_diff_z[1].grid()
subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.96, wspace=0.2, hspace=0.2)

fig_KP, ax_KP = subplots(2, 1, sharey=False)
ax_KP[0].plot(t, K)
ax_KP[0].set(ylabel="Kalman gain",
        title="Scalar Kalman gain over time")
ax_KP[0].grid()
ax_KP[1].plot(t, P)
ax_KP[1].set(xlabel="Time [s]", ylabel="Covariance",
        title="Scalar covariance matrix over time")
ax_KP[1].grid()
subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.96, wspace=0.2, hspace=0.2)
show()
