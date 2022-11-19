from matplotlib import pyplot as plt
import numpy as np
import math
from scipy import interpolate
from scipy.integrate import odeint
#from mpl_toolkits.mplot3d import art3d
import mpl_toolkits.mplot3d.axes3d as p3
from controller_fdbklin_2D import Controller


class parameters:
    def __init__(self):
        self.m_q = 0.468
        self.m_l = 0.01
        self.Ixx = 4.856*1e-3
        self.Iyy = 4.856*1e-3
        self.Izz = 8.801*1e-3
        self.g = 9.81
        self.l = 0.225
        self.cable_l = 0.5
        self.K = 2.980*1e-6
        self.b = 1.14*1e-7
        self.Ax = 0.25*0
        self.Ay = 0.25*0
        self.Az = 0.25*0
        self.pause = 0.01
        self.fps = 30
        self.K_z = np.array([1, 2])
        self.K_phi = np.array([4, 4])
        self.K_theta = np.array([10, 5])
        self.K_psi = np.array([-8, -8])
        self.Kp = np.array([5, 5, 5])
        self.Kd = np.array([4, 4, 4])
        self.Kdd = np.array([1.00, 1.00, 1.00])
        self.Ki = np.array([1.5, 1.5, 1.5])
        self.Kx = np.array([2, 1])
        self.Kz = np.array([2, 1])

        omega = 1
        speed = omega*np.sqrt(1/self.K)
        dspeed = 0.05*speed
        self.u1 = 0
        self.theta_d = 0
        self.tau_theta = 0

def cos(angle):
    return np.cos(angle)

def sin(angle):
    return np.sin(angle);

def rotation(phi,theta,psi):

    R_x = np.array([
        [1,            0,         0],
        [0,     cos(phi), -sin(phi)],
        [0,     sin(phi),  cos(phi)]

    ])

    R_y = np.array([
        [cos(theta),  0, sin(theta)],
        [0,           1,          0],
        [-sin(theta),  0, cos(theta)]
    ])

    R_z = np.array( [
        [cos(psi), -sin(psi), 0],
        [sin(psi),  cos(psi), 0],
        [0,            0,         1]
    ])

    R_temp = np.matmul(R_y,R_x)
    R = np.matmul(R_z,R_temp)
    return R


def animate(t,Xpos,Xang,Phi_l,parms) -> None:
    """
    This function animates the drone simulation.
    :param t: integration time step
    :param Xpos: Gets the 3-D linear positions of the drone: [x, y, z]
    :param Xang: Gets the 3-D angular positions of the drone: [theta, psi, phi] in radians
    :return: The animation of drone hovering
    """
    t_interp = np.arange(t[0],t[len(t)-1], 1 / parms.fps)
    [m, n] = np.shape(Xpos)
    shape = (len(t_interp), n)
    Xpos_interp = np.zeros(shape)
    Xang_interp = np.zeros(shape)
    phi_l_interp = np.zeros(shape)
    l = parms.l

    for i in range(0, n):
        fpos = interpolate.interp1d(t, Xpos[:,i])
        Xpos_interp[:,i] = fpos(t_interp)
        fang = interpolate.interp1d(t, Xang[:,i])
        Xang_interp[:,i] = fang(t_interp)
    f_phil = interpolate.interp1d(t, Phi_l)
    phi_l_interp = f_phil(t_interp)


    axle_x = np.array([[-l/2, 0, 0],
                       [l/2, 0, 0]])
    axle_y = np.array([[0, -l/2,  0],
                       [0, l/2,   0]])


    [p2,q2] = np.shape(axle_x)
    new_load_pos = np.zeros((2, 1))

    for ii in range(0,len(t_interp)):
        x = Xpos_interp[ii,0]
        y = Xpos_interp[ii,1]
        theta = Xang_interp[ii,0]
        ang = np.array([theta])
        phi_l = phi_l_interp[ii]
        # R = self.Rotation_matrix(ang)
        R = np.array([
            [cos(theta), -sin(theta)],
            [sin(theta),  cos(theta)]])

        R_y_phil = np.array([
            [cos(phi_l), -sin(phi_l)],
            [sin(phi_l),  cos(phi_l)]])

        new_axle_x = np.zeros((p2,q2))
        for i in range(0,p2):
            r_body = axle_x[i,:]
            r_world = R.dot(r_body)
            new_axle_x[i, :] = r_world

        new_axle_x = np.array([x, y, z]) + new_axle_x
        # print(new_axle_x)

        new_axle_y = np.zeros((p2,q2))
        for i in range(0,p2):
            r_body = axle_y[i,:]
            r_world = R.dot(r_body)
            new_axle_y[i, :] = r_world

        new_axle_y = np.array([x, y, z]) + new_axle_y
        # print(new_axle_y)
        # R_load = R_y_phil.dot(R_x_thetal)
        # length = np.array([0, 0, -cable_l])
        # length = length.reshape(len(length),1)
        # pos = np.array([x, y, z]).reshape(3,1)
        # new_load_pos = pos + R_load.dot(length)
        new_load_pos[0] = x - cable_l*sin(phi_l)
        new_load_pos[1] = z - cable_l*cos(phi_l)

        ax = p3.Axes3D(fig, auto_add_to_figure=False)
        fig.add_axes(ax)
        axle1, = ax.plot(new_axle_x[:, 0],new_axle_x[:, 1],new_axle_x[:, 2], 'ro-', linewidth=3)
        axle2, = ax.plot(new_axle_y[:, 0], new_axle_y[:, 1], new_axle_y[:, 2], 'bo-', linewidth=3)
        track, = plt.plot(x, y, z, color='black',marker='o',markersize=2)
        load, = plt.plot([x, new_load_pos[0]], [z, new_load_pos[1]], color='black',marker='o',markersize=5)

        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_zlim(-1, 2.5)
        ax.view_init(azim=-72, elev=20)

        plt.pause(parms.pause)

    plt.close()


def eom(X,t,m_q,m_l,Ixx,Iyy,Izz,g,l,cable_l,K,b,Ax,Ay,Az,u1,theta_d,tau_theta):

    i = 0;
    x = X[i]; i +=1;
    z = X[i]; i +=1;
    phi_l = X[i]; i += 1
    theta = X[i]; i +=1;
    vx = X[i]; i +=1;
    vz = X[i]; i +=1;
    phi_ldot = X[i]; i += 1
    thetadot = X[i]; i +=1;
    ax = cable_l*g*m_l*sin(phi_l)*cos(phi_l)/(cable_l*m_l*sin(phi_l)**2 + cable_l*m_l*cos(phi_l)**2 + cable_l*m_l + cable_l*m_q) + m_l*(-cable_l*m_l*phi_ldot**2*cos(phi_l) - g*(m_l + m_q) + u1)*sin(phi_l)*cos(phi_l)/(m_l**2*sin(phi_l)**2 + m_l**2*cos(phi_l)**2 + m_l**2 + m_l*m_q*sin(phi_l)**2 + m_l*m_q*cos(phi_l)**2 + 2*m_l*m_q + m_q**2) + (-cable_l*m_l*phi_ldot**2*sin(phi_l) + theta_d*u1)*(m_l*sin(phi_l)**2 + m_l + m_q)/(m_l**2*sin(phi_l)**2 + m_l**2*cos(phi_l)**2 + m_l**2 + m_l*m_q*sin(phi_l)**2 + m_l*m_q*cos(phi_l)**2 + 2*m_l*m_q + m_q**2)
    az = -cable_l*g*m_l*sin(phi_l)**2/(cable_l*m_l*sin(phi_l)**2 + cable_l*m_l*cos(phi_l)**2 + cable_l*m_l + cable_l*m_q) + m_l*(-cable_l*m_l*phi_ldot**2*sin(phi_l) + theta_d*u1)*sin(phi_l)*cos(phi_l)/(m_l**2*sin(phi_l)**2 + m_l**2*cos(phi_l)**2 + m_l**2 + m_l*m_q*sin(phi_l)**2 + m_l*m_q*cos(phi_l)**2 + 2*m_l*m_q + m_q**2) + (m_l*cos(phi_l)**2 + m_l + m_q)*(-cable_l*m_l*phi_ldot**2*cos(phi_l) - g*(m_l + m_q) + u1)/(m_l**2*sin(phi_l)**2 + m_l**2*cos(phi_l)**2 + m_l**2 + m_l*m_q*sin(phi_l)**2 + m_l*m_q*cos(phi_l)**2 + 2*m_l*m_q + m_q**2)

    phi_lddot = -u1*sin(phi_l - theta)/(cable_l*(2*m_l + m_q))
    thetaddot = tau_theta/Iyy

    dXdt = np.array([vx, vz, phi_ldot, thetadot, ax, az, phi_lddot, thetaddot])
    return dXdt


def figure8(x0,y0,h,t0,tN):

    N = int((tN-t0)/h) + 1;
    t = np.linspace(t0, tN,N)
    # print(len(t))
    T = t[N-1];
    A = 0.5;
    B = A;
    a = 2;
    b = 1;
    pi = np.pi
    tau = 2*pi*(-15*(t/T)**4+6*(t/T)**5+10*(t/T)**3);
    taudot = 2*pi*(-15*4*(1/T)*(t/T)**3+6*5*(1/T)*(t/T)**4+10*3*(1/T)*(t/T)**2);
    tauddot = 2*pi*(-15*4*3*(1/T)**2*(t/T)**2 + 6*5*4*(1/T)**2*(t/T)**3+10*3*2*(1/T)**2*(t/T));

    x = x0+A*sin(a*tau);
    y = y0+B*cos(b*tau);
    xdot =  A*a*cos(a*tau)*taudot;
    ydot = -B*b*sin(b*tau)*taudot;
    xddot = -A*a*a*sin(a*tau)*taudot+A*a*cos(a*tau)*tauddot;
    yddot = -B*b*b*sin(b*tau)*taudot-B*b*sin(b*tau)*tauddot;

    if (0): #code to check the curve
        plt.figure(1)
        plt.plot(x,y)
        plt.ylabel("y")
        plt.xlabel("x");
        plt.title("Plot of trajectory")
        plt.show(block=False)
        plt.pause(2)
        plt.close()
    return t, x,y,xdot,ydot,xddot,yddot


parms = parameters()
h = 0.005
t0 = 0
tN = 10
N = int((tN-t0)/h) + 1
t = np.linspace(t0, tN, N)
T = t[N-1]

x0_l = 0;
y0_l = 0;

t, x_ref,z_ref,v_x_ref,v_z_ref, \
a_x_ref,a_z_ref  = figure8(x0_l,y0_l,h,t0,tN)

x0 = x_ref[0]; z0 = z_ref[0]
vx0 = 0; vz0 = 0
theta0 = np.deg2rad(0)
thetadot0 = 0
phi_l0 = 0
phidot_l0 = 0

# t = np.linspace(0, 1, 101)
X0 = np.array([x0, z0, phi_l0, theta0, vx0, vz0, phidot_l0, thetadot0], dtype='float64')
# X = odeint(eom, X0, t, args=all_parms)

mm = len(X0)
X_VAL = [X0[0]]; Z = [X0[1]]
THETA = [X0[3]]
VX = [X0[4]]; VZ = [X0[5]]
THETADOT = [X0[7]]
PHI_L = [X0[2]]
PHI_LDOT = [X0[6]]
KE = []; PE = []; TE = []
omega_x=[]; omega_y=[]; omega_z = []
omega_body_x=[]; omega_body_y=[]; omega_body_z = []
m_q = parms.m_q
m_l = parms.m_l
cable_l = parms.cable_l
g = parms.g
Ixx = parms.Ixx
Iyy = parms.Iyy
Izz = parms.Izz
# X_pos = [np.array([X0[0], X0[1], X0[2]])]
# X_ang = [np.array([X0[6], X0[7], X0[8]])]
# x_ref = np.zeros(N)
# v_x_ref = np.zeros(N)
# a_x_ref = np.zeros(N)
# z_ref = np.zeros(N)
# v_z_ref = np.zeros(N)
# a_z_ref = np.zeros(N)
phi_ldes = np.zeros(N)
theta_des = np.zeros(N)
thetadot_des = np.zeros(N)
K_z = parms.K_z
K_phi = parms.K_phi
K_theta = parms.K_theta
K_psi = parms.K_psi
Kp = parms.Kp
Kd = parms.Kd
Kdd = parms.Kdd
Ki = parms.Ki
Kx = parms.Kx
Kz = parms.Kz
A = np.array([parms.Ax, parms.Ay, parms.Az])
k = parms.K
b_drag_const = parms.b
l = parms.l
# omega = np.array([parms.omega1, parms.omega2, parms.omega3, parms.omega4])
OMEGA = np.zeros((len(t), 4))
# OMEGA[0] = omega
X_POS = np.zeros((len(t), 2))
X_ANG = np.zeros((len(t), 2))
X_POS[0, 0] = X0[0]
X_POS[0, 1] = X0[1]

X_ANG[0, 0] = X0[2]
X_ANG[0, 1] = X0[3]


for i in range(0, N-1):
    j = 0
    control = Controller(K_z, K_psi, K_theta, K_phi, Kp, Kd, Kdd, Ki, Kx, Kz, A, k, l, b_drag_const)
    t_temp = np.array([t[i], t[i+1]], dtype='float64')
    # desired_state = control.get_desired_positions(t_temp, desired_traj_values)
    desired_state = np.array([x_ref[i], z_ref[i], v_x_ref[i], v_z_ref[i], a_x_ref[i], a_z_ref[i]])
    all_parms = (parms.m_q,parms.m_l,parms.Ixx,parms.Iyy,parms.Izz,parms.g,parms.l,parms.cable_l,
                 parms.K,parms.b,parms.Ax,parms.Ay,parms.Az,
                 parms.u1,parms.theta_d,parms.tau_theta)
    X = odeint(eom, X0, t_temp, args=all_parms)
    x_des, z_des, xdot_des, zdot_des, xddot_des, zddot_des = desired_state
    x = X[1][0]
    z = X[1][1]
    phil = X[1][2]
    vx = X[1][4]
    vz = X[1][5]
    phildot = X[1][6]
    theta = X[1][3]
    thetadot = X[1][7]
    # x, vx, z, vz, x_des, z_des, xdot_des, zdot_des, xddot_des, zddot_des, phi_l, phi_ldot

    # ang = np.array([[X[1][3], X[1][9]], [X[1][4], X[1][10]], [X[1][5], X[1][11]]], dtype='float64')
    # translation = np.array([[X[1][0], X[1][6]], [X[1][1], X[1][7]], [X[1][2], X[1][8]]], dtype='float64')
    # vertical = translation[2]
    # torques, theta_temp, psi_temp = control.get_forces(vertical, ang, desired_state)
    # omega = control.get_action(desired_state, ang, translation)
    # omega_temp = omega.reshape(4,)
    # parms.omega1, parms.omega2, parms.omega3, parms.omega4 = omega
    theta_d, u1 = control.get_desired_positions(x, vx, z, vz, x_des, z_des, xdot_des, zdot_des, xddot_des, zddot_des, phil, phildot)
    parms.theta_d = theta_d
    parms.u1 = u1
    ang = np.array([theta, thetadot], dtype='float64')
    THETA_DES_STATE = np.array([theta_d], dtype='float64')
    T_theta = control.choose_action(ang, THETA_DES_STATE)
    parms.tau_theta = T_theta
    X0 = X[1]
    X_VAL.append(X[1][j]); j+=1;
    Z.append(X[1][j]); j+=1;
    PHI_L.append(X[1][j]); j+=1;
    THETA.append(X[1][j]); j+=1;
    VX.append(X[1][j]); j+=1;
    VZ.append(X[1][j]); j+=1;
    PHI_LDOT.append(X[1][j]); j+=1;
    THETADOT.append(X[1][j]); j+=1;
    X_POS[i+1, 0] = X[1][0]
    X_POS[i+1, 1] = X[1][1]

    X_ANG[i+1, 0] = X[1][2]
    X_ANG[i+1, 1] = X[1][3]
    # X_pos.append(np.array([X_VAL[i], Y[i], Z[i]]))
    # X_ang.append(np.array([PHI[i], THETA[i], PSI[i]]))
    # OMEGA[i+1] = omega_temp

plt.figure(1)
plt.subplot(3,1,1)
plt.plot(t,X_VAL);
plt.plot(t,Z)
plt.ylabel('linear position');
plt.legend(['x_pos', 'y_pos', 'z_pos'])
plt.subplot(3,1,2)
plt.plot(t,VX);
plt.plot(t,VZ)
plt.xlabel('time')
plt.ylabel('linear velocity');
plt.legend(['x_vel', 'y_vel', 'z_vel'])
plt.subplot(3,1,3)
plt.plot(t,THETA)
plt.xlabel('time')
plt.ylabel('angular position');
plt.legend(['Angle phi', 'Angle theta', 'Angle psi'])

# plt.figure(2)
# plt.subplot(2,1,1)
# plt.plot(t,phi);
# plt.plot(t,theta)
# plt.plot(t,psi)
# plt.ylabel('angular position');
# plt.subplot(2,1,2)
# plt.plot(t,phidot);
# plt.plot(t,thetadot)
# plt.plot(t,psidot)
# plt.xlabel('time')
# plt.ylabel('angular velocity');
# #

#
# ax=plt.figure(4)
# plt.subplot(1,1,1)
# plt.plot(t,OMEGA);
# plt.ylabel('omega');
# # plt.subplot(2,1,2)
# # plt.plot(t,omega_body_x);
# # plt.plot(t,omega_body_y);
# # plt.plot(t,omega_body_z);
# ax.legend(['rotor1', 'rotor2','rotor3', 'rotor4'])
# # plt.ylabel('omega body');
# plt.xlabel('time')
#
#
fig = plt.figure(5)
plt.subplot(1, 1, 1)
plt.plot(t, PHI_L)
plt.legend(['Theta_l', 'Phi_l'])
plt.ylabel('Load swing angles')
plt.xlabel('time'
           )
#plt.show()
plt.show(block=False)
plt.pause(60)
plt.close()
#
fig = plt.figure(6)
animate(t,X_POS,X_ANG,PHI_L,parms)
