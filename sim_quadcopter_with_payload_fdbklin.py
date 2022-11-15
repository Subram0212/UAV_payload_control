from matplotlib import pyplot as plt
import numpy as np
import math
from scipy import interpolate
from scipy.integrate import odeint
#from mpl_toolkits.mplot3d import art3d
import mpl_toolkits.mplot3d.axes3d as p3
from controller import Controller


class parameters:
    def __init__(self):
        self.m_q = 0.468
        self.m_l = 0.01
        self.Ixx = 4.856*1e-3
        self.Iyy = 4.856*1e-3
        self.Izz = 8.801*1e-3
        self.g = 9.81
        self.l = 0.225
        self.cable_l = 0.2
        self.K = 2.980*1e-6
        self.b = 1.14*1e-7
        self.Ax = 0.25
        self.Ay = 0.25
        self.Az = 0.25
        self.pause = 0.01
        self.fps = 30
        self.K_z = np.array([1, 2])
        self.K_phi = np.array([4, 4])
        self.K_theta = np.array([4, 4])
        self.K_psi = np.array([-8, -8])
        self.Kp = np.array([5, 5, 5])
        self.Kd = np.array([4, 4, 4])
        self.Kdd = np.array([1.00, 1.00, 1.00])
        self.Ki = np.array([1.5, 1.5, 1.5])

        omega = 1
        speed = omega*np.sqrt(1/self.K)
        dspeed = 0.05*speed
        self.omega1 = 600
        self.omega2 = 600
        self.omega3 = 600
        self.omega4 = 600

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


def animate(t,Xpos,Xang,Theta_l, Phi_l, parms) -> None:
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
    theta_l_interp = np.zeros(shape)
    phi_l_interp = np.zeros(shape)
    l = parms.l

    for i in range(0, n):
        fpos = interpolate.interp1d(t, Xpos[:,i])
        Xpos_interp[:,i] = fpos(t_interp)
        fang = interpolate.interp1d(t, Xang[:,i])
        Xang_interp[:,i] = fang(t_interp)
    f_thetal = interpolate.interp1d(t, Theta_l)
    theta_l_interp = f_thetal(t_interp)
    f_phil = interpolate.interp1d(t, Phi_l)
    phi_l_interp = f_phil(t_interp)


    axle_x = np.array([[-l/2, 0, 0],
                       [l/2, 0, 0]])
    axle_y = np.array([[0, -l/2,  0],
                       [0, l/2,   0]])


    [p2,q2] = np.shape(axle_x)
    new_load_pos = np.zeros((3, 1))

    for ii in range(0,len(t_interp)):
        x = Xpos_interp[ii,0]
        y = Xpos_interp[ii,1]
        z = Xpos_interp[ii,2]
        theta = Xang_interp[ii,0]
        psi = Xang_interp[ii,1]
        phi = Xang_interp[ii,2]
        ang = np.array([theta, psi, phi])
        theta_l = theta_l_interp[ii]
        phi_l = phi_l_interp[ii]
        # R = self.Rotation_matrix(ang)
        R = rotation(phi,theta,psi)
        R_y_phil = np.array([
            [cos(phi_l),  0, sin(phi_l)],
            [0,           1,          0],
            [-sin(phi_l),  0, cos(phi_l)]
        ])

        R_x_thetal = np.array([
            [1,            0,         0],
            [0,     cos(theta_l), -sin(theta_l)],
            [0,     sin(theta_l),  cos(theta_l)]

        ])

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
        new_load_pos[0] = x - cable_l*cos(theta_l)*sin(phi_l)
        new_load_pos[1] = y + cable_l*sin(theta_l)
        new_load_pos[2] = z - cable_l*cos(theta_l)*cos(phi_l)

        ax = p3.Axes3D(fig, auto_add_to_figure=False)
        fig.add_axes(ax)
        axle1, = ax.plot(new_axle_x[:, 0],new_axle_x[:, 1],new_axle_x[:, 2], 'ro-', linewidth=3)
        axle2, = ax.plot(new_axle_y[:, 0], new_axle_y[:, 1], new_axle_y[:, 2], 'bo-', linewidth=3)
        track, = plt.plot(x, y, z, color='black',marker='o',markersize=2)
        load, = plt.plot([x, new_load_pos[0]], [y, new_load_pos[1]], [z, new_load_pos[2]], color='black',marker='o',markersize=5)

        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_zlim(-1, 2.5)
        ax.view_init(azim=-72, elev=20)

        plt.pause(parms.pause)

    plt.close()


def eom(X,t,m_q,m_l,Ixx,Iyy,Izz,g,l,cable_l,K,b,Ax,Ay,Az,omega1,omega2,omega3,omega4):

    i = 0;
    x = X[i]; i +=1;
    y = X[i]; i +=1;
    z = X[i]; i +=1;
    theta_l = X[i]; i += 1
    phi_l = X[i]; i += 1
    phi = X[i]; i +=1;
    theta = X[i]; i +=1;
    psi = X[i]; i+=1;
    vx = X[i]; i +=1;
    vy = X[i]; i +=1;
    vz = X[i]; i +=1;
    theta_ldot = X[i]; i+=1
    phi_ldot = X[i]; i += 1
    phidot = X[i]; i +=1;
    thetadot = X[i]; i +=1;
    psidot = X[i]; i+=1;

    A = np.zeros((8,8))
    B = np.zeros((8,1))
    B_ctrl = np.zeros((8,4))

    A[ 0 , 0 ]= 1.0*m_l + 1.0*m_q
    A[ 0 , 1 ]= 0
    A[ 0 , 2 ]= 0
    A[ 0 , 3 ]= 1.0*cable_l*m_l*sin(phi_l)*sin(theta_l)
    A[ 0 , 4 ]= -1.0*cable_l*m_l*cos(phi_l)*cos(theta_l)
    A[ 0 , 5 ]= 0
    A[ 0 , 6 ]= 0
    A[ 0 , 7 ]= 0
    A[ 1 , 0 ]= 0
    A[ 1 , 1 ]= 1.0*m_l + 1.0*m_q
    A[ 1 , 2 ]= 0
    A[ 1 , 3 ]= 1.0*cable_l*m_l*cos(theta_l)
    A[ 1 , 4 ]= 0
    A[ 1 , 5 ]= 0
    A[ 1 , 6 ]= 0
    A[ 1 , 7 ]= 0
    A[ 2 , 0 ]= 0
    A[ 2 , 1 ]= 0
    A[ 2 , 2 ]= 1.0*m_l + 1.0*m_q
    A[ 2 , 3 ]= 1.0*cable_l*m_l*sin(theta_l)*cos(phi_l)
    A[ 2 , 4 ]= 1.0*cable_l*m_l*sin(phi_l)*cos(theta_l)
    A[ 2 , 5 ]= 0
    A[ 2 , 6 ]= 0
    A[ 2 , 7 ]= 0
    A[ 3 , 0 ]= 1.0*cable_l*m_l*sin(phi_l)*sin(theta_l)
    A[ 3 , 1 ]= 1.0*cable_l*m_l*cos(theta_l)
    A[ 3 , 2 ]= 1.0*cable_l*m_l*sin(theta_l)*cos(phi_l)
    A[ 3 , 3 ]= 1.0*cable_l**2*m_l
    A[ 3 , 4 ]= 0
    A[ 3 , 5 ]= 0
    A[ 3 , 6 ]= 0
    A[ 3 , 7 ]= 0
    A[ 4 , 0 ]= -1.0*cable_l*m_l*cos(phi_l)*cos(theta_l)
    A[ 4 , 1 ]= 0
    A[ 4 , 2 ]= 1.0*cable_l*m_l*sin(phi_l)*cos(theta_l)
    A[ 4 , 3 ]= 0
    A[ 4 , 4 ]= 1.0*cable_l**2*m_l*cos(theta_l)**2
    A[ 4 , 5 ]= 0
    A[ 4 , 6 ]= 0
    A[ 4 , 7 ]= 0
    A[ 5 , 0 ]= 0
    A[ 5 , 1 ]= 0
    A[ 5 , 2 ]= 0
    A[ 5 , 3 ]= 0
    A[ 5 , 4 ]= 0
    A[ 5 , 5 ]= 1.0*Ixx
    A[ 5 , 6 ]= 0
    A[ 5 , 7 ]= -1.0*Ixx*sin(theta)
    A[ 6 , 0 ]= 0
    A[ 6 , 1 ]= 0
    A[ 6 , 2 ]= 0
    A[ 6 , 3 ]= 0
    A[ 6 , 4 ]= 0
    A[ 6 , 5 ]= 0
    A[ 6 , 6 ]= 1.0*Iyy*cos(phi)**2 + 1.0*Izz*sin(phi)**2
    A[ 6 , 7 ]= 0.25*(Iyy - Izz)*(sin(2*phi - theta) + sin(2*phi + theta))
    A[ 7 , 0 ]= 0
    A[ 7 , 1 ]= 0
    A[ 7 , 2 ]= 0
    A[ 7 , 3 ]= 0
    A[ 7 , 4 ]= 0
    A[ 7 , 5 ]= -1.0*Ixx*sin(theta)
    A[ 7 , 6 ]= 0.25*(Iyy - Izz)*(sin(2*phi - theta) + sin(2*phi + theta))
    A[ 7 , 7 ]= 1.0*Ixx*sin(theta)**2 + 1.0*Iyy*sin(phi)**2*cos(theta)**2 + 1.0*Izz*cos(phi)**2*cos(theta)**2
    B[ 0 ]= -Ax*vx + K*(sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi))*(omega1**2 + omega2**2 + omega3**2 + omega4**2) - 1.0*cable_l*m_l*phi_ldot*(phi_ldot*sin(phi_l)*cos(theta_l) + theta_ldot*sin(theta_l)*cos(phi_l)) - 1.0*cable_l*m_l*theta_ldot*(phi_ldot*sin(theta_l)*cos(phi_l) + theta_ldot*sin(phi_l)*cos(theta_l))
    B[ 1 ]= -Ay*vy - K*(sin(phi)*cos(psi) - sin(psi)*sin(theta)*cos(phi))*(omega1**2 + omega2**2 + omega3**2 + omega4**2) + 1.0*cable_l*m_l*theta_ldot**2*sin(theta_l)
    B[ 2 ]= -Az*vz + K*(omega1**2 + omega2**2 + omega3**2 + omega4**2)*cos(phi)*cos(theta) - 1.0*cable_l*m_l*phi_ldot*(phi_ldot*cos(phi_l)*cos(theta_l) - theta_ldot*sin(phi_l)*sin(theta_l)) + 1.0*cable_l*m_l*theta_ldot*(phi_ldot*sin(phi_l)*sin(theta_l) - theta_ldot*cos(phi_l)*cos(theta_l)) - g*m_l - g*m_q
    B[ 3 ]= -1.0*m_l*(cable_l**2*phi_ldot**2*cos(theta_l) + g*l*cos(phi_l))*sin(theta_l)
    B[ 4 ]= m_l*(2.0*cable_l**2*phi_ldot*theta_ldot*sin(theta_l) - 1.0*g*l*sin(phi_l))*cos(theta_l)
    B[ 5 ]= Ixx*psidot*thetadot*cos(theta) - K*l*(omega2**2 - omega4**2) + psidot*(psidot*(Iyy - Izz)*(sin(2*phi - theta) + sin(2*phi + theta)) - 2.0*thetadot*(-Iyy + Izz)*cos(2*phi))*cos(theta)/4 - thetadot*(psidot*(-Iyy + Izz)*cos(2*phi)*cos(theta) + thetadot*(Iyy - Izz)*sin(2*phi))/2
    B[ 6 ]= -0.5*Ixx*phidot*psidot*cos(theta) - K*l*(omega1**2 - omega3**2) - 1.0*phidot*(psidot*(Iyy - Izz)*cos(2*phi)*cos(theta) - thetadot*(Iyy - Izz)*sin(2*phi)) + 0.125*psidot*thetadot*(Iyy - Izz)*(cos(2*phi - theta) - cos(2*phi + theta)) - psidot*(0.5*Ixx*phidot*cos(theta) + psidot*(-Ixx + Iyy*sin(phi)**2 + Izz*cos(phi)**2)*sin(theta)*cos(theta) + 0.125*thetadot*(Iyy - Izz)*(cos(2*phi - theta) - cos(2*phi + theta)))
    B[ 7 ]= b*(omega1**2 - omega2**2 + omega3**2 - omega4**2) - phidot*(0.5*psidot*(Iyy - Izz)*(sin(2*phi - theta) + sin(2*phi + theta)) - 1.0*thetadot*(-Iyy + Izz)*cos(2*phi))*cos(theta) + thetadot*(1.0*Ixx*phidot*cos(theta) + 2.0*psidot*(-Ixx + Iyy*sin(phi)**2 + Izz*cos(phi)**2)*sin(theta)*cos(theta) + 0.25*thetadot*(Iyy - Izz)*(cos(2*phi - theta) - cos(2*phi + theta)))
    B_ctrl[ 0 , 0 ]= sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi)
    B_ctrl[ 0 , 1 ]= 0
    B_ctrl[ 0 , 2 ]= 0
    B_ctrl[ 0 , 3 ]= 0
    B_ctrl[ 1 , 0 ]= -sin(phi)*cos(psi) + sin(psi)*sin(theta)*cos(phi)
    B_ctrl[ 1 , 1 ]= 0
    B_ctrl[ 1 , 2 ]= 0
    B_ctrl[ 1 , 3 ]= 0
    B_ctrl[ 2 , 0 ]= cos(phi)*cos(theta)
    B_ctrl[ 2 , 1 ]= 0
    B_ctrl[ 2 , 2 ]= 0
    B_ctrl[ 2 , 3 ]= 0
    B_ctrl[ 3 , 0 ]= 0
    B_ctrl[ 3 , 1 ]= 0
    B_ctrl[ 3 , 2 ]= 0
    B_ctrl[ 3 , 3 ]= 0
    B_ctrl[ 4 , 0 ]= 0
    B_ctrl[ 4 , 1 ]= 0
    B_ctrl[ 4 , 2 ]= 0
    B_ctrl[ 4 , 3 ]= 0
    B_ctrl[ 5 , 0 ]= 0
    B_ctrl[ 5 , 1 ]= 1
    B_ctrl[ 5 , 2 ]= 0
    B_ctrl[ 5 , 3 ]= 0
    B_ctrl[ 6 , 0 ]= 0
    B_ctrl[ 6 , 1 ]= 0
    B_ctrl[ 6 , 2 ]= 1
    B_ctrl[ 6 , 3 ]= 0
    B_ctrl[ 7 , 0 ]= 0
    B_ctrl[ 7 , 1 ]= 0
    B_ctrl[ 7 , 2 ]= 0
    B_ctrl[ 7 , 3 ]= 1
    invA = np.linalg.inv(A)
    val = invA.dot(B_ctrl)
    inv_val = np.linalg.pinv(val)
    print(inv_val)
    Xddot = invA.dot(B)
    i = 0
    ax = Xddot[i,0]; i+=1
    ay = Xddot[i,0]; i+=1
    az = Xddot[i,0]; i+=1
    theta_lddot = Xddot[i,0]; i += 1
    phi_lddot = Xddot[i,0]; i += 1
    phiddot = Xddot[i,0]; i+=1
    thetaddot = Xddot[i,0]; i+=1
    psiddot = Xddot[i,0]; i+=1


    dXdt = np.array([vx, vy, vz, theta_ldot, phi_ldot, phidot, thetadot, psidot, ax, ay, az, theta_lddot, phi_lddot, phiddot,thetaddot,psiddot])
    return dXdt


parms = parameters()
h = 0.005
t0 = 0
tN = 5
N = int((tN-t0)/h) + 1
t = np.linspace(t0, tN, N)
T = t[N-1]

x0 = 0; y0 = 0; z0 = 1
vx0 = 0; vy0 = 0; vz0 = 0
phi0 = np.deg2rad(0); theta0 = np.deg2rad(0); psi0 = np.deg2rad(0)
phidot0 = 0; thetadot0 = 0; psidot0 = 0
theta_l0 = 0; phi_l0 = 0
thetadot_l0 = 0; phidot_l0 = 0

# t = np.linspace(0, 1, 101)
X0 = np.array([x0, y0, z0, theta_l0, phi_l0, phi0, theta0, psi0, vx0, vy0, vz0, thetadot_l0, phidot_l0, phidot0, thetadot0, psidot0], dtype='float64')
# X = odeint(eom, X0, t, args=all_parms)

mm = len(X0)
X_VAL = [X0[0]]; Y = [X0[1]]; Z = [X0[2]]
PHI = [X0[5]]; THETA = [X0[6]]; PSI = [X0[7]]
VX = [X0[8]]; VY = [X0[9]]; VZ = [X0[10]]
PHIDOT = [X0[13]]; THETADOT = [X0[14]]; PSIDOT = [X0[15]]
THETA_L = [X0[3]]; PHI_L = [X0[4]]
THETA_LDOT = [X0[11]]; PHI_LDOT = [X0[12]]
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
z_ref = np.zeros(N)
v_z = np.zeros(N)
a_z = np.zeros(N)
theta_ldes = np.zeros(N)
phi_ldes = np.zeros(N)
phi = np.zeros(N)
phidot = np.zeros(N)
theta = np.zeros(N)
thetadot = np.zeros(N)
psi = np.zeros(N)
psidot = np.zeros(N)
K_z = parms.K_z
K_phi = parms.K_phi
K_theta = parms.K_theta
K_psi = parms.K_psi
Kp = parms.Kp
Kd = parms.Kd
Kdd = parms.Kdd
Ki = parms.Ki
A = np.array([parms.Ax, parms.Ay, parms.Az])
k = parms.K
b_drag_const = parms.b
l = parms.l
omega = np.array([parms.omega1, parms.omega2, parms.omega3, parms.omega4])
OMEGA = np.zeros((len(t), 4))
OMEGA[0] = omega
X_POS = np.zeros((len(t), 6))
X_ANG = np.zeros((len(t), 6))
X_POS[0, 0] = X0[0]
X_POS[0, 1] = X0[1]
X_POS[0, 2] = X0[2]
X_POS[0, 3] = X0[8]
X_POS[0, 4] = X0[9]
X_POS[0, 5] = X0[10]

X_ANG[0, 0] = X0[5]
X_ANG[0, 1] = X0[6]
X_ANG[0, 2] = X0[7]
X_ANG[0, 3] = X0[13]
X_ANG[0, 4] = X0[14]
X_ANG[0, 5] = X0[15]


for i in range(0, N-1):
    j = 0
    control = Controller(K_z, K_psi, K_theta, K_phi, Kp, Kd, Kdd, Ki, A, k, l, b_drag_const)
    t_temp = np.array([t[i], t[i+1]], dtype='float64')
    # desired_state = control.get_desired_positions(t_temp, desired_traj_values)
    desired_state = np.array([z_ref[i], v_z[i], theta[i], thetadot[i], psi[i], psidot[i], phi[i], phidot[i]])
    all_parms = (parms.m_q,parms.m_l,parms.Ixx,parms.Iyy,parms.Izz,parms.g,parms.l,parms.cable_l,
                 parms.K,parms.b,parms.Ax,parms.Ay,parms.Az,
                 parms.omega1,parms.omega2,parms.omega3,parms.omega4)
    X = odeint(eom, X0, t_temp, args=all_parms)
    ang = np.array([[X[1][3], X[1][9]], [X[1][4], X[1][10]], [X[1][5], X[1][11]]], dtype='float64')
    translation = np.array([[X[1][0], X[1][6]], [X[1][1], X[1][7]], [X[1][2], X[1][8]]], dtype='float64')
    vertical = translation[2]
    torques, theta_temp, psi_temp = control._get_torques(vertical, ang, desired_state)
    omega = control.get_action(desired_state, ang, translation)
    omega_temp = omega.reshape(4,)
    parms.omega1, parms.omega2, parms.omega3, parms.omega4 = omega
    X0 = X[1]
    X_VAL.append(X[1][j]); j+=1;
    Y.append(X[1][j]); j+=1;
    Z.append(X[1][j]); j+=1;
    THETA_L.append(X[1][j]); j+=1;
    PHI_L.append(X[1][j]); j+=1;
    PHI.append(X[1][j]); j+=1;
    THETA.append(X[1][j]); j+=1;
    PSI.append(X[1][j]); j+=1;
    VX.append(X[1][j]); j+=1;
    VY.append(X[1][j]); j+=1;
    VZ.append(X[1][j]); j+=1;
    THETA_LDOT.append(X[1][j]); j+=1;
    PHI_LDOT.append(X[1][j]); j+=1;
    PHIDOT.append(X[1][j]); j+=1;
    THETADOT.append(X[1][j]); j+=1;
    PSIDOT.append(X[1][j]); j+=1;
    X_POS[i+1, 0] = X[1][0]
    X_POS[i+1, 1] = X[1][1]
    X_POS[i+1, 2] = X[1][2]
    X_POS[i+1, 3] = X[1][8]
    X_POS[i+1, 4] = X[1][9]
    X_POS[i+1, 5] = X[1][10]

    X_ANG[i+1, 0] = X[1][5]
    X_ANG[i+1, 1] = X[1][6]
    X_ANG[i+1, 2] = X[1][7]
    X_ANG[i+1, 3] = X[1][13]
    X_ANG[i+1, 4] = X[1][14]
    X_ANG[i+1, 5] = X[1][15]
    # X_pos.append(np.array([X_VAL[i], Y[i], Z[i]]))
    # X_ang.append(np.array([PHI[i], THETA[i], PSI[i]]))
    OMEGA[i+1] = omega_temp

plt.figure(1)
plt.subplot(3,1,1)
plt.plot(t,X_VAL);
plt.plot(t,Y)
plt.plot(t,Z)
plt.ylabel('linear position');
plt.legend(['x_pos', 'y_pos', 'z_pos'])
plt.subplot(3,1,2)
plt.plot(t,VX);
plt.plot(t,VY)
plt.plot(t,VZ)
plt.xlabel('time')
plt.ylabel('linear velocity');
plt.legend(['x_vel', 'y_vel', 'z_vel'])
plt.subplot(3,1,3)
plt.plot(t,PHI);
plt.plot(t,THETA)
plt.plot(t,PSI)
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
ax=plt.figure(4)
plt.subplot(1,1,1)
plt.plot(t,OMEGA);
plt.ylabel('omega');
# plt.subplot(2,1,2)
# plt.plot(t,omega_body_x);
# plt.plot(t,omega_body_y);
# plt.plot(t,omega_body_z);
ax.legend(['rotor1', 'rotor2','rotor3', 'rotor4'])
# plt.ylabel('omega body');
plt.xlabel('time')
#
#
fig = plt.figure(5)
plt.subplot(1, 1, 1)
plt.plot(t, THETA_L)
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
animate(t,X_POS,X_ANG,THETA_L,PHI_L,parms)
