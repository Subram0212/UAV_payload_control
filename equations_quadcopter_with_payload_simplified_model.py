import sympy as sy
import numpy as np

def cos(theta):
    return sy.cos(theta)

def sin(theta):
    return sy.sin(theta)

x,y,z  = sy.symbols('x y z', real=True)
vx,vy,vz  = sy.symbols('vx vy vz', real=True)
ax,ay,az  = sy.symbols('ax ay az', real=True)
xl, yl, zl = sy.symbols('xl yl zl', real=True)
vxl, vyl, vzl = sy.symbols('vxl vyl vzl', real=True)
phi,theta,psi  = sy.symbols('phi theta psi', real=True)
phidot,thetadot,psidot = sy.symbols('phidot thetadot psidot', real=True)
phiddot,thetaddot,psiddot = sy.symbols('phiddot thetaddot psiddot', real=True)
theta_l, phi_l = sy.symbols('theta_l phi_l', real=True)
theta_ldot, phi_ldot = sy.symbols('theta_ldot phi_ldot', real=True)
theta_lddot, phi_lddot = sy.symbols('theta_lddot phi_lddot', real=True)
m_q,m_l,g,Ixx,Iyy,Izz = sy.symbols('m_q m_l g Ixx Iyy Izz', real=True)
K,l,b, cable_l, Ax,Ay,Az = sy.symbols('K l b cable_l Ax Ay Az', real=True)
omega1,omega2,omega3,omega4 = sy.symbols('omega1 omega2 omega3 omega4', real=True)
u1, tau_theta, theta_d = sy.symbols('u1 tau_theta theta_d', real=True)


# %%%%%%% unit vectors %%%%%%%
i = sy.Matrix([1, 0, 0])
j = sy.Matrix([0, 1, 0])
k = sy.Matrix([0, 0, 1])

# 1) position and angles
R_x = sy.Matrix([
    [1,            0,         0],
    [0,     cos(phi), -sin(phi)],
    [0,     sin(phi),  cos(phi)]

])

R_y = sy.Matrix([
    [cos(theta),  0, sin(theta)],
    [0,           1,          0],
    [-sin(theta),  0, cos(theta)]
])

R_z = sy.Matrix( [
    [cos(psi), -sin(psi), 0],
    [sin(psi),  cos(psi), 0],
    [0,            0,         1]
])

R = R_z*R_y*R_x

R_y_phil = sy.Matrix([
    [cos(phi_l),  0, sin(phi_l)],
    [0,           1,          0],
    [-sin(phi_l),  0, cos(phi_l)]
])

R_x_thetal = sy.Matrix([
    [1,            0,         0],
    [0,     cos(theta_l), -sin(theta_l)],
    [0,     sin(theta_l),  cos(theta_l)]

])

# Calculate the position of the payload
X = sy.Matrix([x, y, z])
X_l = X + R_y_phil*R_x_thetal*sy.Matrix([0, 0, -cable_l])

#2) angular velocity and energy
om_b = phidot*i + R_x.transpose()*(thetadot*j) + R_x.transpose()*R_y.transpose()*(psidot*k)
I = sy.Matrix( [
    [Ixx, 0, 0],
    [0, Iyy, 0],
    [0,  0, Izz]
])

v = sy.Matrix([vx, vy, vz])
q = sy.Matrix([x, y, z, theta_l, phi_l, phi, theta, psi])
qdot = sy.Matrix([vx, vy, vz, theta_ldot, phi_ldot, phidot, thetadot, psidot])
J_11 = sy.eye(3)
J_12 = sy.Matrix([
    [cable_l*sin(phi_l)*sin(theta_l), -cable_l*cos(theta_l)*cos(phi_l)],
    [cable_l*cos(theta_l), 0],
    [cable_l*cos(phi_l)*sin(theta_l), cable_l*cos(theta_l)*sin(phi_l)]
])
J_1 = sy.BlockMatrix([J_11, J_12])
vel = sy.Matrix([vx, vy, vz, theta_ldot, phi_ldot])
V_l = J_1*vel
J_2 = sy.Matrix([
    [1, 0, -sin(theta)],
    [0, cos(phi), cos(theta)*sin(phi)],
    [0, -sin(phi), cos(theta)*cos(phi)]
])
M_Q = sy.diag(m_q, m_q, m_q)
M_l = sy.diag(m_l, m_l, m_l)
# M11 = M_Q+J_11.transpose()*M_l*J_11
# M12 = J_11.transpose()*M_l*J_12
# M13 = sy.zeros(3,3)
# M21 = J_12.transpose()*M_l*J_11
# M22 = J_12.transpose()*M_l*J_12
# M23 = sy.zeros(3,3)
# M31 = sy.zeros(3,3)
# M32 = sy.zeros(3,2)
# M33 = J_12.transpose()*I*J_2
# M = sy.BlockMatrix([
#     [M11, M12, M13],
#     [M21, M22, M23],
#     [M31, M32, M33]
# ])
M11 = M_Q
M12 = sy.zeros(3,3)
M13 = sy.zeros(3,3)
M21 = sy.zeros(3,3)
M22 = M_l
M23 = sy.zeros(3,3)
M31 = sy.zeros(3,3)
M32 = sy.zeros(3,3)
M33 = I
M = sy.Matrix([
    [M_Q, sy.zeros(3,3), sy.zeros(3,3)],
    [sy.zeros(3,3), M_l, sy.zeros(3,3)],
    [sy.zeros(3,3), sy.zeros(3,3), I]
])
J = sy.Matrix([
    [sy.eye(3), sy.zeros(3,2), sy.zeros(3,3)],
    [J_11, J_12, sy.zeros(3,3)],
    [sy.zeros(3,3), sy.zeros(3,2), J_2]
])
T = 0.5*qdot.transpose()*(J.transpose()*M*J)*qdot
V = m_q*g*z + m_l*g*(z - l*cos(theta_l)*cos(phi_l))
L = T[0] - V
# print('KE=',T);
# print('PE=',V);
# print('TE= PE+KE');

'''A simplified model of the Quadcopter with payload control - 2D model in X-Z plane'''
om_b_simple = thetaddot

#external force
Fz = K*(omega1**2+omega2**2+omega3**2+omega4**2)
Thrust = sy.Matrix([0, 0, Fz])
Drag = sy.Matrix([Ax*vx, Ay*vy, Az*vz])
F_ext = R*Thrust-Drag
tau_phi = K*l*(omega4**2 - omega2**2)
# tau_theta = K*l*(omega3**2 - omega1**2)
tau_psi = b*(omega1**2-omega2**2+omega3**2-omega4**2)
tau_ext = sy.Matrix([0, 0, tau_phi,tau_theta,tau_psi])
# T_ext = np.concatenate(F_ext,tau_ext)
# print(T_ext)
T_ext = F_ext.col_join(tau_ext)

# T_ext = [F_ext; tau_ext];


# #3) Derive equations
q_smplfd = sy.Matrix([x, z, phi_l, theta])
qdot = sy.Matrix([vx, vz, phi_ldot, thetadot])
qddot = sy.Matrix([ax, az, phi_lddot, thetaddot])
dLdqdot = []
ddt_dLdqdot = []
dLdq = []
# EOM = []
mm = len(qddot)
# for ii in range(0,mm):
#     dLdqdot.append(sy.diff(L,qdot[ii]))
#     tmp = 0
#     for jj in range(0,mm):
#         tmp += sy.diff(dLdqdot[ii],q[jj])*qdot[jj]+ sy.diff(dLdqdot[ii],qdot[jj])*qddot[jj]
#     ddt_dLdqdot.append(tmp)
#     dLdq.append(sy.diff(L,q[ii]))
#     EOM.append(ddt_dLdqdot[ii] - dLdq[ii]-T_ext[ii])

ndof = len(q_smplfd)
EOM1 = (m_q + m_l)*(az + g) + m_l*cable_l*(cos(phi_l)*(phi_ldot**2) - sin(phi_l)*phi_lddot) - u1*cos(theta)
EOM2 = (m_q + m_l)*ax + m_l*cable_l*(sin(phi_l)*(phi_ldot**2) + cos(phi_l)*phi_lddot) -u1*sin(theta)
EOM3 = m_l*cable_l**2*phi_lddot + m_l*cable_l*sin(phi_l)*az - m_l*cable_l*cos(phi_l)*ax + m_l*g*cable_l*sin(phi_l)
EOM4 = Iyy*om_b_simple - tau_theta

EOM = sy.Matrix([EOM1,EOM2,EOM3,EOM4])
# print(EOM)
# print(len(EOM))
# print(type(qddot))
# print(type(EOM))

A = EOM.jacobian(qddot)
B = []
for ii in range(0,ndof):
    B_temp = -EOM[ii].subs([(ax,0), (az,0), (phi_lddot, 0), (thetaddot,0)])
    B.append(B_temp)

[mm,nn]=np.shape(A)
for ii in range(0,mm):
    for jj in range(0,nn):
        print('A[',ii,',',jj,']=',sy.simplify(A[ii,jj]))

mm = len(B)
for ii in range(0,mm):
    print('B[',ii,']=',sy.simplify(B[ii]))


#world frame velocity
angdot = sy.Matrix([phidot, thetadot, psidot])

om  = psidot*k  + R_z*(thetadot*j) + R_z*R_y*(phidot*i)
R_we = om.jacobian(angdot)
[mm,nn]=np.shape(R_we)
for ii in range(0,mm):
    for jj in range(0,nn):
        print('R_we[',ii,',',jj,']=',sy.simplify(R_we[ii,jj]))

R_be = om_b.jacobian(angdot)
[mm,nn]=np.shape(R_be)
for ii in range(0,mm):
    for jj in range(0,nn):
        print('R_be[',ii,',',jj,']=',sy.simplify(R_be[ii,jj]))

Ainv = A.inv()
xddot = Ainv.dot(B)
print('xddot',sy.simplify(xddot[0]))
print('zddot',sy.simplify(xddot[1]))
print('philddot = ',sy.simplify(xddot[2]))
print('thetaddot',sy.simplify(xddot[3]))

x_ddot = xddot[0]
zddot = xddot[1]

xddot_approx = x_ddot.subs([(cos(theta), 1), (sin(theta), theta_d)])
zddot_approx = zddot.subs([(cos(theta), 1), (sin(theta), theta_d)])
print("xddot_app = ",xddot_approx)
print("zddot_app = ",zddot_approx)
# # print('philddot',sy.simplify(xddot[4]))
# # print('phiddot',sy.simplify(xddot[5]))
# # print('thetaddot',sy.simplify(xddot[6]))
# # print('psiddot',sy.simplify(xddot[7]))
#
# x_ddot = xddot[0]
# yddot = xddot[1]
# zddot = xddot[2]
# thetalddot = xddot[3]
# philddot = xddot[4]
# phiddot = xddot[5]
# thetaddot = xddot[6]
# psiddot = xddot[7]
#
# xddot_lin = x_ddot.subs([(sin(phi), phi), (sin(theta), theta), (cos(phi), 1), (cos(theta), 1), (psi, 0), (theta_ldot, 0), (phi_ldot, 0), (sin(theta_l), theta_l), (sin(phi_l), phi_l), (cos(theta_l), 1), (cos(phi_l), 1), (phi_l*theta_l, 0), (phi_l**2, 0), (phi_l**4, 0), (theta_l**4, 0), (theta_l**2, 0), (phi*theta, 0), (phi**2, 0), (theta**2, 0), (phi**4, 0), (theta**4, 0)])
# yddot_lin = yddot.subs([(sin(phi), phi), (sin(theta), theta), (cos(phi), 1), (cos(theta), 1), (psi, 0), (theta_ldot, 0), (phi_ldot, 0), (sin(theta_l), theta_l), (sin(phi_l), phi_l), (cos(theta_l), 1), (cos(phi_l), 1), (phi_l*theta_l, 0), (phi_l**2, 0), (phi_l**4, 0), (theta_l**4, 0), (theta_l**2, 0), (phi*theta, 0), (phi**2, 0), (theta**2, 0), (phi**4, 0), (theta**4, 0)])
# zddot_lin = zddot.subs([(sin(phi), phi), (sin(theta), theta), (cos(phi), 1), (cos(theta), 1), (psi, 0), (theta_ldot, 0), (phi_ldot, 0), (sin(theta_l), theta_l), (sin(phi_l), phi_l), (cos(theta_l), 1), (cos(phi_l), 1), (phi_l*theta_l, 0), (phi_l**2, 0), (phi_l**4, 0), (theta_l**4, 0), (theta_l**2, 0), (phi*theta, 0), (phi**2, 0), (theta**2, 0), (phi**4, 0), (theta**4, 0)])
# thetalddot_lin = thetalddot.subs([(sin(phi), phi), (sin(theta), theta), (cos(phi), 1), (cos(theta), 1), (psi, 0), (theta_ldot, 0), (phi_ldot, 0), (sin(theta_l), theta_l), (sin(phi_l), phi_l), (cos(theta_l), 1), (cos(phi_l), 1), (phi_l*theta_l, 0), (phi_l**2, 0), (phi_l**4, 0), (theta_l**4, 0), (theta_l**2, 0), (phi*theta, 0), (phi**2, 0), (theta**2, 0), (phi**4, 0), (theta**4, 0)])
# philddot_lin = philddot.subs([(sin(phi), phi), (sin(theta), theta), (cos(phi), 1), (cos(theta), 1), (psi, 0), (theta_ldot, 0), (phi_ldot, 0), (sin(theta_l), theta_l), (sin(phi_l), phi_l), (cos(theta_l), 1), (cos(phi_l), 1), (phi_l*theta_l, 0), (phi_l**2, 0), (phi_l**4, 0), (theta_l**4, 0), (theta_l**2, 0), (phi*theta, 0), (phi**2, 0), (theta**2, 0), (phi**4, 0), (theta**4, 0)])
# phiddot_lin = phiddot.subs([(sin(phi), phi), (sin(theta), theta), (cos(phi), 1), (cos(theta), 1), (psi, 0), (theta_ldot, 0), (phi_ldot, 0), (sin(theta_l), theta_l), (sin(phi_l), phi_l), (cos(theta_l), 1), (cos(phi_l), 1), (phidot, 0), (thetadot, 0), (psidot, 0), (phi_l*theta_l, 0), (phi_l**2, 0), (phi_l**4, 0), (theta_l**4, 0), (theta_l**2, 0), (phi*theta, 0), (phi**2, 0), (theta**2, 0), (phi**4, 0), (theta**4, 0)])
# thetaddot_lin = thetaddot.subs([(sin(phi), phi), (sin(theta), theta), (cos(phi), 1), (cos(theta), 1), (psi, 0), (theta_ldot, 0), (phi_ldot, 0), (sin(theta_l), theta_l), (sin(phi_l), phi_l), (cos(theta_l), 1), (cos(phi_l), 1), (phidot, 0), (thetadot, 0), (psidot, 0), (phi_l*theta_l, 0), (phi_l**2, 0), (phi_l**4, 0), (theta_l**4, 0), (theta_l**2, 0), (phi*theta, 0), (phi**2, 0), (theta**2, 0), (phi**4, 0), (theta**4, 0)])
# psiddot_lin = psiddot.subs([(sin(phi), phi), (sin(theta), theta), (cos(phi), 1), (cos(theta), 1), (psi, 0), (theta_ldot, 0), (phi_ldot, 0), (sin(theta_l), theta_l), (sin(phi_l), phi_l), (cos(theta_l), 1), (cos(phi_l), 1), (phidot, 0), (thetadot, 0), (psidot, 0), (phi_l*theta_l, 0), (phi_l**2, 0), (phi_l**4, 0), (theta_l**4, 0), (theta_l**2, 0), (phi*theta, 0), (phi**2, 0), (theta**2, 0), (phi**4, 0), (theta**4, 0)])
#
# print('xddot_lin',sy.simplify(xddot_lin))
# print('yddot_lin',sy.simplify(yddot_lin))
# print('zddot_lin',sy.simplify(zddot_lin))
# print('thetalddot_lin',sy.simplify(thetalddot_lin))
# print('philddot_lin',sy.simplify(philddot_lin))
# print('phiddot_lin',sy.simplify(phiddot_lin))
# print('thetaddot_lin',sy.simplify(thetaddot_lin))
# print('psiddot_lin',sy.simplify(psiddot_lin))
#
# print(EOM[0])
# print('\n')
# print(EOM[1])
# print('\n')
# print(EOM[2])
# print('\n')
# print(EOM[3])
# print('\n')
# print(EOM[4])
# print('\n')
# print(EOM[5])
# print('\n')
# print(EOM[6])
# print('\n')
# print(EOM[7])
# print('\n')
