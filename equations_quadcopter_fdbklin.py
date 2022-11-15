import sympy as sy
import numpy as np

def cos(theta):
    return sy.cos(theta)

def sin(theta):
    return sy.sin(theta)

x,y,z, xd, yd, zd  = sy.symbols('x y z xd yd zd', real=True)
vx,vy,vz, vxd, vyd, vzd  = sy.symbols('vx vy vz vxd vyd vzd', real=True)
ax,ay,az,axd,ayd,azd  = sy.symbols('ax ay az axd ayd azd', real=True)
xl, yl, zl = sy.symbols('xl yl zl', real=True)
vxl, vyl, vzl = sy.symbols('vxl vyl vzl', real=True)
phi,theta,psi,psid  = sy.symbols('phi theta psi psid', real=True)
phidot,thetadot,psidot,psidesdot = sy.symbols('phidot thetadot psidot psidesdot', real=True)
phiddot,thetaddot,psiddot = sy.symbols('phiddot thetaddot psiddot', real=True)
theta_l, phi_l = sy.symbols('theta_l phi_l', real=True)
theta_ldot, phi_ldot = sy.symbols('theta_ldot phi_ldot', real=True)
theta_lddot, phi_lddot = sy.symbols('theta_lddot phi_lddot', real=True)
m_q,m_l,g,Ixx,Iyy,Izz = sy.symbols('m_q m_l g Ixx Iyy Izz', real=True)
K,l,b, cable_l, Ax,Ay,Az = sy.symbols('K l b cable_l Ax Ay Az', real=True)
u1,u2,u3,u4 = sy.symbols('u1 u2 u3 u4', real=True)
kpx, kpy, kpz, kppsi = sy.symbols('kpx kpy kpz kppsi', real=True)
kdx, kdy, kdz, kdpsi = sy.symbols('kdx kdy kdz kdpsi', real=True)


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

#external force
Fz = u1
Thrust = sy.Matrix([0, 0, Fz])
Drag = sy.Matrix([Ax*vx, Ay*vy, Az*vz])
F_ext = R*Thrust-Drag
tau_phi = u2
tau_theta = u3
tau_psi = u4
tau_ext = sy.Matrix([0, 0, tau_phi,tau_theta,tau_psi])
# T_ext = np.concatenate(F_ext,tau_ext)
# print(T_ext)
T_ext = F_ext.col_join(tau_ext)
# T_ext = sy.Matrix(T_ext)

# T_ext = [F_ext; tau_ext];


# #3) Derive equations
qdot = sy.Matrix([vx, vy, vz, theta_ldot, phi_ldot, phidot, thetadot, psidot])
qddot = sy.Matrix([ax, ay, az, theta_lddot, phi_lddot, phiddot, thetaddot, psiddot])
dLdqdot = []
ddt_dLdqdot = []
dLdq = []
EOM = []
mm = len(qddot)
for ii in range(0,mm):
    dLdqdot.append(sy.diff(L,qdot[ii]))
    tmp = 0
    for jj in range(0,mm):
        tmp += sy.diff(dLdqdot[ii],q[jj])*qdot[jj]+ sy.diff(dLdqdot[ii],qdot[jj])*qddot[jj]
    ddt_dLdqdot.append(tmp)
    dLdq.append(sy.diff(L,q[ii]))
    EOM.append(ddt_dLdqdot[ii] - dLdq[ii] - T_ext[ii])

ndof = len(q)
EOM = sy.Matrix([EOM[0],EOM[1],EOM[2],EOM[3],EOM[4],EOM[5],EOM[6],EOM[7]])
# print(EOM)
# print(len(EOM))
# print(type(qddot))
# print(type(EOM))

A = EOM.jacobian(qddot)
B = []
for ii in range(0,ndof):
    B_temp = -EOM[ii].subs([(ax,0), (ay,0), (az,0), (theta_lddot, 0), (phi_lddot, 0), (phiddot,0), (thetaddot,0), (psiddot,0)])
    B.append(B_temp)
B_sym = sy.Matrix(B)

[mm,nn]=np.shape(A)
for ii in range(0,mm):
    for jj in range(0,nn):
        print('A[',ii,',',jj,']=',sy.simplify(A[ii,jj]))

mm = len(B)
for ii in range(0,mm):
    print('B[',ii,']=',sy.simplify(B[ii]))

D = []
for ii in range(0,ndof):
    D_temp = -B_sym[ii].subs([(u1,0), (u2,0), (u3,0), (u4, 0)])
    D.append(D_temp)

D_sym = sy.Matrix(D)
Bu_mat = []
for ii in range(0,ndof):
    Bu_mat_temp = B_sym[ii] - D_sym[ii]
    Bu_mat.append(Bu_mat_temp)

u = sy.Matrix([u1, u2, u3, u4])
Bu_mat_sym = sy.Matrix(Bu_mat)

B_control = Bu_mat_sym.jacobian(u)
[mm,nn]=np.shape(B_control)
for ii in range(0,mm):
    for jj in range(0,nn):
        print('B_ctrl[',ii,',',jj,']=',sy.simplify(B_control[ii,jj]))

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

# u = sy.Matrix([u1, u2, u3, u4])
# B_control = sy.Matrix([[0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [K, K, K, K],
#                [0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [0, -K*l, 0, K*l],
#                [-K*l, 0, K*l, 0],
#                [b, -b, b, -b]])
# D_mat = sy.Matrix(B)
#
# T_ext = B_control*u
# control_ip = T_ext + D_mat
Ainv = A.inv()
# xddot = Ainv.dot(control_ip)  # Have a check on this kind of manipulation in case solution isn't obtained.


# print('xddot',sy.simplify(xddot[0]))
# print('yddot',sy.simplify(xddot[1]))
# print('zddot',sy.simplify(xddot[2]))
# print('thetalddot',sy.simplify(xddot[3]))
# print('philddot',sy.simplify(xddot[4]))
# print('phiddot',sy.simplify(xddot[5]))
# print('thetaddot',sy.simplify(xddot[6]))
# print('psiddot',sy.simplify(xddot[7]))

selection_matrix = sy.diag(1, 1, 1, 0, 0, 0, 0, 1)
selected_A = selection_matrix*A
selected_A = sy.Matrix([[selected_A[0], selected_A[1], selected_A[2], selected_A[7]],
                        [selected_A[8], selected_A[9], selected_A[10], selected_A[15]],
                        [selected_A[16], selected_A[17], selected_A[18], selected_A[23]],
                        [selected_A[56], selected_A[57], selected_A[58], selected_A[63]]])
selected_B = selection_matrix*B_control
selected_B = sy.Matrix([[selected_B[0], selected_B[1], selected_B[2], selected_B[3]],
                        [selected_B[4], selected_B[5], selected_B[6], selected_B[7]],
                        [selected_B[8], selected_B[9], selected_B[10], selected_B[11]],
                        [selected_B[28], selected_B[29], selected_B[30], selected_B[31]]])
selected_D = selection_matrix*D_sym
selected_D = sy.Matrix([selected_D[0], selected_D[1], selected_D[2], selected_D[7]])

v_x = axd + kdx*(vxd - vx) + kpx*(xd - x)
v_y = ayd + kdy*(vyd - vy) + kpy*(yd - y)
v_z = azd + kdz*(vzd - vz) + kpz*(zd - z)
v_thetal = 0
v_phil = 0
v_phi = 0
v_theta = 0
v_psi = kdpsi*(psidesdot - psidot) + kppsi*(psid - psi)

val = Ainv*B_control
inverse_val = val.pinv(method='ED')

V = sy.Matrix([v_x, v_y, v_z, v_psi])
selected_Ainv = selected_A.inv()
selected_Binv = selected_B.inv()
val = selected_Ainv*selected_B
inverse_val = val.inv()
AinvD = Ainv*selected_D

U = selection_matrix*inverse_val*V + selection_matrix*inverse_val*AinvD
