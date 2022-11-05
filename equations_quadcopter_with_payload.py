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
m,g,Ixx,Iyy,Izz = sy.symbols('m g Ixx Iyy Izz', real=True)
K,l,b,Ax,Ay,Az = sy.symbols('K l b Ax Ay Az', real=True)
omega1,omega2,omega3,omega4 = sy.symbols('omega1 omega2 omega3 omega4', real=True)


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


#2) angular velocity and energy
om_b = phidot*i + R_x.transpose()*(thetadot*j) + R_x.transpose()*R_y.transpose()*(psidot*k)
I = sy.Matrix( [
    [Ixx, 0, 0],
    [0, Iyy, 0],
    [0,  0, Izz]
]) #body frame inertia
v = sy.Matrix([vx, vy, vz])

T = 0.5*m*v.dot(v) + 0.5*om_b.dot(I*om_b)
V = m*g*z
L = T - V
# print('KE=',T);
# print('PE=',V);
# print('TE= PE+KE');

#external force
Thrust = sy.Matrix([0, 0,K*(omega1**2+omega2**2+omega3**2+omega4**2)])
Drag = sy.Matrix([Ax*vx, Ay*vy, Az*vz])
F_ext = R*Thrust-Drag
tau_phi = K*l*(omega4**2 - omega2**2)
tau_theta = K*l*(omega3**2 - omega1**2)
tau_psi = b*(omega1**2-omega2**2+omega3**2-omega4**2)
tau_ext = sy.Matrix([tau_phi,tau_theta,tau_psi])
# T_ext = np.concatenate(F_ext,tau_ext)
# print(T_ext)
T_ext = F_ext.col_join(tau_ext)

# T_ext = [F_ext; tau_ext];


# #3) Derive equations
q = sy.Matrix([x,y,z,phi,theta,psi])
qdot = sy.Matrix([vx,vy,vz,phidot,thetadot,psidot])
qddot = sy.Matrix([ax, ay, az, phiddot, thetaddot, psiddot])
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
    EOM.append(ddt_dLdqdot[ii] - dLdq[ii]-T_ext[ii])

ndof = len(q)
EOM = sy.Matrix([EOM[0],EOM[1],EOM[2],EOM[3],EOM[4],EOM[5]])
# print(EOM)
# print(len(EOM))
# print(type(qddot))
# print(type(EOM))

A = EOM.jacobian(qddot)
B = []
for ii in range(0,ndof):
    B_temp = -EOM[ii].subs([(ax,0), (ay,0), (az,0), (phiddot,0), (thetaddot,0), (psiddot,0)])
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
