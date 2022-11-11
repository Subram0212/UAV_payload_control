import sympy as sy
import numpy as np
import control


x,y,z, xd, yd, zd  = sy.symbols('x y z xd yd zd', real=True)
vx,vy,vz, vxd, vyd, vzd  = sy.symbols('vx vy vz vxd vyd vzd', real=True)
ax,ay,az  = sy.symbols('ax ay az', real=True)
theta_l,phi_l,theta_ldot,phi_ldot  = sy.symbols('theta_l phi_l theta_ldot phi_ldot', real=True)
phi,theta,psi, phid, thetad, psid  = sy.symbols('phi theta psi phid thetad psid', real=True)
phidot,thetadot,psidot, phidotd, thetadotd, psidotd  = sy.symbols('phidot thetadot psidot phidotd thetadotd psidotd', real=True)
phiddot,thetaddot,psiddot = sy.symbols('phiddot thetaddot psiddot', real=True)
m_l,m_q,g,Ixx,Iyy,Izz = sy.symbols('m_l m_q g Ixx Iyy Izz', real=True)
k,l,cable_l,b,Ax,Ay,Az = sy.symbols('k l cable_l b Ax Ay Az', real=True)
omega1,omega2,omega3,omega4 = sy.symbols('omega1 omega2 omega3 omega4', real=True)
K_td, K_tp, K_phid, K_phip, K_thetad, K_thetap, K_psid, K_psip = sy.symbols('K_td K_tp K_phid K_phip K_thetad K_thetap K_psid K_psip', real=True)


A = sy.Matrix([[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, -(m_l/m_q)*(l/cable_l)*g, 0, 0, 0, -((m_l+m_q)/m_q)*g, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, (m_l/m_q)*(l/cable_l)*g, 0, 0, 0, -((m_l+m_q)/m_q)*g, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g/cable_l, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g/cable_l, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Izz, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

B = sy.Matrix([[0, 0, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 0],
               [-g+(k/(m_l+m_q)), -g+(k/(m_l+m_q)), -g+(k/(m_l+m_q)), -g+(k/(m_l+m_q))],
               [0, 0, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 0],
               [0, -l*k/(Ixx*theta), 0, l*k/(Ixx*theta)],
               [(b/Iyy-b/Izz-l*k/(Iyy*phi)), (-b/Izz-b/Iyy), (b/Iyy+b/Izz+l*k/(Iyy*phi)), (-b/Izz-b/Iyy)],
               [(b/(Izz*phi)-l*k/Izz-l*k/Iyy), (-b/(Izz*phi)-l*k/Izz), (b/(Izz*phi)+l*k/Izz+l*k/Iyy), (-b/(Izz*phi)+l*k/Izz)]])


T = (m_l+m_q)*g + (m_l+m_q)*K_td*(vzd - vz) + (m_l+m_q)*K_tp*(zd - z)
T_phi = Ixx*K_phid*(phidotd - phidot) + Ixx*K_phip*(phid - phi)
T_theta = Iyy*K_thetad*(thetadotd - thetadot) + Iyy*K_thetap*(thetad - theta)
T_psi = Izz*K_psid*(psidotd - psidot) + Izz*K_psip*(psid - psi)

omega1_sq = T/(4*k) - T_theta/(2*k*l) - T_psi/(4*b)
omega2_sq = T/(4*k) - T_theta/(2*k*l) + T_psi/(4*b)
omega3_sq = T/(4*k) + T_theta/(2*k*l) - T_psi/(4*b)
omega4_sq = T/(4*k) + T_theta/(2*k*l) + T_psi/(4*b)


X = sy.Matrix([x, y, z, vx, vy, vz, theta_l, phi_l, theta_ldot, phi_ldot, phi, theta, psi, phidot, thetadot, psidot])
u = sy.Matrix([omega1_sq, omega2_sq, omega3_sq, omega4_sq])  # All are square values!!!!

first_term = A*X
second_term = B*u

xdot = first_term + second_term
xdot = sy.simplify(xdot)
print(xdot)

stability_test = A.eigenvals()
print(stability_test)

g = 9.81
l = 0.225  # in m
cable_l = 0.2
k = 2.980*1e-6  # this is to be found via calculation
b = 1.140*1e-2  # this is to be found via calculation
m_q = 0.468
m_l = 0.1
m = m_q + m_l
Ixx = 4.856*1e-3
Iyy = 4.856*1e-3
Izz = 8.801*1e-3
phi = 0.5
theta = 0.5

A_np = np.array([[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, -(m_l/m_q)*(l/cable_l)*g, 0, 0, 0, -(m/m_q)*g, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, (m_l/m_q)*(l/cable_l)*g, 0, 0, 0, -(m/m_q)*g, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g/cable_l, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g/cable_l, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Izz, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

B_np = np.array([[0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [-g+(k/(m_l+m_q)), -g+(k/(m_l+m_q)), -g+(k/(m_l+m_q)), -g+(k/(m_l+m_q))],
                 [0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, -l*k/Ixx, 0, l*k/Ixx],
                 [-l*k/Iyy, 0, l*k/Iyy, 0],
                 [b/Izz, -b/Izz, b/Izz, -b/Izz]])

controllability_test = control.ctrb(A_np, B_np)
print(controllability_test)
ctrb_rank = np.linalg.matrix_rank(controllability_test)
print(ctrb_rank)
Q = np.eye(16)
R = 1e-5*np.eye(4)
K, P, E = control.lqr(A_np, B_np, Q, R)
print("Feedback gain values are: {}".format(K))

