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
K,l,cable_l,b,Ax,Ay,Az = sy.symbols('k l cable_l b Ax Ay Az', real=True)
u1,u2,u3,u4 = sy.symbols('u1 u2 u3 u4', real=True)
K_td, K_tp, K_phid, K_phip, K_thetad, K_thetap, K_psid, K_psip = sy.symbols('K_td K_tp K_phid K_phip K_thetad K_thetap K_psid K_psip', real=True)


# A = sy.Matrix([[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                [0, 0, 0, 0, 0, 0, 0, -(m_l/m_q)*(l/cable_l)*g, 0, 0, 0, -((m_l+m_q)/m_q)*g, 0, 0, 0, 0],
#                [0, 0, 0, 0, 0, 0, (m_l/m_q)*(l/cable_l)*g, 0, 0, 0, -((m_l+m_q)/m_q)*g, 0, 0, 0, 0, 0],
#                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
#                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
#                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g/cable_l, 0, 0, 0, 0, 0],
#                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g/cable_l, 0, 0, 0, 0],
#                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
#                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
#                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
#                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Izz, 0, 0, 0, 0],
#                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
#
# B = sy.Matrix([[0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [-g+(k/(m_l+m_q)), -g+(k/(m_l+m_q)), -g+(k/(m_l+m_q)), -g+(k/(m_l+m_q))],
#                [0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [0, 0, 0, 0],
#                [0, -l*k/(Ixx*theta), 0, l*k/(Ixx*theta)],
#                [(b/Iyy-b/Izz-l*k/(Iyy*phi)), (-b/Izz-b/Iyy), (b/Iyy+b/Izz+l*k/(Iyy*phi)), (-b/Izz-b/Iyy)],
#                [(b/(Izz*phi)-l*k/Izz-l*k/Iyy), (-b/(Izz*phi)-l*k/Izz), (b/(Izz*phi)+l*k/Izz+l*k/Iyy), (-b/(Izz*phi)+l*k/Izz)]])
#
#
# T = (m_l+m_q)*g + (m_l+m_q)*K_td*(vzd - vz) + (m_l+m_q)*K_tp*(zd - z)
# T_phi = Ixx*K_phid*(phidotd - phidot) + Ixx*K_phip*(phid - phi)
# T_theta = Iyy*K_thetad*(thetadotd - thetadot) + Iyy*K_thetap*(thetad - theta)
# T_psi = Izz*K_psid*(psidotd - psidot) + Izz*K_psip*(psid - psi)
#
# omega1_sq = T/(4*k) - T_theta/(2*k*l) - T_psi/(4*b)
# omega2_sq = T/(4*k) - T_theta/(2*k*l) + T_psi/(4*b)
# omega3_sq = T/(4*k) + T_theta/(2*k*l) - T_psi/(4*b)
# omega4_sq = T/(4*k) + T_theta/(2*k*l) + T_psi/(4*b)
#
#
# X = sy.Matrix([x, y, z, vx, vy, vz, theta_l, phi_l, theta_ldot, phi_ldot, phi, theta, psi, phidot, thetadot, psidot])
# u = sy.Matrix([omega1_sq, omega2_sq, omega3_sq, omega4_sq])  # All are square values!!!!
#
# first_term = A*X
# second_term = B*u
#
# xdot = first_term + second_term
# xdot = sy.simplify(xdot)
# print(xdot)
#
# stability_test = A.eigenvals()
# print(stability_test)

# g = 9.81
# l = 0.225  # in m
# cable_l = 0.2
# k = 2.980*1e-6  # this is to be found via calculation
# b = 1.140*1e-2  # this is to be found via calculation
# m_q = 0.468
# m_l = 0.1
# m = m_q + m_l
# Ixx = 4.856*1e-3
# Iyy = 4.856*1e-3
# Izz = 8.801*1e-3
# phi = 0.5
# theta = 0.5
#
# A_np = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
#                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
#                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
#                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
#                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
#                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
#                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
#                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
#                  [0, 0, 0, 0, -(m_l/m_q)*(l/cable_l)*g, 0, -(m/m_q)*g, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                  [0, 0, 0, (m_l/m_q)*(l/cable_l)*g, 0, -(m/m_q)*g, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                  [0, 0, 0, 0, 0, g/cable_l, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                  [0, 0, 0, 0, 0, 0, g/cable_l, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                  [0, 0, 0, 0, 0, 0, 1/Izz, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
#
# B_np = np.array([[0, 0, 0, 0],
#                  [0, 0, 0, 0],
#                  [0, 0, 0, 0],
#                  [0, 0, 0, 0],
#                  [0, 0, 0, 0],
#                  [0, 0, 0, 0],
#                  [0, 0, 0, 0],
#                  [0, 0, 0, 0],
#                  [0, 0, 0, 0],
#                  [0, 0, 0, 0],
#                  [-g+(k/(m_l+m_q)), -g+(k/(m_l+m_q)), -g+(k/(m_l+m_q)), -g+(k/(m_l+m_q))],
#                  [0, 0, 0, 0],
#                  [0, 0, 0, 0],
#                  [0, -l*k/Ixx, 0, l*k/Ixx],
#                  [-l*k/Iyy, 0, l*k/Iyy, 0],
#                  [b/Izz, -b/Izz, b/Izz, -b/Izz]])
#
# controllability_test = control.ctrb(A_np, B_np)
# print(controllability_test)
# ctrb_rank = np.linalg.matrix_rank(controllability_test)
# print(ctrb_rank)
# Q = np.eye(16)
# R = 1e-5*np.eye(4)
# K, P, E = control.lqr(A_np, B_np, Q, R)
# print("Feedback gain values are: {}".format(K))

sin = sy.sin
cos = sy.cos
A = sy.zeros(8,8)
B = sy.zeros(8,1)

# A[ 0 , 0 ]= 1.0*m_l + 1.0*m_q
# A[ 0 , 1 ]= 0
# A[ 0 , 2 ]= 0
# A[ 0 , 3 ]= 1.0*cable_l*m_l*sin(phi_l)*sin(theta_l)
# A[ 0 , 4 ]= -1.0*cable_l*m_l*cos(phi_l)*cos(theta_l)
# A[ 0 , 5 ]= 0
# A[ 0 , 6 ]= 0
# A[ 0 , 7 ]= 0
# A[ 1 , 0 ]= 0
# A[ 1 , 1 ]= 1.0*m_l + 1.0*m_q
# A[ 1 , 2 ]= 0
# A[ 1 , 3 ]= 1.0*cable_l*m_l*cos(theta_l)
# A[ 1 , 4 ]= 0
# A[ 1 , 5 ]= 0
# A[ 1 , 6 ]= 0
# A[ 1 , 7 ]= 0
# A[ 2 , 0 ]= 0
# A[ 2 , 1 ]= 0
# A[ 2 , 2 ]= 1.0*m_l + 1.0*m_q
# A[ 2 , 3 ]= 1.0*cable_l*m_l*sin(theta_l)*cos(phi_l)
# A[ 2 , 4 ]= 1.0*cable_l*m_l*sin(phi_l)*cos(theta_l)
# A[ 2 , 5 ]= 0
# A[ 2 , 6 ]= 0
# A[ 2 , 7 ]= 0
# A[ 3 , 0 ]= 1.0*cable_l*m_l*sin(phi_l)*sin(theta_l)
# A[ 3 , 1 ]= 1.0*cable_l*m_l*cos(theta_l)
# A[ 3 , 2 ]= 1.0*cable_l*m_l*sin(theta_l)*cos(phi_l)
# A[ 3 , 3 ]= 1.0*cable_l**2*m_l
# A[ 3 , 4 ]= 0
# A[ 3 , 5 ]= 0
# A[ 3 , 6 ]= 0
# A[ 3 , 7 ]= 0
# A[ 4 , 0 ]= -1.0*cable_l*m_l*cos(phi_l)*cos(theta_l)
# A[ 4 , 1 ]= 0
# A[ 4 , 2 ]= 1.0*cable_l*m_l*sin(phi_l)*cos(theta_l)
# A[ 4 , 3 ]= 0
# A[ 4 , 4 ]= 1.0*cable_l**2*m_l*cos(theta_l)**2
# A[ 4 , 5 ]= 0
# A[ 4 , 6 ]= 0
# A[ 4 , 7 ]= 0
# A[ 5 , 0 ]= 0
# A[ 5 , 1 ]= 0
# A[ 5 , 2 ]= 0
# A[ 5 , 3 ]= 0
# A[ 5 , 4 ]= 0
# A[ 5 , 5 ]= 1.0*Ixx
# A[ 5 , 6 ]= 0
# A[ 5 , 7 ]= -1.0*Ixx*sin(theta)
# A[ 6 , 0 ]= 0
# A[ 6 , 1 ]= 0
# A[ 6 , 2 ]= 0
# A[ 6 , 3 ]= 0
# A[ 6 , 4 ]= 0
# A[ 6 , 5 ]= 0
# A[ 6 , 6 ]= 1.0*Iyy*cos(phi)**2 + 1.0*Izz*sin(phi)**2
# A[ 6 , 7 ]= 0.25*(Iyy - Izz)*(sin(2*phi - theta) + sin(2*phi + theta))
# A[ 7 , 0 ]= 0
# A[ 7 , 1 ]= 0
# A[ 7 , 2 ]= 0
# A[ 7 , 3 ]= 0
# A[ 7 , 4 ]= 0
# A[ 7 , 5 ]= -1.0*Ixx*sin(theta)
# A[ 7 , 6 ]= 0.25*(Iyy - Izz)*(sin(2*phi - theta) + sin(2*phi + theta))
# A[ 7 , 7 ]= 1.0*Ixx*sin(theta)**2 + 1.0*Iyy*sin(phi)**2*cos(theta)**2 + 1.0*Izz*cos(phi)**2*cos(theta)**2
# B[ 0 ]= -Ax*vx + K*(sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi))*(omega1**2 + omega2**2 + omega3**2 + omega4**2) - 1.0*cable_l*m_l*phi_ldot*(phi_ldot*sin(phi_l)*cos(theta_l) + theta_ldot*sin(theta_l)*cos(phi_l)) - 1.0*cable_l*m_l*theta_ldot*(phi_ldot*sin(theta_l)*cos(phi_l) + theta_ldot*sin(phi_l)*cos(theta_l))
# B[ 1 ]= -Ay*vy - K*(sin(phi)*cos(psi) - sin(psi)*sin(theta)*cos(phi))*(omega1**2 + omega2**2 + omega3**2 + omega4**2) + 1.0*cable_l*m_l*theta_ldot**2*sin(theta_l)
# B[ 2 ]= -Az*vz + K*(omega1**2 + omega2**2 + omega3**2 + omega4**2)*cos(phi)*cos(theta) - 1.0*cable_l*m_l*phi_ldot*(phi_ldot*cos(phi_l)*cos(theta_l) - theta_ldot*sin(phi_l)*sin(theta_l)) + 1.0*cable_l*m_l*theta_ldot*(phi_ldot*sin(phi_l)*sin(theta_l) - theta_ldot*cos(phi_l)*cos(theta_l)) - g*m_l - g*m_q
# B[ 3 ]= -1.0*m_l*(cable_l**2*phi_ldot**2*cos(theta_l) + g*l*cos(phi_l))*sin(theta_l)
# B[ 4 ]= m_l*(2.0*cable_l**2*phi_ldot*theta_ldot*sin(theta_l) - 1.0*g*l*sin(phi_l))*cos(theta_l)
# B[ 5 ]= Ixx*psidot*thetadot*cos(theta) - K*l*(omega2**2 - omega4**2) + psidot*(psidot*(Iyy - Izz)*(sin(2*phi - theta) + sin(2*phi + theta)) - 2.0*thetadot*(-Iyy + Izz)*cos(2*phi))*cos(theta)/4 - thetadot*(psidot*(-Iyy + Izz)*cos(2*phi)*cos(theta) + thetadot*(Iyy - Izz)*sin(2*phi))/2
# B[ 6 ]= -0.5*Ixx*phidot*psidot*cos(theta) - K*l*(omega1**2 - omega3**2) - 1.0*phidot*(psidot*(Iyy - Izz)*cos(2*phi)*cos(theta) - thetadot*(Iyy - Izz)*sin(2*phi)) + 0.125*psidot*thetadot*(Iyy - Izz)*(cos(2*phi - theta) - cos(2*phi + theta)) - psidot*(0.5*Ixx*phidot*cos(theta) + psidot*(-Ixx + Iyy*sin(phi)**2 + Izz*cos(phi)**2)*sin(theta)*cos(theta) + 0.125*thetadot*(Iyy - Izz)*(cos(2*phi - theta) - cos(2*phi + theta)))
# B[ 7 ]= b*(omega1**2 - omega2**2 + omega3**2 - omega4**2) - phidot*(0.5*psidot*(Iyy - Izz)*(sin(2*phi - theta) + sin(2*phi + theta)) - 1.0*thetadot*(-Iyy + Izz)*cos(2*phi))*cos(theta) + thetadot*(1.0*Ixx*phidot*cos(theta) + 2.0*psidot*(-Ixx + Iyy*sin(phi)**2 + Izz*cos(phi)**2)*sin(theta)*cos(theta) + 0.25*thetadot*(Iyy - Izz)*(cos(2*phi - theta) - cos(2*phi + theta)))
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
B[ 0 ]= -1.0*Ax*vx - 1.0*cable_l*m_l*phi_ldot**2*sin(phi_l)*cos(theta_l) - 2.0*cable_l*m_l*phi_ldot*theta_ldot*sin(theta_l)*cos(phi_l) - 1.0*cable_l*m_l*theta_ldot**2*sin(phi_l)*cos(theta_l) + 1.0*u1*sin(phi)*sin(psi) + 1.0*u1*sin(theta)*cos(phi)*cos(psi)
B[ 1 ]= -Ay*vy + 1.0*cable_l*m_l*theta_ldot**2*sin(theta_l) - u1*(sin(phi)*cos(psi) - sin(psi)*sin(theta)*cos(phi))
B[ 2 ]= -1.0*Az*vz - 1.0*cable_l*m_l*phi_ldot**2*cos(phi_l)*cos(theta_l) + 2.0*cable_l*m_l*phi_ldot*theta_ldot*sin(phi_l)*sin(theta_l) - 1.0*cable_l*m_l*theta_ldot**2*cos(phi_l)*cos(theta_l) - 1.0*g*m_l - 1.0*g*m_q + 1.0*u1*cos(phi)*cos(theta)
B[ 3 ]= -1.0*m_l*(cable_l**2*phi_ldot**2*cos(theta_l) + g*l*cos(phi_l))*sin(theta_l)
B[ 4 ]= m_l*(2.0*cable_l**2*phi_ldot*theta_ldot*sin(theta_l) - 1.0*g*l*sin(phi_l))*cos(theta_l)
B[ 5 ]= Ixx*psidot*thetadot*cos(theta) + psidot*(psidot*(Iyy - Izz)*(sin(2*phi - theta) + sin(2*phi + theta)) - 2.0*thetadot*(-Iyy + Izz)*cos(2*phi))*cos(theta)/4 - thetadot*(psidot*(-Iyy + Izz)*cos(2*phi)*cos(theta) + thetadot*(Iyy - Izz)*sin(2*phi))/2 + u2
B[ 6 ]= -1.0*Ixx*phidot*psidot*cos(theta) + Ixx*psidot**2*sin(theta)*cos(theta) - 1.0*Iyy*phidot*psidot*cos(2*phi)*cos(theta) + 1.0*Iyy*phidot*thetadot*sin(2*phi) - Iyy*psidot**2*sin(phi)**2*sin(theta)*cos(theta) + 1.0*Izz*phidot*psidot*cos(2*phi)*cos(theta) - 1.0*Izz*phidot*thetadot*sin(2*phi) - Izz*psidot**2*sin(theta)*cos(phi)**2*cos(theta) + u3
B[ 7 ]= -phidot*(0.5*psidot*(Iyy - Izz)*(sin(2*phi - theta) + sin(2*phi + theta)) - 1.0*thetadot*(-Iyy + Izz)*cos(2*phi))*cos(theta) + thetadot*(1.0*Ixx*phidot*cos(theta) + 2.0*psidot*(-Ixx + Iyy*sin(phi)**2 + Izz*cos(phi)**2)*sin(theta)*cos(theta) + 0.25*thetadot*(Iyy - Izz)*(cos(2*phi - theta) - cos(2*phi + theta))) + u4
Ainv = A.inv()
xddot = Ainv.dot(B)

XDOT = sy.Matrix([vx, vy, vz, theta_ldot, phi_ldot, phidot, thetadot, psidot, xddot[0], xddot[1], xddot[2], xddot[3], xddot[4], xddot[5], xddot[6], xddot[7]])
X = sy.Matrix([x, y, z, theta_l, phi_l, phi, theta, psi, vx, vy, vz, theta_ldot, phi_ldot, phidot, thetadot, psidot])
u = sy.Matrix([u1, u2, u3, u4])

A_lin = XDOT.jacobian(X)
B_lin = XDOT.jacobian(u)

A_lin = A_lin.subs([(Ax, 0), (Ay, 0), (Az, 0), (phi**3, 0), (sin(2*phi + theta), 2*phi+theta), (sin(2*phi-theta), 2*phi-theta), (cos(2*phi + theta), 1), (cos(2*phi-theta), 1), (sin(phi), phi), (sin(theta), theta), (cos(phi), 1), (cos(theta), 1), (psi, 0), (theta_ldot, 0), (phi_ldot, 0), (sin(theta_l), theta_l), (sin(phi_l), phi_l), (cos(theta_l), 1), (cos(phi_l), 1), (phi_l*theta_l, 0), (phi_l**2, 0), (phi_l**4, 0), (theta_l**4, 0), (theta_l**2, 0), (phi*theta, 0), (phi**2, 0), (theta**2, 0), (phi**4, 0), (theta**4, 0), (phidot, 0), (thetadot, 0), (psidot, 0)])
B_lin = B_lin.subs([(Ax, 0), (Ay, 0), (Az, 0), (phi**3, 0), (sin(2*phi + theta), 2*phi+theta), (sin(2*phi-theta), 2*phi-theta), (cos(2*phi + theta), 1), (cos(2*phi-theta), 1), (sin(phi), phi), (sin(theta), theta), (cos(phi), 1), (cos(theta), 1), (psi, 0), (theta_ldot, 0), (phi_ldot, 0), (sin(theta_l), theta_l), (sin(phi_l), phi_l), (cos(theta_l), 1), (cos(phi_l), 1), (phi_l*theta_l, 0), (phi_l**2, 0), (phi_l**4, 0), (theta_l**4, 0), (theta_l**2, 0), (phi*theta, 0), (phi**2, 0), (theta**2, 0), (phi**4, 0), (theta**4, 0), (phidot, 0), (thetadot, 0), (psidot, 0)])


A_lin = A_lin.subs([((2*phi + theta)**2, 0), ((2*phi - theta)**2, 0), ((2*phi + theta)*(2*phi - theta), 0), ((2*phi - theta)*(2*phi + theta), 0)])
B_lin = B_lin.subs([((2*phi + theta)**2, 0), ((2*phi - theta)**2, 0), ((2*phi + theta)*(2*phi - theta), 0), ((2*phi - theta)*(2*phi + theta), 0)])

[mm,nn]=np.shape(A_lin)
for ii in range(0,mm):
    for jj in range(0,nn):
        print('A[',ii,',',jj,']=',A_lin[ii,jj])


[mm, nn] = np.shape(B_lin)
for ii in range(0,mm):
    for jj in range(0,nn):
        print('B[',ii,',',jj,']=',B_lin[ii,jj])
