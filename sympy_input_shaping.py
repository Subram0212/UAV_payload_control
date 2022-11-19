import sympy as sy
from sympy.abc import s
import numpy as np
from sympy.integrals import inverse_laplace_transform
import matplotlib.pyplot as plt

# Inputs for this step function
A1 = 0.27
A2 = 0.5
A3 = 0.23
t1 = 0
t2 = 0.55
t3 = 1.1

tvar = sy.symbols('tvar', real=True)

# Transfer function
G1 = 0.27*sy.exp(-0*s)
G2 = 0.5*sy.exp(-0.55*s)
G3 = 0.23*sy.exp(-1.11*s)

G = G1 + G2 + G3

# Taking inverse Laplace for time domain results
# G_t = inverse_laplace_transform(G, s, tvar)
G_t1 = inverse_laplace_transform(G1, s, tvar)
G_t2 = inverse_laplace_transform(G2, s, tvar)
G_t3 = inverse_laplace_transform(G3, s, tvar)
G_t = G_t2 + G_t3

h = 0.01
t0 = 0
tN = 240
N = int((tN-t0)/h) + 1
t = np.linspace(t0, tN, N)
T = t[N-1]

# x_0 = 0
# y_0 = 0.5
A_const = 0.5
B = 0.5
a = 2
b = 1
mpi = np.pi
# phi = np.zeros(N)
# phidot = np.zeros(N)
# phiddot = np.zeros(N)
# phitdot = np.zeros(N)

tau = 2*mpi*(-15*(t/T)**4+6*(t/T)**5+10*(t/T)**3)
taudot = 2*mpi*(-15*4*(1/T)*(t/T)**3+6*5*(1/T)*(t/T)**4+10*3*(1/T)*(t/T)**2)
tauddot = 2*mpi*(-15*4*3*(1/T)**2*(t/T)**2 + 6*5*4*(1/T)**2*(t/T)**3+10*3*2*(1/T)**2*(t/T))
tautdot = 2*mpi*(-15*4*3*2*(1/T)**3*(t/T) + 6*5*4*3*(1/T)**3*(t/T)**2 + 10*3*2*(1/T)**3)
# taujdot = 2*mpi*(-15*4*3*2*(1/T)**3 + 6*5*4*3*2*(1/T)**4*(t/T))

'''# Trajectory: Lemniscate'''
x_ref = A_const*np.sin(a*tau)
y_ref = B*np.cos(b*tau)
v_x = A_const*a*np.cos(a*tau)*taudot
v_y = -B*b*np.sin(b*tau)*taudot
a_x = -A_const*a*a*np.sin(a*tau)*taudot+A_const*a*np.cos(a*tau)*tauddot
a_y = -B*b*b*np.sin(b*tau)*taudot-B*b*np.sin(b*tau)*tauddot
psi_des = np.zeros(N)
psi_desdot = np.zeros(N)
psi_desddot = np.zeros(N)
A1 = 0.27
A2 = 0.5
A3 = 0.23
t1 = 0
t2 = 0.55
t3 = 1.11

for i in range(len(t)):
    if t[i] < 0.55:
        x_ref[i] = (0.27) * x_ref[i]
        y_ref[i] = (0.27) * y_ref[i]
        psi_des[i] = (0.27) * psi_des[i]
        v_x[i] = (0.27) * v_x[i]
        v_y[i] = (0.27) * v_y[i]
        psi_desdot[i] = (0.27) * psi_desdot[i]
        a_x[i] = (0.27) * a_x[i]
        a_y[i] = (0.27) * a_y[i]
        psi_desddot[i] = (0.27) * psi_desddot[i]
    elif 0.55 <= t[i] < 1.11:
        x_ref[i] = (0.27 + 0.5) * x_ref[i]
        y_ref[i] = (0.27 + 0.5) * y_ref[i]
        psi_des[i] = (0.27 + 0.5) * psi_des[i]
        v_x[i] = (0.27 + 0.5) * v_x[i]
        v_y[i] = (0.27 + 0.5) * v_y[i]
        psi_desdot[i] = (0.27 + 0.5) * psi_desdot[i]
        a_x[i] = (0.27 + 0.5) * a_x[i]
        a_y[i] = (0.27 + 0.5) * a_y[i]
        psi_desddot[i] = (0.27 + 0.5) * psi_desddot[i]
    else:
        x_ref[i] = (0.27 + 0.5 + 0.23) * x_ref[i]
        y_ref[i] = (0.27 + 0.5 + 0.23) * y_ref[i]
        psi_des[i] = (0.27 + 0.5 + 0.23) * psi_des[i]
        v_x[i] = (0.27 + 0.5 + 0.23) * v_x[i]
        v_y[i] = (0.27 + 0.5 + 0.23) * v_y[i]
        psi_desdot[i] = (0.27 + 0.5 + 0.23) * psi_desdot[i]
        a_x[i] = (0.27 + 0.5 + 0.23) * a_x[i]
        a_y[i] = (0.27 + 0.5 + 0.23) * a_y[i]
        psi_desddot[i] = (0.27 + 0.5 + 0.23) * psi_desddot[i]
