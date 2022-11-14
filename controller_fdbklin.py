import numpy as np
# import numpy.typing as npt


class Controller:
    def __init__(self, K_z, K_psi, K_theta, K_phi, Kp, Kd, Kdd, Ki, Kx, Kz, A, k: float,
                 l: float, b: float) -> None:
        self.K_z = K_z
        self.K_psi = K_psi
        self.K_theta = K_theta
        self.K_phi = K_phi
        self.A = A[0]  # Considering Ax = Ay = Az = 0.25
        self.k = k
        self.l = l
        self.cable_l = 0.2
        self.b = b
        self.last_ang_pos = np.zeros(3)
        self.integral_error_x = 0
        self.integral_error_y = 0
        self.integral_error_z = 0

        # Constants
        self.g = 9.81
        self.m_q = 0.468
        self.m_l = 0.1
        self.I = np.array([4.856*1e-3, 4.856*1e-3, 8.801*1e-3])
        self.Kx, self.Ky, self.Kz = Kp
        self.Kdx, self.Kdy, self.Kdz = Kd
        self.Kddx, self.Kddy, self.Kddz = Kdd
        self.Kix, self.Kiy, self.Kiz = Ki
        self.kpx, self.kdx = Kx
        self.kpz, self.kdz = Kz

    def get_desired_positions(self, x, vx, z, vz, x_des, z_des, xdot_des, zdot_des, xddot_des, zddot_des, phi_l, phi_ldot):
        u_z = -self.kpz*(z - z_des) - self.kdz*(vz - zdot_des) + zddot_des - 1.5
        u_x = -self.kpx*(x - x_des) - self.kdx*(vx - xdot_des) + xddot_des

        theta_d = (u_x + ((self.m_l*self.cable_l*np.sin(phi_l)*(phi_ldot**2))/(self.m_l+self.m_q)) + ((self.m_l*self.g*np.sin(2*phi_l))/(2*self.m_q))) / (self.g + ((self.g*self.m_l)/(2*self.m_q))*(1 - np.cos(2*phi_l)))
        u1 = ((self.m_q + self.m_l)*(self.g + (self.m_l*self.cable_l*np.cos(phi_l)*phi_ldot**2)/(self.m_q + self.m_l) + u_z)) / (1 + (self.m_l/2*self.m_q)*(1 - np.cos(2*phi_l) - theta_d*np.sin(2*phi_l)))

        return theta_d, u1

    def choose_action(self, ang, desired_state):
        e_theta = ang[0] - desired_state[0]
        omega_yd = -self.K_theta[0]*e_theta
        e_omega = ang[1] - omega_yd
        T_theta = (-self.K_theta[1]*e_omega) * self.I[1]
        return T_theta

