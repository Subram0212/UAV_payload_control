import numpy as np
from scipy.optimize import fsolve


class Controller_fdbklin:
    """
    Doing feedback linearization.
    Inputs:- are the positions, orientations and their derivatives, whereas
    output:- is the control input values (omega1, omega2, omega3, omega4) or (u1, u2, u3, u4)
    """
    def __init__(self, Kp, Kd, Kdd, Ki, K_phi, K_theta, K_psi, Ixx, Iyy, Izz):
        self.kpx, self.kpy, self.kpz, self.kppsi = Kp
        self.kdx, self.kdy, self.kdz, self.kdpsi = Kd
        self.kddx, self.kddy, self.kddz, self.kddpsi = Kdd
        self.K_phi, self.K_phid = K_phi
        self.K_theta, self.K_thetad = K_theta
        self.K_psi, self.K_psid = K_psi
        self.Ixx = Ixx
        self.Iyy = Iyy
        self.Izz = Izz
        self.integral_error_x = 0
        self.integral_error_y = 0
        self.integral_error_z = 0
        self.integral_error_psi = 0
        self.Kxi, self.Kyi, self.Kzi, self.Kpsii = Ki

    def equations(self, variables, v_x, v_y, v_z, m_q,m_l,Ixx,Iyy,Izz,g,l,cable_l,K,b,Ax,Ay,Az,u2,u3,u4,vx, vy, vz,theta_l,phi_l,theta_ldot, phi_ldot, psi):
        """
        Solves a system of nonlinear equations to get desired roll and pitch angles. The desired roll and pitch
        angles mean that the drone flight should be maintained such a way that those angles (roll, pitch and yaw) should
        be as small as possible, and have a constant thrust throughout its flight to maintain stability.
        :param vx, vy
        :return: phi_d, theta_d
        """
        phi, theta = variables
        sin = np.sin
        cos = np.cos
        return ((-Ax*cable_l*m_l*vx*cos(phi_l)**2*cos(theta_l)**2 + Ax*cable_l*m_l*vx*cos(theta_l)**2 - Ax*cable_l*m_l*vx - Ax*cable_l*m_q*vx - Ay*cable_l*m_l*vy*sin(phi_l)*sin(theta_l)*cos(theta_l) + Az*cable_l*m_l*vz*sin(phi_l)*cos(phi_l)*cos(theta_l)**2 - cable_l**2*m_l*m_q*phi_ldot**2*sin(phi_l)*cos(theta_l)**3 - cable_l**2*m_l*m_q*theta_ldot**2*sin(phi_l)*cos(theta_l) + cable_l*g*m_l**2*sin(phi_l)*cos(phi_l)*cos(theta_l)**2 + cable_l*g*m_l*m_q*sin(phi_l)*cos(phi_l)*cos(theta_l)**2 - cable_l*m_l*phi*(m_l+m_q)*g*sin(phi_l)*sin(theta_l)*cos(psi)*cos(theta_l) + cable_l*m_l*phi*(m_l+m_q)*g*sin(psi)*cos(phi_l)**2*cos(theta_l)**2 - cable_l*m_l*phi*(m_l+m_q)*g*sin(psi)*cos(theta_l)**2 + cable_l*m_l*phi*(m_l+m_q)*g*sin(psi) + cable_l*m_l*theta*(m_l+m_q)*g*sin(phi_l)*sin(psi)*sin(theta_l)*cos(theta_l) + cable_l*m_l*theta*(m_l+m_q)*g*cos(phi_l)**2*cos(psi)*cos(theta_l)**2 - cable_l*m_l*theta*(m_l+m_q)*g*cos(psi)*cos(theta_l)**2 + cable_l*m_l*theta*(m_l+m_q)*g*cos(psi) - cable_l*m_l*(m_l+m_q)*g*sin(phi_l)*cos(phi_l)*cos(theta_l)**2 + cable_l*m_q*phi*(m_l+m_q)*g*sin(psi) + cable_l*m_q*theta*(m_l+m_q)*g*cos(psi) - g*l*m_l**2*sin(phi_l)*cos(phi_l)*cos(theta_l)**2 - g*l*m_l*m_q*sin(phi_l)*cos(phi_l)*cos(theta_l)**2)/(cable_l*m_q*(m_l + m_q)) - v_x,
                (-Ax*cable_l*m_l*vx*cos(phi_l - 2*theta_l)/4 + Ax*cable_l*m_l*vx*cos(phi_l + 2*theta_l)/4 - Ay*cable_l*m_l*vy*cos(2*theta_l)/2 - Ay*cable_l*m_l*vy/2 - Ay*cable_l*m_q*vy + Az*cable_l*m_l*vz*sin(phi_l - 2*theta_l)/4 - Az*cable_l*m_l*vz*sin(phi_l + 2*theta_l)/4 + cable_l**2*m_l*m_q*phi_ldot**2*sin(theta_l)/4 + cable_l**2*m_l*m_q*phi_ldot**2*sin(3*theta_l)/4 + cable_l**2*m_l*m_q*theta_ldot**2*sin(theta_l) + cable_l*g*m_l**2*sin(phi_l - 2*theta_l)/4 - cable_l*g*m_l**2*sin(phi_l + 2*theta_l)/4 + cable_l*g*m_l*m_q*sin(phi_l - 2*theta_l)/4 - cable_l*g*m_l*m_q*sin(phi_l + 2*theta_l)/4 + cable_l*m_l*phi*(m_l+m_q)*g*sin(-phi_l + psi + 2*theta_l)/8 + cable_l*m_l*phi*(m_l+m_q)*g*sin(phi_l - psi + 2*theta_l)/8 + cable_l*m_l*phi*(m_l+m_q)*g*sin(phi_l + psi - 2*theta_l)/8 - cable_l*m_l*phi*(m_l+m_q)*g*sin(phi_l + psi + 2*theta_l)/8 - cable_l*m_l*phi*(m_l+m_q)*g*cos(psi)/2 - cable_l*m_l*phi*(m_l+m_q)*g*cos(psi - 2*theta_l)/4 - cable_l*m_l*phi*(m_l+m_q)*g*cos(psi + 2*theta_l)/4 + cable_l*m_l*theta*(m_l+m_q)*g*sin(psi)/2 + cable_l*m_l*theta*(m_l+m_q)*g*sin(psi - 2*theta_l)/4 + cable_l*m_l*theta*(m_l+m_q)*g*sin(psi + 2*theta_l)/4 + cable_l*m_l*theta*(m_l+m_q)*g*cos(-phi_l + psi + 2*theta_l)/8 - cable_l*m_l*theta*(m_l+m_q)*g*cos(phi_l - psi + 2*theta_l)/8 + cable_l*m_l*theta*(m_l+m_q)*g*cos(phi_l + psi - 2*theta_l)/8 - cable_l*m_l*theta*(m_l+m_q)*g*cos(phi_l + psi + 2*theta_l)/8 - cable_l*m_l*(m_l+m_q)*g*sin(phi_l - 2*theta_l)/4 + cable_l*m_l*(m_l+m_q)*g*sin(phi_l + 2*theta_l)/4 - cable_l*m_q*phi*(m_l+m_q)*g*cos(psi) + cable_l*m_q*theta*(m_l+m_q)*g*sin(psi) - g*l*m_l**2*sin(phi_l - 2*theta_l)/4 + g*l*m_l**2*sin(phi_l + 2*theta_l)/4 - g*l*m_l*m_q*sin(phi_l - 2*theta_l)/4 + g*l*m_l*m_q*sin(phi_l + 2*theta_l)/4)/(cable_l*m_q*(m_l + m_q)) - v_y)

    def get_desired_positions(self, t_step, invA, B_control, D_sym, X, Xddot, xd, yd, zd, psid, vxd, vyd, vzd, psidesdot, axd, ayd, azd, psidesddot, m_q,m_l,Ixx,Iyy,Izz,g,l,cable_l,K,b,Ax,Ay,Az,u1,u2,u3,u4):
        """
        To get the values of the control inputs to feed it into the dynamic equations such that XDDOT = V (feedback linearization)
        :param invA: inverse matrix of the Mass inertia matrix A
        :param B_control: B matrix of the EOM: AXDDOT + D = B*u
        :param D_sym: D matrix of the EOM: AXDDOT + D = B*u
        :return: u values (the nonlinear control inputs) -> specifically its 'u1' and 'u4' are being applied to the dynamics.
        """
        t = t_step[1] - t_step[0]
        sin = np.sin
        cos = np.cos
        x, y, z, theta_l, phi_l, phi, theta, psi, vx, vy, vz, theta_ldot, phi_ldot, phidot, thetadot, psidot = X
        vx, vy, vz, theta_ldot, phi_ldot, phidot, thetadot, psidot, ax, ay, az, theta_lddot, phi_lddot, phiddot, thetaddot, psiddot = Xddot
        ex = (xd - x)
        ey = (yd - y)
        ez = (zd - z)
        epsi = (psid - psi)
        exd = (vxd - vx)
        eyd = (vyd - vy)
        ezd = (vzd - vz)
        exdd = (axd - ax)
        eydd = (ayd - ay)
        ezdd = (azd - az)
        self.integral_error_x += ex * t
        self.integral_error_y += ey * t
        self.integral_error_z += ez * t
        self.integral_error_psi += epsi * t
        v_x = self.kddx*(axd - ax) + self.Kxi*self.integral_error_x + self.kdx*(vxd - vx) + self.kpx*(xd - x)
        v_y = self.kddx*(ayd - ay) + self.Kyi*self.integral_error_y + self.kdy*(vyd - vy) + self.kpy*(yd - y)
        v_z = self.kddx*(azd - az) + self.Kzi*self.integral_error_z + self.kdz*(vzd - vz) + self.kpz*(zd - z)
        v_thetal = 0
        v_phil = 0
        v_phi = 0
        v_theta = 0
        v_psi = self.kddpsi*(psidesddot - psiddot) + self.Kpsii*self.integral_error_psi + self.kdpsi*(psidesdot - psidot) + self.kppsi*(psid - psi)

        V = np.array([v_x, v_y, v_z, v_thetal, v_phil, v_phi, v_theta, v_psi])
        V = V.reshape(len(V),1)

        val = invA.dot(B_control)
        inverse_val = np.linalg.pinv(val)
        AinvD = invA.dot(D_sym)

        U = inverse_val.dot(V) + inverse_val.dot(AinvD)
        U = U.reshape(4,)
        print("*************************************************************************\n")
        print("Thrust value is: {}".format(U[0]))
        print("Roll torque is: {}".format(U[1]))
        print("Pitch torque is: {}".format(U[2]))
        print("Yaw torque is: {}".format(U[3]))
        print("*************************************************************************\n")
        parms = v_x, v_y, v_z, m_q,m_l,Ixx,Iyy,Izz,g,l,cable_l,K,b,Ax,Ay,Az,u2,u3,u4,vx, vy, vz,theta_l,phi_l,theta_ldot, phi_ldot, psi
        variables = 0, 0
        phi_d, theta_d = fsolve(self.equations, (variables), args=parms)

        return U, phi_d, theta_d

    def control_action(self, ang, desired_state):
        """
        Does the torque control of quadcopter (pitch and roll torques)
        :return: Torques for roll, pitch, yaw
        """
        phi, phidot, theta, thetadot, psi, psidot = ang
        phi_des, theta_des, psi_des = desired_state

        # e_phi = phi - phi_des
        # omega_xd = -self.K_phi*e_phi
        # e_omega_phi = phidot - omega_xd
        # T_phi = (-self.K_thetad*e_omega_phi) * self.Ixx
        #
        # e_theta = theta - theta_des
        # omega_yd = -self.K_theta*e_theta
        # e_omega_theta = thetadot - omega_yd
        # T_theta = (-self.K_thetad*e_omega_theta) * self.Iyy

        T_phi = (self.K_phid*(-phidot) + self.K_phi*(phi_des-phi))*self.Ixx
        T_theta = (self.K_thetad*(-thetadot) + self.K_theta*(theta_des-theta))*self.Iyy
        # e_psi = psi - psi_des
        # omega_zd = -self.K_psi*e_psi
        # e_omega_psi = psidot - omega_zd
        # T_psi = (-self.K_thetad*e_omega_psi) * self.Izz

        T_psi = (self.K_psid*(-psidot) + self.K_psi*(psi_des-psi))*self.Izz

        return T_phi, T_theta



