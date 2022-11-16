import numpy as np
from scipy.optimize import fsolve


class Controller_fdbklin:
    """
    Doing feedback linearization.
    Inputs:- are the positions, orientations and their derivatives, whereas
    output:- is the control input values (omega1, omega2, omega3, omega4) or (u1, u2, u3, u4)
    """
    def __init__(self, Kp, Kd):
        self.kpx, self.kpy, self.kpz, self.kppsi = Kp
        self.kdx, self.kdy, self.kdz, self.kdpsi = Kd

    def equations(self, parms):
        """
        Solves a system of nonlinear equations
        :param parms: vx, vy, vz
        :return: phi_d, theta_d, u1
        """
        v_x, v_y, v_z, m_q,m_l,Ixx,Iyy,Izz,g,l,cable_l,K,b,Ax,Ay,Az,u1,u2,u3,u4,vx, vy, vz,theta_l,phi_l,theta_ldot, phi_ldot, phi, theta, psi = parms
        sin = np.sin
        cos = np.cos
        return ((-Ax*cable_l*m_l*vx*cos(phi_l)**2*cos(theta_l)**2 + Ax*cable_l*m_l*vx*cos(theta_l)**2 - Ax*cable_l*m_l*vx - Ax*cable_l*m_q*vx - Ay*cable_l*m_l*vy*sin(phi_l)*sin(theta_l)*cos(theta_l) + Az*cable_l*m_l*vz*sin(phi_l)*cos(phi_l)*cos(theta_l)**2 - cable_l**2*m_l*m_q*phi_ldot**2*sin(phi_l)*cos(theta_l)**3 - cable_l**2*m_l*m_q*theta_ldot**2*sin(phi_l)*cos(theta_l) + cable_l*g*m_l**2*sin(phi_l)*cos(phi_l)*cos(theta_l)**2 + cable_l*g*m_l*m_q*sin(phi_l)*cos(phi_l)*cos(theta_l)**2 - cable_l*m_l*phi*u1*sin(phi_l)*sin(theta_l)*cos(psi)*cos(theta_l) + cable_l*m_l*phi*u1*sin(psi)*cos(phi_l)**2*cos(theta_l)**2 - cable_l*m_l*phi*u1*sin(psi)*cos(theta_l)**2 + cable_l*m_l*phi*u1*sin(psi) + cable_l*m_l*theta*u1*sin(phi_l)*sin(psi)*sin(theta_l)*cos(theta_l) + cable_l*m_l*theta*u1*cos(phi_l)**2*cos(psi)*cos(theta_l)**2 - cable_l*m_l*theta*u1*cos(psi)*cos(theta_l)**2 + cable_l*m_l*theta*u1*cos(psi) - cable_l*m_l*u1*sin(phi_l)*cos(phi_l)*cos(theta_l)**2 + cable_l*m_q*phi*u1*sin(psi) + cable_l*m_q*theta*u1*cos(psi) - g*l*m_l**2*sin(phi_l)*cos(phi_l)*cos(theta_l)**2 - g*l*m_l*m_q*sin(phi_l)*cos(phi_l)*cos(theta_l)**2)/(cable_l*m_q*(m_l + m_q)) - v_x,
                (-Ax*cable_l*m_l*vx*cos(phi_l - 2*theta_l)/4 + Ax*cable_l*m_l*vx*cos(phi_l + 2*theta_l)/4 - Ay*cable_l*m_l*vy*cos(2*theta_l)/2 - Ay*cable_l*m_l*vy/2 - Ay*cable_l*m_q*vy + Az*cable_l*m_l*vz*sin(phi_l - 2*theta_l)/4 - Az*cable_l*m_l*vz*sin(phi_l + 2*theta_l)/4 + cable_l**2*m_l*m_q*phi_ldot**2*sin(theta_l)/4 + cable_l**2*m_l*m_q*phi_ldot**2*sin(3*theta_l)/4 + cable_l**2*m_l*m_q*theta_ldot**2*sin(theta_l) + cable_l*g*m_l**2*sin(phi_l - 2*theta_l)/4 - cable_l*g*m_l**2*sin(phi_l + 2*theta_l)/4 + cable_l*g*m_l*m_q*sin(phi_l - 2*theta_l)/4 - cable_l*g*m_l*m_q*sin(phi_l + 2*theta_l)/4 + cable_l*m_l*phi*u1*sin(-phi_l + psi + 2*theta_l)/8 + cable_l*m_l*phi*u1*sin(phi_l - psi + 2*theta_l)/8 + cable_l*m_l*phi*u1*sin(phi_l + psi - 2*theta_l)/8 - cable_l*m_l*phi*u1*sin(phi_l + psi + 2*theta_l)/8 - cable_l*m_l*phi*u1*cos(psi)/2 - cable_l*m_l*phi*u1*cos(psi - 2*theta_l)/4 - cable_l*m_l*phi*u1*cos(psi + 2*theta_l)/4 + cable_l*m_l*theta*u1*sin(psi)/2 + cable_l*m_l*theta*u1*sin(psi - 2*theta_l)/4 + cable_l*m_l*theta*u1*sin(psi + 2*theta_l)/4 + cable_l*m_l*theta*u1*cos(-phi_l + psi + 2*theta_l)/8 - cable_l*m_l*theta*u1*cos(phi_l - psi + 2*theta_l)/8 + cable_l*m_l*theta*u1*cos(phi_l + psi - 2*theta_l)/8 - cable_l*m_l*theta*u1*cos(phi_l + psi + 2*theta_l)/8 - cable_l*m_l*u1*sin(phi_l - 2*theta_l)/4 + cable_l*m_l*u1*sin(phi_l + 2*theta_l)/4 - cable_l*m_q*phi*u1*cos(psi) + cable_l*m_q*theta*u1*sin(psi) - g*l*m_l**2*sin(phi_l - 2*theta_l)/4 + g*l*m_l**2*sin(phi_l + 2*theta_l)/4 - g*l*m_l*m_q*sin(phi_l - 2*theta_l)/4 + g*l*m_l*m_q*sin(phi_l + 2*theta_l)/4)/(cable_l*m_q*(m_l + m_q)) - v_y,
                (1.0*Ax*cable_l*m_l*vx*sin(phi_l)*cos(phi_l)*cos(theta_l)**2 - 1.0*Ay*cable_l*m_l*vy*sin(theta_l)*cos(phi_l)*cos(theta_l) + 1.0*Az*cable_l*m_l*vz*cos(phi_l)**2*cos(theta_l)**2 - 1.0*Az*cable_l*m_l*vz - 1.0*Az*cable_l*m_q*vz - 1.0*cable_l**2*m_l*m_q*phi_ldot**2*cos(phi_l)*cos(theta_l)**3 - 1.0*cable_l**2*m_l*m_q*theta_ldot**2*cos(phi_l)*cos(theta_l) + 1.0*cable_l*g*m_l**2*cos(phi_l)**2*cos(theta_l)**2 - 1.0*cable_l*g*m_l**2 + 1.0*cable_l*g*m_l*m_q*cos(phi_l)**2*cos(theta_l)**2 - 2.0*cable_l*g*m_l*m_q - 1.0*cable_l*g*m_q**2 - 1.0*cable_l*m_l*phi*u1*sin(phi_l)*sin(psi)*cos(phi_l)*cos(theta_l)**2 - 1.0*cable_l*m_l*phi*u1*sin(theta_l)*cos(phi_l)*cos(psi)*cos(theta_l) - 1.0*cable_l*m_l*theta*u1*sin(phi_l)*cos(phi_l)*cos(psi)*cos(theta_l)**2 + 1.0*cable_l*m_l*theta*u1*sin(psi)*sin(theta_l)*cos(phi_l)*cos(theta_l) - 1.0*cable_l*m_l*u1*cos(phi_l)**2*cos(theta_l)**2 + 1.0*cable_l*m_l*u1 + 1.0*cable_l*m_q*u1 - 1.0*g*l*m_l**2*cos(phi_l)**2*cos(theta_l)**2 + 1.0*g*l*m_l**2 - 1.0*g*l*m_l*m_q*cos(phi_l)**2*cos(theta_l)**2 + 1.0*g*l*m_l*m_q)/(cable_l*m_q*(m_l + m_q)) - v_z)

    def get_desired_positions(self, invA, B_control, D_sym, X, xd, yd, zd, psid, vxd, vyd, vzd, psidesdot, axd, ayd, azd, theta_l, theta_ldot, phi_l, phi_ldot, m_q,m_l,Ixx,Iyy,Izz,g,l,cable_l,K,b,Ax,Ay,Az,u1,u2,u3,u4):
        """
        To get the values of the control inputs to feed it into the dynamic equations such that XDDOT = V (feedback linearization)
        :param invA: inverse matrix of the Mass inertia matrix A
        :param B_control: B matrix of the EOM: AXDDOT + D = B*u
        :param D_sym: D matrix of the EOM: AXDDOT + D = B*u
        :return: u values (the nonlinear control inputs)
        """
        sin = np.sin
        cos = np.cos
        x, y, z, theta_l, phi_l, phi, theta, psi, vx, vy, vz, theta_ldot, phi_ldot, phidot, thetadot, psidot = X
        v_x = axd + self.kdx*(vxd - vx) + self.kpx*(xd - x)
        v_y = ayd + self.kdy*(vyd - vy) + self.kpy*(yd - y)
        v_z = azd + self.kdz*(vzd - vz) + self.kpz*(zd - z)
        v_thetal = 0
        v_phil = 0
        v_phi = 0
        v_theta = 0
        v_psi = self.kdpsi*(psidesdot - psidot) + self.kppsi*(psid - psi)

        V = np.array([v_x, v_y, v_z, v_thetal, v_phil, v_phi, v_theta, v_psi])
        V = V.reshape(len(V),1)

        val = invA.dot(B_control)
        inverse_val = np.linalg.pinv(val)
        AinvD = invA.dot(D_sym)

        U = inverse_val.dot(V) + inverse_val.dot(AinvD)
        print("*************************************************************************\n")
        print("Thrust value is: {}".format(U[0]))
        print("Roll torque is: {}".format(U[1]))
        print("Pitch torque is: {}".format(U[2]))
        print("Yaw torque is: {}".format(U[3]))
        print("*************************************************************************\n")
        parms = v_x, v_y, v_z, m_q,m_l,Ixx,Iyy,Izz,g,l,cable_l,K,b,Ax,Ay,Az,u1,u2,u3,u4,vx, vy, vz,theta_l,phi_l,theta_ldot, phi_ldot, phi, theta, psi
        phi_d, theta_d, u1 = fsolve(self.equations(parms), (v_x, v_y, v_z))

        return U, phi_d, theta_d



    def control_action(self, ang, desired_state):
        """
        Does the torque control of quadcopter (pitch and roll torques)
        :return:
        """

