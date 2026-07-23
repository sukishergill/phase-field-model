# This file contains helper functions for the SSAV scheme

import numpy as np
from numpy.fft import fftn, ifftn, fftfreq

def compute_E(ufft, w, Grid, Para, D, dim):
    # Compute free energy
    #En = w**22 - Para.B

    #e = 0.5*Para.eps**2
    ux = np.real(ifftn(-1j * Grid.kxx * ufft))
    uy = np.real(ifftn(-1j * Grid.kyy * ufft))

    du = ux**2 + uy**2

    return np.sum(0.5 * Para.eps**2 * du)*Grid.dV + w**2 - Para.B
    
    #En += np.sum(e*du) * Grid.dV
    # if Para.model != "PFC":

    #     ux = np.real(ifftn(-1j * Grid.kxx * ufft))
    #     du = ux**2

    #     if dim != 1:
    #         uy = np.real(ifftn(-1j * Grid.kyy * ufft))
    #         du += uy**2

    #         if dim == 3:
    #             uz = np.real(ifftn(-1j * Grid.kzz * ufft))
    #             du += uz**2

    #     En += np.sum(e*du) * Grid.dV

    #     if Para.model == "OK":
    #         um_fft = ufft - Para.m
    #         v = -Grid.inv_k * um_fft
    #         v[Grid.k == 0] = 0;

    #         vx = np.real(ifftn(-1j*Grid.kxx*v))
    #         dv = vx**2

    #         if dim != 1:
    #             vy = np.real(ifftn(-1j*Grid.kyy*v))
    #             dv += vy**2
    
    #             if dim == 3:
    #                 vz = np.real(ifftn(-1j*Grid.kzz*v))
    #                 dv += vz**2
    
    #         En += np.sum(e*Para.alpha * dv) * Grid.dV

    #return En

def compute_us(u_curr, u_prev, gamma):
    return (gamma+1) * u_curr - gamma * u_prev

def compute_F(u, beta):
    F = 0.25 * (u**2 - beta)**2
    f = u**3 - beta*u

    return F, f


def compute_w(int_F, B):
    return np.sqrt(int_F + B)

def compute_H(f, w):
    return  f / w

def compute_dt_coeffs(dt, dt_new):

    a = 1/ dt_new + 1 / (dt_new + dt)
    b = -1 / dt_new - 1 / dt
    c = 1 / dt - 1 / (dt_new + dt)

    return a, b, c


def compute_r(us, us_fft, u_curr, u_prev, w_curr, w_prev, Hs, H2, a, b, c, Grid, Para, G, t, dt):

    r = -(b*w_curr + c*w_prev)/a + 0.5 * np.sum(Hs * (b*u_curr + c*u_prev)/a) * Grid.dV;

    r_tilde = -(b*u_curr + c*u_prev) - Para.S*np.real(ifftn(G*us_fft)) + r*np.real(ifftn(G*H2))

    m_est = (Para.alpha*Para.eps**2) / (a*Grid.V) * np.sum(r_tilde)*Grid.dV

    r_hat = r_tilde + m_est

    return r, r_hat

def compute_ip(Hs, psi_r, psi_H, Grid):

    innprod_Hu = np.sum(Hs * psi_r)*Grid.dV / (1 - 0.5*np.sum(Hs*psi_H)*Grid.dV)

    return innprod_Hu


def compute_unew(u_curr, ucurr_fft, u_prev, uprev_fft, w_curr, w_prev, dt, dt_new, Para, Grid, G, P1, t):

    us = compute_us(u_curr, u_prev, dt_new/dt)

    F, f = compute_F(us, Para.beta)

    int_F = np.sum(F)*Grid.dV

    ws = compute_w(int_F, Para.B)

    Hs = compute_H(f, ws)

    a,b,c = compute_dt_coeffs(dt, dt_new)

    H2 = fftn(Hs)

    r, r_hat = compute_r(us, 2*ucurr_fft - uprev_fft, u_curr, u_prev, w_curr, w_prev, Hs, H2, a, b, c, Grid, Para, G, t + dt, dt)

    P = a + P1

    r_hat = fftn(r_hat)
    psi_r = r_hat / P
    psi_r = np.real(ifftn(psi_r))

    psi_H = (G * H2) / P
    psi_H = np.real(ifftn(psi_H))

    innprod_Hu = compute_ip(Hs, psi_r, psi_H, Grid)

    w_new = 0.5*innprod_Hu + r

    u_new = 0.5*innprod_Hu*psi_H + psi_r

    return u_new, w_new


def compute_dtnew(E_t, err_tol, delta, dt, dt_min, dt_max):

    dt_1 = err_tol[0] / (E_t + delta)
    dt_2 = err_tol[1] / np.sqrt(E_t + delta**2)

    dt_prop = min(dt_1, dt_2)

    min_dt = min(dt_max, 2*dt, dt_prop)

    dt_new = max(dt_min, min_dt, 0.1*dt)

    return dt_new
