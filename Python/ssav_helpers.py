# This file contains helper functions for the SSAV scheme

import numpy as np
from numpy.fft import fftn, ifftn, fftfreq

def compute_E(ufft, w, Grid, Para):

    ux = np.real(ifftn(-1j * Grid.kxx * ufft))
    uy = np.real(ifftn(-1j * Grid.kyy * ufft))

    du = ux**2 + uy**2

    return np.sum(0.5 * Para.eps**2 * du)*np.prod(Grid.d) + w**2 - Para.B

def compute_us(u_curr, u_prev):
    return 2 * u_curr - u_prev

def compute_F(u):
    F = 0.25 * (u**2 - 1)**2
    f = u**3 - u

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

    r = -(b*w_curr + c*w_prev)/a + 0.5 * np.sum(Hs * (b*u_curr + c*u_prev)/a) * np.prod(Grid.d);

    r_tilde = -(b*u_curr + c*u_prev) - Para.S*np.real(ifftn(G*us_fft)) + r*np.real(ifftn(G*H2))

    # m_est

    r_hat = r_tilde

    return r, r_hat

def compute_ip(Hs, psi_r, psi_H, Grid):

    innprod_Hu = np.sum(Hs * psi_r)*np.prod(Grid.d) / (1 - 0.5*np.sum(Hs*psi_H)*np.prod(Grid.d))

    return innprod_Hu


def compute_unew(u_curr, ucurr_fft, u_prev, uprev_fft, w_curr, w_prev, dt, Para, Grid, G, P1, t):

    us = compute_us(u_curr, u_prev)

    F, f = compute_F(us)

    int_F = np.sum(F)*np.prod(Grid.d)

    ws = compute_w(int_F, Para.B)

    Hs = compute_H(f, ws)

    a,b,c = compute_dt_coeffs(dt, dt)

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