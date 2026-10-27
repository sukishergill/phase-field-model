# This file contains helper functions for the SSAV scheme

import numpy as np
from numpy.fft import fftn, ifftn, fftfreq

def compute_E(ufft, w, Grid, Para, D, dim):
    # Compute free energy
    En = w**2 - Para.B

    e = 0.5*Para.eps**2

    if Para.model == "PFC":
        du = np.real(ifftn(D*ufft))
        En += np.sum(e*du) * Grid.dV
    
    elif Para.model in ("CH", "OK", "FCH"):

        ux = np.real(ifftn(-1j * Grid.kxx * ufft))
        du = ux**2

        if dim != 1:
            uy = np.real(ifftn(-1j * Grid.kyy * ufft))
            du += uy**2

        if dim == 3:
            uz = np.real(ifftn(-1j * Grid.kzz * ufft))
            du += uz**2

        if Para.model != "FCH":
            En += np.sum(e*du) * Grid.dV

        if Para.model == "OK":
            um_fft = ufft - Para.m
            v = -Grid.inv_k * um_fft
            v = np.where(Grid.k == 0, 0, v)

            vx = np.real(ifftn(-1j*Grid.kxx*v))
            dv = vx**2

            if dim != 1:
                vy = np.real(ifftn(-1j*Grid.kyy*v))
                dv += vy**2
    
            if dim == 3:
                vz = np.real(ifftn(-1j*Grid.kzz*v))
                dv += vz**2
    
            En += np.sum(e*Para.alpha * dv) * Grid.dV

        if Para.model == "FCH":

            En += np.sum(0.5 * Para.eps**4 * np.real(ifftn(-Grid.k**2 * ufft))**2
                          - (Para.tau**2 / 3 + Para.eta1 / 2 + 1) * 
                         Para.eps**2 * du) * Grid.dV

    return En/Grid.V

def compute_us(u_curr, u_prev, gamma):
    return (gamma+1) * u_curr - gamma * u_prev

def compute_F(u, beta, tau):
    if tau == 0:
        F = 0.25 * (u**2 - beta)**2
        F1 = u**3 - beta*u
        F2 = 3*u**2 - beta

    else:
        F = 0.5 * (u + tau/3 + 1)**2 * (0.5 * (u + tau/3 - 1)**2
        + (2 * tau/3) * (u + tau/3 - 2))

        F1 = (u**3- (1 + tau**2/3) * u - (2 * tau/3) * (1 + tau/3) * (1 - tau/3))
        
        F2 = 3 * u**2 - (1 + tau**2/3)

    return F, F1, F2

def compute_G(u, ufft, Para, F, F1, Grid, dim):

    if Para.model == "FCH":

        if dim == 1:
            ux = np.real(ifftn(-1j*Grid.k*ufft))
            
            G = np.sum(3*(Para.eps*u*ux)**2 + 0.5*F1**2 - Para.eta2*F)*Grid.dV

        else:
            ux = np.real(ifftn(-1j*Grid.kxx*ufft))
            uy = np.real(ifftn(-1j*Grid.kyy*ufft))

            G = np.sum(3*(Para.eps*u)**2 * (ux**2 + uy**2) + 0.5* F1**2 - 
                       Para.eta2*F) * Grid.dV

    else:
        G = np.sum(F) * Grid.dV

    return G


def compute_w(G, B):
    return np.sqrt(G + B)

def compute_H(f, w):

    return  f / w

def compute_H_FCH(u, w, Para, F1, F2, Grid, dim):

    v = fftn(u)
    v2 = fftn(u**2)

    if dim == 1:
        ux = np.real(ifftn(-1j * Grid.k * v))
        ux2 = np.real(ifftn(-1j * Grid.k * v2))

        u_grad = ux**2

        H = (6*Para.eps**2 * (u*u_grad - ux*ux2 - (u**2)*np.real(ifftn(-Grid.k**2*v))) 
             + F1*F2 - Para.eta2 * F1) / w
    else:
        ux = np.real(ifftn(-1j * Grid.kxx * v))
        uy = np.real(ifftn(-1j * Grid.kyy * v))

        ux2 = np.real(ifftn(-1j * Grid.kxx * v2))
        uy2 = np.real(ifftn(-1j * Grid.kyy * v2))

        u_grad = ux**2 + uy**2

        H = (6*Para.eps**2 * (u*u_grad - ux*ux2 - uy*uy2 - 
                              u**2 * np.real(ifftn(-Grid.k**2*v))) + 
             F1*F2 - Para.eta2*F1) / w

    return H

def compute_dt_coeffs(dt, dt_new):

    a = 1/ dt_new + 1 / (dt_new + dt)
    b = -1 / dt_new - 1 / dt
    c = 1 / dt - 1 / (dt_new + dt)

    return a, b, c


def compute_r(us, us_fft, u_curr, u_prev, w_curr, w_prev, Hs, H2, a, b, c, Grid, Para, L1, t, dt):

    r = -(b*w_curr + c*w_prev)/a + 0.5 * np.sum(Hs * (b*u_curr + c*u_prev)/a) * Grid.dV;

    r_tilde = -(b*u_curr + c*u_prev) + np.real(ifftn(L1*us_fft)) + r*np.real(ifftn(-Grid.k2*H2))

    m_est = (Para.alpha*Para.eps**2) / (a*Grid.V) * np.sum(r_tilde)*Grid.dV

    r_hat = r_tilde + m_est

    return r, r_hat

def compute_ip(Hs, psi_r, psi_H, Grid):

    innprod_Hu = np.sum(Hs * psi_r)*Grid.dV / (1 - 0.5*np.sum(Hs*psi_H)*Grid.dV)

    return innprod_Hu


def compute_unew(u_curr, ucurr_fft, u_prev, uprev_fft, w_curr, w_prev, dt, dt_new, Para, Grid, L0, L1, t, dim):

    gamma = dt_new/dt
    us = compute_us(u_curr, u_prev, gamma)

    F, F1, F2 = compute_F(us, Para.beta, Para.tau)

    # int_F = np.sum(F)*Grid.dV
    G = compute_G(u_curr, ucurr_fft, Para, F, F1, Grid, dim)

    ws = compute_w(G, Para.B)

    if Para.model == "FCH":
        Hs = compute_H_FCH(us, ws, Para, F1, F2, Grid, dim)

    else:
        Hs = compute_H(F1, ws)

    a,b,c = compute_dt_coeffs(dt, dt_new)

    H2 = fftn(Hs)

    r, r_hat = compute_r(us, (1 + gamma) * ucurr_fft - gamma*uprev_fft, u_curr, u_prev, w_curr, w_prev, Hs, H2, a, b, c, Grid, Para, L1, t + dt, dt)

    Ln = a + L0

    r_hat = fftn(r_hat)
    psi_r = r_hat / Ln
    psi_r = np.real(ifftn(psi_r))

    psi_H = (-Grid.k2 * H2) / Ln
    psi_H = np.real(ifftn(psi_H))

    innprod_Hu = compute_ip(Hs, psi_r, psi_H, Grid)

    w_new = 0.5*innprod_Hu + r

    u_new = 0.5*innprod_Hu*psi_H + psi_r

    return u_new, w_new


def compute_dtnew(E_t, err_tol, delta, dt, dt_min, dt_max):

    dt_1 = err_tol[0] / (E_t + delta)
    dt_2 = err_tol[1] / np.sqrt(E_t + delta**2)

    dt_prop = min(dt_1, dt_2)
    min_dt = min(dt_max, max(dt_min, dt_prop))
    dt_new = min(2*dt, max(0.1*dt, min_dt))

    return dt_new

