# SSAV

# import libraries
import ssav_helpers
from ssav_helpers import *
import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
from dataclasses import dataclass

# Create a data class for the output
@dataclass
class Output:
    u: np.ndarray
    Eu: np.ndarray
    tt: np.ndarray


def ssav(Grid, Time, Para, u, dim):

    dt = Time.dt_min
    tt = []
    tt.append(0)

    Nt = int(Time.tf / dt)

    if Para.model == "PFC":
        D = (-Grid.k2 + 1)

    else:
        D = -1j*Grid.k

    if Para.model == "AC":
        G = -1

    else:
        G = -Grid.k2

    ufft = fftn(u)

    ##################################################
    ######## Compute u^1 using backward Euler ########
    ##################################################
    F, f = compute_F(u, Para.beta)
    int_F = np.sum(F) * Grid.dV

    # Define this as w_old b/c we will need it in the loop below
    w_old = compute_w(int_F, Para.B)

    E = []
    E.append(compute_E(ufft, w_old, Grid, Para))

    H = compute_H(f, w_old)

    r = -0.5 * np.sum(H*u)*Grid.dV + w_old
    H2 = fftn(H)

    r_tilde = u/dt - Para.S*np.real(ifftn(G*ufft)) - r*np.real(ifftn(G*H2))

    r_hat = r_tilde

    # Define linear operator
    if Para.model == "PFC":
        P1 = -Para.eps**2*G*(D**2) - Para.S*G

    else:
        P1 = Para.eps**2*G*(D**2) - Para.S*G

    P = 1/dt + P1

    r_hat = fftn(r_hat)
    psi_r = r_hat / P
    psi_r = np.real(ifftn(psi_r))

    psi_H = (G * H2) / P
    psi_H = np.real(ifftn(psi_H))

    innprod_Hu = compute_ip(H, psi_r, psi_H, Grid)

    w = 0.5*innprod_Hu + r

    u_old = u
    uold_fft = ufft

    u = 0.5*innprod_Hu * psi_H + psi_r
    ufft = fftn(u)

    E_curr = compute_E(ufft, w, Grid, Para)
    E.append(E_curr)

    t = dt
    tt.append(t)
    dt_new = dt
    
    ##################################################
    ################# Main time loop #################
    ##################################################
    for i in range(Nt):

        u_new, w_new = compute_unew(u, ufft, u_old, uold_fft, w, w_old, dt, dt_new,
                                    Para, Grid,G, P1, t)

        unew_fft = fftn(u_new)

        E_new = compute_E(unew_fft, w_new, Grid, Para)

        E_t = (E_new - E_curr) / dt

        if Time.dt_min == Time.dt_max:
            dt_new = dt

        else:
            dt_new = compute_dtnew(abs(E_t), Para.err_tol, Para.delta, dt,
                                   Time.dt_min, Time.dt_max)

        E.append(E_new)
        E_curr = E_new
        
        w_old = w
        w = w_new

        u_old = u
        u = u_new

        uold_fft = ufft
        ufft = unew_fft

        t += dt
        tt.append(t)

    Results = Output(u=u, Eu=E, tt=tt)

    return Results
