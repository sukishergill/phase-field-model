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
    dt: np.ndarray


def ssav(Grid, Time, Para, u, dim):

    dt = Time.dt_min
    tt = []
    tt.append(0)

    Nt = int(Time.tf / dt)

    if Para.model == "PFC":
        D = (-Grid.k2 + 1)

    else:
        D = -1j*Grid.k

    # Define linear operators
    if Para.model == "FCH":
        L0 = Para.eps**4*Grid.k6 + Para.S*Grid.k4
        L1 = (Para.eps**2 * (2 + Para.eta1 + 2*Para.tau**2/3) + Para.S) * Grid.k4

    else:
        L1 = Para.S*Grid.k2
        
        if Para.model == "PFC":
            L0 = -Para.eps**2*Grid.k2*(D**2) + Para.S*Grid.k2

        else:
            L0 = Para.alpha*Para.eps**2 - Para.eps**2*Grid.k2*(D**2) + Para.S*Grid.k2


    ufft = fftn(u)

    ##################################################
    ######## Compute u^1 using backward Euler ########
    ##################################################
    F, F1, F2 = compute_F(u, Para.beta, Para.tau)
    G = compute_G(u, ufft, Para, F, F1, Grid, dim)

    # Define this as w_old b/c we will need it in the loop below
    w_old = compute_w(G, Para.B)

    E = []
    E.append(compute_E(ufft, w_old, Grid, Para, D, dim))

    if Para.model == "FCH":
        H = compute_H_FCH(u, w_old, Para, F1, F2, Grid, dim)

    else:
        H = compute_H(F1, w_old)

    r = -0.5 * np.sum(H*u)*Grid.dV + w_old
    H2 = fftn(H)

    r_tilde = u/dt + np.real(ifftn(L1*ufft)) + r*np.real(ifftn(-Grid.k2*H2))
    
    r_hat = r_tilde + (Para.alpha*Para.eps**2) / (Grid.V) * np.sum(r_tilde)*Grid.dV
    
    Ln = 1/dt + L0

    r_hat = fftn(r_hat)
    psi_r = r_hat / Ln
    psi_r = np.real(ifftn(psi_r))

    psi_H = (-Grid.k2 * H2) / Ln
    psi_H = np.real(ifftn(psi_H))

    innprod_Hu = compute_ip(H, psi_r, psi_H, Grid)

    w = 0.5*innprod_Hu + r

    u_old = u
    uold_fft = ufft

    u = 0.5*innprod_Hu * psi_H + psi_r
    ufft = fftn(u)

    E_curr = compute_E(ufft, w, Grid, Para, D, dim)
    E.append(E_curr)

    t = dt
    tt.append(t)
    dt_new = dt
    
    ##################################################
    ################# Main time loop #################
    ##################################################
    while t < Time.tf:

        step_dt = dt_new

        u_new, w_new = compute_unew(u, ufft, u_old, uold_fft, w, w_old, dt, step_dt,
                                    Para, Grid, L0, L1, t, dim)

        unew_fft = fftn(u_new)

        E_new = compute_E(unew_fft, w_new, Grid, Para, D, dim)

        E_t = (E_new - E_curr) / step_dt

        if Time.dt_min == Time.dt_max:
            dt_new = dt

        else:
            dt = step_dt
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

        t += step_dt
        tt.append(t)


    t = np.array(tt[0:-1])
    t_new = np.array(tt[1:])

    dt_vals = t_new - t
    
    Results = Output(u=u, Eu=E, tt=tt, dt=dt_vals)

    return Results
