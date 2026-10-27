# SSAV

# import libraries
import ssav_helpers
from ssav_helpers import *
import numpy as np
from numpy.fft import fftn, ifftn, fftfreq
from dataclasses import dataclass

# Create a data class to store the solution and simulation history
@dataclass
class Output:
    u: np.ndarray
    Eu: np.ndarray
    tt: np.ndarray
    dt: np.ndarray
    reject: float


def ssav(Grid, Time, Para, u, dim):
    """
    Solve the phase field model equation using the SSAV method

    Parameters
    ----------
    Grid : object
        Spatial grid, including wave numbers, volume, and grid spacing.
    Time : object
        Time-stepping settings, including final time and time-step bounds.
    Para : object
        Model parameters and numerical tolerances.
    u : ndarray
        Initial condition.
    dim : int
        Spatial dimension (1, 2 or 3).

    Returns
    -------
    Output
        Final solution, energy history, time points, and time-step.
    """

    # Initialize the time step and record the initial time
    dt = Time.dt_min
    tt = []
    tt.append(0)

    # Keep track of rejected steps
    rejected_steps = 0

    Nt = int(Time.tf / dt)

    if Para.model == "PFC":
        D = (-Grid.k2 + 1)

    else:
        D = -1j*Grid.k

    # Define linear operators L0 and L1
    # L0 is used in the 2 linear solves
    # L1 is applied to the extrapolated u when computing r
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

    #==================================================
    #======== Compute u^1 using backward Euler ========
    #==================================================
    F, F1, F2 = compute_F(u, Para.beta, Para.tau)
    G = compute_G(u, ufft, Para, F, F1, Grid, dim)

    # Define this as w_old b/c we will need it in the loop below
    w_old = compute_w(G, Para.B)

    # Initialize energy
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

    # Define linear operator
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

    E_SSAV = compute_E(2*ufft - uold_fft, 2*w - w_old,
                       Grid, Para, D, dim)

    E_sym = []
    E_sym.append(0.5*(E_curr + E_SSAV))

    t = dt
    tt.append(t)
    dt_new = dt
    
    #==================================================
    #================= Main time loop =================
    #==================================================
    while t < Time.tf:

        step_dt = dt_new

        u_new, w_new = compute_unew(u, ufft, u_old, uold_fft, w, w_old, dt, step_dt,
                                    Para, Grid, L0, L1, t, dim)

        unew_fft = fftn(u_new)

        E_new = compute_E(unew_fft, w_new, Grid, Para, D, dim)

        E_SSAV = compute_E(2*unew_fft - ufft, 2*w_new - w,
                       Grid, Para, D, dim)
        E_sym.append(0.5*(E_new + E_SSAV))
        
        E_t = (E_sym[-1] - E_sym[-2]) / step_dt

        # E_t = (E_new - E_curr) / step_dt

        # Check if energy is decreasing
        # If not, then decrease the time step and recompute the solution
        l = 0
        while E_t > 1e-8 and dt_new > Time.dt_min and l <= 5:
            
            dt_floor = min(Time.dt_min, 0.1 * dt)
            reject_steps += 1

            step_dt = max(0.5 * dt_new, dt_floor)

            gamma = dt_new / dt

            u_new, w_new = compute_unew(u, ufft, u_old, uold_fft, w, w_old, dt, step_dt,
                                    Para, Grid, L0, L1, t, dim)

            unew_fft = fftn(u_new)

            E_new = compute_E(unew_fft, w_new, Grid, Para, D, dim)

            E_SSAV = compute_E(2*unew_fft - ufft, 2*w_new - w,
                       Grid, Para, D, dim)
            E_sym[-1] = 0.5*(E_new + E_SSAV)
        
            E_t = (E_sym[-1] - E_sym[-2]) / step_dt
            
            # E_t = (E_new - E_curr) / step_dt

        
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

        if (t + dt_new) > Time.tf and Time.dt_max != Time.dt_min:
            dt_new = Time.tf - t
            
            if dt_new < Time.dt_min:
                break
            


    t = np.array(tt[0:-1])
    t_new = np.array(tt[1:])

    dt_vals = t_new - t
    
    Results = Output(u=u, Eu=E, tt=tt, dt=dt_vals, reject=rejected_steps)

    return Results
