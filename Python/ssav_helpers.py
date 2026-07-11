from dataclasses import dataclass
from typing import Optional
import numpy as np
from numpy.fft import fftn, ifftn, fftfreq

@dataclass
class Grid:

    L: np.ndarray
    N: np.ndarray
    d: np.ndarray

    x: Optional[list] = None
    k: Optional[np.ndarray] = None

    xx: Optional[np.ndarray] = None
    yy: Optional[np.ndarray] = None
    zz: Optional[np.ndarray] = None

    kxx: Optional[np.ndarray] = None
    kyy: Optional[np.ndarray] = None
    kzz: Optional[np.ndarray] = None

    inv_k: Optional[np.ndarray] = None
    k2: Optional[np.ndarray] = None
    k4: Optional[np.ndarray] = None
    k6: Optional[np.ndarray] = None


def generate_Grid(L, N, dim):

    # Generate spatial and Fourier grids

    L = np.asarray(L, dtype=float)
    N = np.asarray(N, dtype=int)

    d = L / N

    grid = Grid(L=L, N=N, d=d)

    x = []
    k = []

    for i in range(dim):
        xi = L[i] * (np.arange(1, N[i] + 1)) / N[i] - L[i]/2

        ki = np.concatenate((
            np.arange(0, N[i] // 2),
            np.array([0.0]),
            np.arange(-N[i] // 2 + 1, 0)
        )) / (L[i] / np.pi / 2)

        x.append(xi)
        k.append(ki)

    if dim == 1:

        grid.x = x[0]
        grid.k = k[0]

    elif dim == 2:
        grid.x = x

        grid.xx, grid.yy = np.meshgrid(x[0], x[1], indexing="xy")

        grid.kxx, grid.kyy = np.meshgrid(k[0], k[1], indexing="xy")

        grid.k = np.sqrt(grid.kxx**2 + grid.kyy**2)

    elif dim == 3:
        grid.x = x

        grid.xx, grid.yy, grid.zz = np.meshgrid(x[0], x[1], x[2], indexing="ij")

        grid.kxx, grid.kyy, grid.kzz = np.meshgrid(k[0], k[1], k[2], indexing="ij")

        grid.k = np.sqrt(grid.kxx**2 + grid.kyy**2 + grid.kzz**2)

    else:
        raise ValueError("dimension must be 1, 2, or 3")


    grid.inv_k = np.zeros_like(grid.k, dtype=float)
    mask = grid.k != 0
    grid.inv_k[mask] = 1.0 / (grid.k[mask]**2)
    grid.inv_k[~mask] = 1.0

    grid.k2 = grid.k**2
    grid.k4 = grid.k2**2
    grid.k6 = grid.k2**3

    return grid


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