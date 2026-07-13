from dataclasses import dataclass
import numpy as np
from typing import Optional

@dataclass
class Parameters:
    eps: float
    M: float
    m: float
    S: float
    B: float

@dataclass
class TimeDomain:
    #dt_min: float
    #dt_max: float
    dt: float
    t0: float
    tf: float

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