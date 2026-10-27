# Define classes needed for model setup

from dataclasses import dataclass
import numpy as np
from typing import Optional

@dataclass
class Parameters:
    model: str
    eps: float
    m: float

    
    alpha: Optional[float] = None
    beta: Optional[float] = None
    delta: Optional[float] = None
    err_tol: Optional[np.ndarray] = None
    S: Optional[float] = None
    B: Optional[float] = None   
    eta1: float = None
    eta2: float = None
    tau: float = 0.0

    def __post_init__(self):

        if self.alpha is None:
            self.alpha = 0.0

        if self.beta is None:
            self.beta = 1.0

        if self.delta is None:
            self.delta = 1e-8

        if self.err_tol is None:
            self.err_tol = np.array([1e-4, 1e-4])

        if self.B is None:
            self.B = 1.0

        if self.S is None:
            self.S = 2.0
            
        if self.eta1 is None:
            self.eta1 = self.eps**2

        if self.eta2 is None:
            self.eta2 = self.eta1

@dataclass
class TimeDomain:
    dt_min: float
    dt_max: float
    t0: float
    tf: float

@dataclass
class Grid:

    L: np.ndarray
    N: np.ndarray
    d: np.ndarray

    V: Optional[float] = None
    dV: Optional[float] = None

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

    grid.V = np.prod(L)
    grid.dV = np.prod(d)

    x = []
    k = []
    k_even = []

    for i in range(dim):
        xi = L[i] * (np.arange(1, N[i] + 1)) / N[i] - L[i]/2

        ki = np.concatenate((
            np.arange(0, N[i] // 2),
            np.array([0.0]),
            np.arange(-N[i] // 2 + 1, 0)
        )) / (L[i] / np.pi / 2)

        k_even_i = np.concatenate([
            np.arange(0, N[i] // 2),
            [N[i] // 2],
            np.arange(-N[i] // 2 + 1, 0)
        ]).reshape(-1, 1) / (L[i] / np.pi / 2)

        x.append(xi)
        k.append(ki)
        k_even.append(k_even_i)

    if dim == 1:

        grid.x = x[0]
        grid.k = k[0]
        k2_full = k_even[0]**2

    elif dim == 2:
        grid.x = x

        grid.xx, grid.yy = np.meshgrid(x[0], x[1], indexing="xy")

        grid.kxx, grid.kyy = np.meshgrid(k[0], k[1], indexing="xy")

        kxx_even, kyy_even = np.meshgrid(k_even[0], k_even[1], indexing="xy")

        grid.k = np.sqrt(grid.kxx**2 + grid.kyy**2)
        k2_full = kxx_even**2 + kyy_even**2

    elif dim == 3:
        grid.x = x

        grid.xx, grid.yy, grid.zz = np.meshgrid(x[0], x[1], x[2], indexing="ij")

        grid.kxx, grid.kyy, grid.kzz = np.meshgrid(k[0], k[1], k[2], indexing="ij")

        kxx_even, kyy_even, kzz_even = np.meshgrid(k_even[0], k_even[1], k_even[2], 
                                                   indexing="xy")

        grid.k = np.sqrt(grid.kxx**2 + grid.kyy**2 + grid.kzz**2)
        k2_full = kxx_even**2 + kyy_even**2 + kzz_even**2

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
