"""
Species transport solver for the counterflow configuration.

This module provides:
- allocation of species fields
- boundary conditions for species
- explicit diffusion/advection update
"""

import numpy as np
from counterflow.config import default_params
params = default_params()
Nx = params["Nx"]

def k(i, j, Nx):
    """
    Convert 2D (i,j) grid indices into a 1D flattened index.

    Parameters
    i : int
        Index along x.
    j : int
        Index along y.
    Nx : int
    Returns
    int
    Flattened index corresponding to (i, j)
    """
    return i + j*Nx

def YN2_0(Nx):
    """
    Define the initial mass fraction field of N2 at t = 0.

    The boundary conditions correspond to:
    - Bottom: N2 = 0.79 for 0 ≤ i < Nx/4, N₂ = 1 for Nx/4 ≤ i < Nx/2
    - Top:    N2 = 1    for Nx/4 ≤ i < Nx/2

    Parameters
    Nx : int
    Number of grid points in each direction.

    Returns
    np.ndarray
    Flattened array (Nx*Nx,) containing the initial N₂ field.
    """
    
    Y0 = np.zeros((Nx*Nx))
    for i in range(Nx):      # x becomes the row index
        for j in range(Nx):  # y becomes the column index
            if j == 0:  # bottom boundary
                if 0 <= i < Nx//4:
                    Y0[k(i, j, Nx)] = 0.79
                if Nx//4 <= i < Nx//2:
                    Y0[k(i, j, Nx)] = 1
            if j == Nx-1:  # top boundary
                if Nx//4 <= i < Nx//2:
                    Y0[k(i, j, Nx)] = 1
    return Y0


def calcul_YN2(u, v, YN2_0, dt, dx, D, Nt, Nx, rho):
    """
    Compute the time evolution of the nitrogen mass fraction YN2 using a 2D advection–diffusion equation using an explicit finite-difference scheme with:
        - central differences for diffusion
        - central differences for advection
        - Dirichlet boundary conditions on bottom/top for imposed YN2
        - Neumann (zero-gradient) on left/right walls

    Parameters
    u : ndarray, shape (Nx*Nx, Nt)
        Horizontal velocity field.
    v : ndarray, shape (Nx*Nx, Nt)
        Vertical velocity field.
    YN2_0 : function
        Function that returns the initial YN2 distribution.
    dt : float
        Time step.
    dx : float
        Grid spacing (dx = dy).
    D : float
        Diffusivity of nitrogen.
    Nt : int
        Number of time iterations.
    Nx : int
        Grid size in x and y direction.
    rho : float
        Density (not used here but kept for uniform interface).

    Returns
    YN2_f : ndarray, shape (Nx*Nx, 4)
        YN2 at four key times [0, step=30, step=100, final].
    t : ndarray, shape (4,)
        Corresponding physical times for each stored field.
    """
    mu  = dt / (2*dx)
    mu2 = dt / dx**2

    YN2 = YN2_0(Nx).copy()
    YN2_f = np.zeros((Nx*Nx, 4))
    t = np.zeros(4)
    interval = [0, 30, 100, Nt-1]

    YN2_f[:, 0] = YN2.copy()
    t[0] = 0.0

    for n in range(1, Nt):
        Y_old = YN2.copy()  

        for i in range(1, Nx-1):
            for j in range(1, Nx-1):
                YN2[k(i,j,Nx)] = (Y_old[k(i,j,Nx)]- mu * ( u[k(i,j,Nx), n-1] * (Y_old[k(i+1,j,Nx)] - Y_old[k(i-1,j,Nx)]) + v[k(i,j,Nx), n-1] * (Y_old[k(i,j+1,Nx)] - Y_old[k(i,j-1,Nx)]))+ mu2 * D * (Y_old[k(i+1,j,Nx)] - 4*Y_old[k(i,j,Nx)] + Y_old[k(i-1,j,Nx)]+ Y_old[k(i,j+1,Nx)] + Y_old[k(i,j-1,Nx)]))

        # Boundary conditions
        for i in range(Nx):
            YN2[k(i,0,Nx)]    = YN2[k(i,1,Nx)]
            YN2[k(i,Nx-1,Nx)] = YN2[k(i,Nx-2,Nx)]
            if 0 <= i < Nx//4:
                YN2[k(i,0,Nx)] = 0.79
            elif Nx//4 <= i < Nx//2:
                YN2[k(i,0,Nx)]    = 1.0
                YN2[k(i,Nx-1,Nx)] = 1.0

        for j in range(Nx):
            YN2[k(0,j,Nx)]    = YN2[k(1,j,Nx)]
            YN2[k(Nx-1,j,Nx)] = YN2[k(Nx-2,j,Nx)]

        
        YN2 = np.clip(YN2, 0.0, 1.0)

        if n in interval:
            idx = interval.index(n)
            YN2_f[:, idx] = YN2
            t[idx] = n * dt

    return YN2_f, t