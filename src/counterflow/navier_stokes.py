""" 
Fractional-step solver for the counterflow velocity field.

This file contains your exact physical method:
- index k(i,j,Nx)
- initial V0 condition
- u* computation
- pressure matrix A
- RHS b
- direct Poisson solver
- fractional step loop
- steady state detection
"""

import numpy as np


# Index mapping 
def k(i, j, Nx):
    """
    Convert 2D indices (i, j) into a 1D flattened index.

    Parameters
    i : int
    Index in the x-direction.
    j : int
    Index in the y-direction.
    Nx : int
    Number of grid points in one direction (square grid).

    Returns
    int
    Flattened index k = i + j * Nx.
    """
    return i + j * Nx


# Initial V condition 
def V0_cond(Nx):
    """
    Build the initial vertical velocity field V0 for the counterflow.
    Boundary condition:
    - Bottom (j = 0):   fast jet (i < Nx/4) = 1, slow jet (Nx/4 < i < Nx/2) = 0.2
    - Top (j = Nx-1):   fast downward jet = -1, slow = -0.2

    Parameters
    Nx : int
    Grid resolution (square domain Nx × Nx).

    Returns
    V0 : ndarray of shape (Nx*Nx,)
    Vector storing the initial vertical velocity.
    """
    V0 = np.zeros(Nx * Nx)
    for i in range(Nx):
        for j in range(Nx):
            if j == 0:  # bottom
                if 0 <= i < Nx//4:
                    V0[k(i, j, Nx)] = 1
                elif Nx//4 <= i < Nx//2:
                    V0[k(i, j, Nx)] = 0.2
            if j == Nx-1:  # top
                if 0 <= i < Nx//4:
                    V0[k(i, j, Nx)] = -1
                elif Nx//4 <= i < Nx//2:
                    V0[k(i, j, Nx)] = -0.2
    return V0


# Compute u* and v* 
def u_etoile(u0, v0, dx, dt, Nx, nu):
    """
    Compute intermediate velocities u* and v* (fractional-step method).

    The formula applies:
    u* = u - dt (u ∂u/∂x + v ∂u/∂y) + ν Δu
    v* = v - dt (u ∂v/∂x + v ∂v/∂y) + ν Δv

    Boundary conditions: values are simply copied on the boundaries
    exactly as in the original student method.

    Parameters
    u0, v0 : ndarray (Nx*Nx,)
    Current velocity field (flattened 1D vectors).
    dx : float
    Grid spacing (dx = dy assumed).
    dt : float
    Time step.
    Nx : int
    Grid resolution.
    nu : float
    Kinematic viscosity.

    Returns
    u, v : ndarray (Nx*Nx,)
    Intermediate velocity fields u* and v*.
    """
    
    u = np.empty(Nx*Nx)
    v = np.empty(Nx*Nx)

    mu = dt / (2*dx)
    mu2 = dt / dx**2

    for i in range(Nx):
        for j in range(Nx):
            idx = k(i, j, Nx)

            if i == 0 or i == Nx-1 or j == 0 or j == Nx-1:
                u[idx] = u0[idx]
                v[idx] = v0[idx]
                continue

            # advection + diffusion
            u[idx] = (u0[idx]- mu * (u0[idx] * (u0[k(i+1, j, Nx)] - u0[k(i-1, j, Nx)])+ v0[idx] * (u0[k(i, j+1, Nx)] - u0[k(i, j-1, Nx)]))
                      + nu * mu2 * (u0[k(i+1, j, Nx)] - 2*u0[idx] + u0[k(i-1, j, Nx)]+ u0[k(i, j+1, Nx)] - 2*u0[idx] + u0[k(i, j-1, Nx)]))
            v[idx] = (v0[idx]- mu * (u0[idx] * (v0[k(i+1, j, Nx)] - v0[k(i-1, j, Nx)])+ v0[idx] * (v0[k(i, j+1, Nx)] - v0[k(i, j-1, Nx)]))
                      + nu * mu2 * (v0[k(i+1, j, Nx)] - 2*v0[idx] + v0[k(i-1, j, Nx)]+ v0[k(i, j+1, Nx)] - 2*v0[idx] + v0[k(i, j-1, Nx)]))
    return u, v


# Pressure matrix A 
def A_matrix(Nx, dx):
    """
    Build the pressure Poisson matrix A (size Nx² × Nx²).
    Boundary conditions:
    - Neumann-type on top/bottom
    - Special left/right BC matching the original method.

    Parameters
    Nx : int
    Grid resolution.
    dx : float
    Cell size.

    Returns
    A : ndarray (Nx*Nx, Nx*Nx)
    Dense Poisson matrix (direct solve).
    """
    
    N = Nx*Nx
    A = np.zeros((N, N))

    for i in range(Nx):
        for j in range(Nx):
            idx = k(i, j, Nx)

            if i == 0:
                A[idx, idx] = -1/dx
                A[idx, k(i+1, j, Nx)] = 1/dx
                continue

            if i == Nx-1:
                A[idx, idx] = 1
                continue

            if j == 0:
                A[idx, idx] = -1/dx
                A[idx, k(i, j+1, Nx)] = 1/dx
                continue

            if j == Nx-1:
                A[idx, idx] = 1/dx
                A[idx, k(i, j-1, Nx)] = -1/dx
                continue

            # interior
            A[idx, idx] = -4/dx**2
            A[idx, k(i+1, j, Nx)] = 1/dx**2
            A[idx, k(i-1, j, Nx)] = 1/dx**2
            A[idx, k(i, j+1, Nx)] = 1/dx**2
            A[idx, k(i, j-1, Nx)] = 1/dx**2

    return A


# RHS b 
def b_secondmember(u_star, v_star, Nx, dx, dt):
    """
    Build the right-hand side b of the pressure Poisson equation.
    b = (1 / (2 dx dt)) * (du*/dx + dv*/dy)

    Parameters
    u_star, v_star : ndarray (Nx*Nx,)
    Intermediate velocities.
    Nx : int
    Grid dimension.
    dx : float
    Cell size.
    dt : float
    Time step.

    Returns
    b : ndarray (Nx*Nx,)
    RHS vector.
    """
    
    b = np.zeros(Nx*Nx)
    mu = 1/(2*dx*dt)

    for i in range(Nx):
        for j in range(Nx):
            idx = k(i, j, Nx)

            if i == 0 or i == Nx-1 or j == 0 or j == Nx-1:
                b[idx] = 0
                continue

            b[idx] = mu * (u_star[k(i+1, j, Nx)] - u_star[k(i-1, j, Nx)]+ v_star[k(i, j+1, Nx)] - v_star[k(i, j-1, Nx)])

    return b


# Direct Poisson solver 
def Poisson_direct(A, b):
    """
    Direct dense solver for the Poisson system.

    Parameters
    A : ndarray (N,N)
    Poisson matrix.
    b : ndarray (N,)
    RHS vector.

    Returns
    ndarray (N,)
    Solution p.
    """
    return np.linalg.solve(A, b)


# Main fractional-step loop 
def fractional_step(U0, V0, P0, dt, dx, nu, Nt, Nx, rho):
    """
    Perform the complete fractional-step algorithm for the 2D velocity field.
    Steps:
    1) Compute u*, v* (intermediate)
    2) Build RHS b
    3) Solve Ap = b
    4) Correct velocities
    5) Apply boundary conditions

    Parameters
    U0, V0, P0 : ndarray (Nx*Nx,)
    Initial velocity and pressure fields.
    dt : float
    Time step.
    dx : float
    Cell size.
    nu : float
    Viscosity.
    Nt : int
    Number of time steps.
    Nx : int
    Grid resolution.
    rho : float
    Density.

    Returns
    p, u, v : ndarray (Nx*Nx, Nt)
    Full time-evolution of pressure and velocity.
    """
    
    u = np.zeros((Nx*Nx, Nt))
    v = np.zeros((Nx*Nx, Nt))
    p = np.zeros((Nx*Nx, Nt))

    u[:, 0] = U0.copy()
    v[:, 0] = V0.copy()
    p[:, 0] = P0.copy()

    A = A_matrix(Nx, dx)
    print("Condition number:", np.linalg.cond(A))

    mu = dt/(2*dx)

    for n in range(1, Nt):

        u_star, v_star = u_etoile(u[:, n-1], v[:, n-1], dx, dt, Nx, nu)
        b = b_secondmember(u_star, v_star, Nx, dx, dt)
        p[:, n] = Poisson_direct(A, b)

        # correction step
        for i in range(1, Nx-1):
            for j in range(1, Nx-1):
                idx = k(i, j, Nx)

                u[idx, n] = u_star[idx] - mu/rho * (p[k(i+1, j, Nx), n] - p[k(i-1, j, Nx), n])
                v[idx, n] = v_star[idx] - mu/rho * (p[k(i, j+1, Nx), n] - p[k(i, j-1, Nx), n])

        # BCs
        for i in range(Nx):
            for j in range(Nx):
                idx = k(i, j, Nx)

                if i == 0:
                    u[idx, n] = 0
                    v[idx, n] = v[k(i+1, j, Nx), n]

                elif i == Nx-1:
                    u[idx, n] = u[k(i, j-1, Nx), n]
                    v[idx, n] = v[k(i, j-1, Nx), n]

                elif j == 0:
                    u[idx, n] = 0
                    v[idx, n] = 0
                    if 0 <= i < Nx//4:
                        v[idx] = 1
                    elif Nx//4 <= i < Nx//2:
                        v[idx] = 0.2

                elif j == Nx-1:
                    u[idx, n] = 0
                    v[idx, n] = 0
                    if 0 <= i < Nx//4:
                        v[idx] = -1
                    elif Nx//4 <= i < Nx//2:
                        v[idx] = -0.2

    return p, u, v

