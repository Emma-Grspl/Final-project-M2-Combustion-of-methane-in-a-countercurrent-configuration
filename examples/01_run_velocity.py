"""
Example 01 — Compute the velocity and pressure fields
"""

import numpy as np
from counterflow.config import default_params
from counterflow.navier_stokes import V0_cond, fractional_step
from counterflow.grid import create_grid, allocate_velocity_fields

# Load config
params = default_params()
Nx = params["Nx"]
Ny = params["Ny"]
dt = params["dt"]
nu = params["nu"]
rho = params["rho"]
Nt = params["Nt"]
Lx = params["Lx"]
Ly = params["Ly"]

# Grid
x, y, dx, dy = create_grid(Lx, Ly, Nx, Ny)

# Initial conditions
U0, V0 = allocate_velocity_fields(Nx, Ny)
V0[:, :] = V0_cond(Nx).reshape(Nx, Ny)
P0 = np.zeros(Nx * Ny)

# Solve
p, u, v = fractional_step(U0.flatten(), V0.flatten(), P0, dt, dx, nu, Nt, Nx, rho)

# Save
np.save("u.npy", u)
np.save("v.npy", v)
np.save("p.npy", p)

print("Velocity example completed.")