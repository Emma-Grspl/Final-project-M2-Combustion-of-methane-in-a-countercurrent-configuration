"""
Example 02 — Compute N2 transport using precomputed velocity field
"""

import numpy as np
from counterflow.config import default_params
from counterflow.transport import calcul_YN2, YN2_0

# --- Load params
params = default_params()
Nx = params["Nx"]
Nt = params["Nt"]
D  = params["D"]
dt = params["dt"]
Lx = params["Lx"]

dx = Lx / (Nx - 1)

# --- Load velocity field (from example 01)
u = np.load("u.npy")
v = np.load("v.npy")

# --- Compute YN2
YN2_calc, times = calcul_YN2(u, v, YN2_0, dt, dx, D, Nt, Nx, rho=1.0)

# --- Save
np.save("YN2.npy", YN2_calc)

print("Transport example completed.")
