"""
Example 03 — Run homogeneous reactor or full combustion
"""

import numpy as np
from counterflow.config import default_params
from counterflow.combustion import reactor

# --- Params
params = default_params()
Nx = params["Nx"]
Nt = params["Nt"]
dt = params["dt"]
rho = params["rho"]
cp  = params["cp"]
D   = params["D"]
Lx  = params["Lx"]
dx  = Lx / (Nx - 1)

# --- Load velocity field
u = np.load("u.npy")
v = np.load("v.npy")

# --- Load steady state (produced via analysis)
Nstat = 1200  # example

# --- Run combustion model
YCH4, YO2, YH2O, YCO2, T, Ntchem, times, Temp_max
