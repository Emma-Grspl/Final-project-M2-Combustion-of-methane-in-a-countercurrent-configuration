import numpy as np
from counterflow.navier_stokes import V0_cond, fractional_step

def test_V0_cond(small_params):
    Nx = small_params["Nx"]
    V0 = V0_cond(Nx)
    assert V0.shape == (Nx*Nx,)

def test_fractional_step_runs(small_params):
    Nx = small_params["Nx"]
    Nt = small_params["Nt"]
    dt = small_params["dt"]
    nu = small_params["nu"]
    rho = small_params["rho"]

    U0 = np.zeros(Nx*Nx)
    V0 = np.zeros(Nx*Nx)
    P0 = np.zeros(Nx*Nx)

    p, u, v = fractional_step(U0, V0, P0, dt, small_params["Lx"]/Nx, nu, Nt, Nx, rho)

    assert p.shape == (Nx*Nx, Nt)
    assert u.shape == (Nx*Nx, Nt)
    assert v.shape == (Nx*Nx, Nt)
