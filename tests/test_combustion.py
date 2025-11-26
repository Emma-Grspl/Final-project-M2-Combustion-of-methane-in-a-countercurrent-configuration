import numpy as np
from counterflow.combustion import Yk_0, T_0, homogeneous_reactor

def test_Yk_0(small_params):
    Nx = small_params["Nx"]
    Y = Yk_0(Nx, "CH4")
    assert Y.shape == (Nx*Nx,)

def test_T_0(small_params):
    Nx = small_params["Nx"]
    T = T_0(Nx)
    assert T.shape == (Nx*Nx,)

def test_homogeneous_reactor_runs(small_params):
    Nx = small_params["Nx"]
    Nt = small_params["Nt"]
    dt = small_params["dt"]
    rho = small_params["rho"]
    cp = small_params["cp"]
    D = small_params["D"]
    Tf = small_params["Tf"]
    Tmax = 1200

    u = np.zeros((Nx*Nx, Nt))
    v = np.zeros((Nx*Nx, Nt))

    YCH4_f, YO2_f, YH2O_f, YCO2_f, T_f, Ntchem, t = homogeneous_reactor(u, v, Nx, Tf, Nt, 2, dt, cp, rho, D, Tmax)

    assert YCH4_f.shape == (Nx*Nx, 4)
    assert YO2_f.shape == (Nx*Nx, 4)
    assert T_f.shape == (Nx*Nx, 4)
