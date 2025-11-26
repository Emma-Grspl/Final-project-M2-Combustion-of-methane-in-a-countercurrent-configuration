import numpy as np
from counterflow.transport import YN2_0, calcul_YN2

def test_YN2_0(small_params):
    Nx = small_params["Nx"]
    Y = YN2_0(Nx)
    assert Y.shape == (Nx*Nx,)

def test_calcul_YN2(small_params):
    Nx = small_params["Nx"]
    Nt = small_params["Nt"]
    dt = small_params["dt"]
    dx = small_params["Lx"]/Nx
    D = small_params["D"]
    rho = small_params["rho"]

    u = np.zeros((Nx*Nx, Nt))
    v = np.zeros((Nx*Nx, Nt))

    YN2_f, t = calcul_YN2(u, v, YN2_0, dt, dx, D, Nt, Nx, rho)

    assert YN2_f.shape == (Nx*Nx, 4)
    assert len(t) == 4
