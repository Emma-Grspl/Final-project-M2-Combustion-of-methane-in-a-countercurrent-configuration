import numpy as np
from counterflow.grid import create_grid, allocate_field, allocate_velocity_fields

def test_create_grid():
    Lx = 1.0
    Nx = 10
    x, y, dx, dy = create_grid(Lx, Nx)

    assert x.shape == (Nx, Nx)
    assert y.shape == (Nx, Nx)
    assert dx == Lx / (Nx-1)
    assert dy == Lx / (Nx-1)

def test_allocate_field():
    f = allocate_field(8, value=3.14)
    assert f.shape == (8, 8)
    assert np.allclose(f, 3.14)

def test_allocate_velocity_fields():
    u, v = allocate_velocity_fields(12)
    assert u.shape == (12, 12)
    assert v.shape == (12, 12)
