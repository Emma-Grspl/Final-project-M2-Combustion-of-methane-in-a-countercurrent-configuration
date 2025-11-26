<<<<<<< HEAD
"""
Creating the grid

This folder allows you to: 
- Create the grid
- Allocate fields for species, temperature, and speed
"""

import numpy as np

def create_grid(Lx, Nx):
=======
import numpy as np

def create_grid(Lx: float, Nx: int):
>>>>>>> origin/main
    """
    Create a 2D Cartesian grid.

<<<<<<< HEAD
    Parameters
    ----------
    Lx : float or tuple(float,float)
        Domain length(s).
    Nx : int or tuple(int,int)
        Number of grid points.
=======
    Parameters :
    Lx : float
        Domain length in x-direction (meters).
    Nx : int
        Number of grid points in x-direction.
>>>>>>> origin/main

    Returns
    -------
    x, y : ndarray
        2D grids of shape (Nx, Ny)
    dx, dy : float
        Grid spacings
    """

    dx = Lx / (Nx - 1)
    dy = Lx / (Nx - 1)

<<<<<<< HEAD
    # 1D coordinates
    x1 = np.linspace(0, Lx, Nx)
    y1 = np.linspace(0, Lx, Nx)

    # 2D mesh
    x, y = np.meshgrid(x1, y1)
=======
    x_vec = np.linspace(0, Lx, Nx)
    y_vec = np.linspace(0, Lx, Nx)
    x, y = np.meshgrid(x_vec, y_vec, indexing='ij')
>>>>>>> origin/main

    return x, y, dx, dy


def allocate_field(Nx: int, value: float = 0.0):
    """
    Allocate a 2D scalar field (e.g., pressure, temperature, species).
    Parameters :
    Nx : int
<<<<<<< HEAD
    Dimensions of the field.
=======
        Dimensions of the field.
>>>>>>> origin/main
    value : float, optional
    Initial value of the field, by default 0.0.

    Returns :
    np.ndarray
<<<<<<< HEAD
    A (Nx, Nx) array filled with `value`.
=======
        A (Nx, Nx) array filled with `value`.
>>>>>>> origin/main
    """
    return value * np.ones((Nx, Nx), dtype=float)

def allocate_velocity_fields(Nx: int):
    """
    Allocate the velocity fields u and v.
    Parameters :
    Nx: int
<<<<<<< HEAD
    Grid dimensions.
    Returns :
    u, v : np.ndarray
    Two (Nx, Nx) arrays initialized to zero.
=======
        Grid dimensions.

    Returns :
    u, v : np.ndarray
        Two (Nx, Nx) arrays initialized to zero.
>>>>>>> origin/main
    """
    u = np.zeros((Nx, Nx), dtype=float)
    v = np.zeros((Nx, Nx), dtype=float)
    return u, v


