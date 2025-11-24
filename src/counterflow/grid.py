"""
grid.py — Grid creation utilities for the counterflow combustion simulation.

This module defines:
- a uniform 2D Cartesian grid
- helper functions to allocate fields (u, v, p, species)
- indexing helpers if needed

All numerical solvers (Navier–Stokes, transport, combustion)
use the grid definitions from this module.
"""

import numpy as np


def create_grid(Lx: float, Ly: float, Nx: int, Ny: int):
    """
    Create a 2D uniform Cartesian grid.

    Parameters :
    Lx : float
        Domain length in x-direction (meters).
    Ly : float
        Domain length in y-direction (meters).
    Nx : int
        Number of grid points in x-direction.
    Ny : int
        Number of grid points in y-direction.

    Returns :
    x, y, dx, dy : arrays and floats
        x and y are 2D arrays of grid coordinates.
        dx and dy are the grid spacings.
    """
    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)

    x_vec = np.linspace(0, Lx, Nx)
    y_vec = np.linspace(0, Ly, Ny)
    x, y = np.meshgrid(x_vec, y_vec, indexing='ij')

    return x, y, dx, dy


def allocate_field(Nx: int, Ny: int, value: float = 0.0):
    """
    Allocate a 2D scalar field (e.g., pressure, temperature, species).

    Parameters :
    Nx, Ny : int
        Dimensions of the field.
    value : float, optional
        Initial value of the field, by default 0.0.

    Returns :
    np.ndarray
        A (Nx, Ny) array filled with `value`.
    """
    return value * np.ones((Nx, Ny), dtype=float)

def allocate_velocity_fields(Nx: int, Ny: int):
    """
    Allocate the velocity fields u and v.

    Parameters :
    Nx, Ny : int
        Grid dimensions.

    Returns :
    u, v : np.ndarray
        Two (Nx, Ny) arrays initialized to zero.
    """
    u = np.zeros((Nx, Ny), dtype=float)
    v = np.zeros((Nx, Ny), dtype=float)
    return u, v


