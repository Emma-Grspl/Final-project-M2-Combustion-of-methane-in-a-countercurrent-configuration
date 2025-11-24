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

    Parameters
    ----------
    Lx : float
        Domain length in x-direction (meters).
    Ly : float
        Domain length in y-direction (meters).
    Nx : int
        Number of grid points in x-direction.
    Ny : int
        Number of grid points in y-direction.

    Returns
    -------
    x, y, dx, dy : arrays and floats
        x and y are 2D arrays of grid coordinates.
        dx and dy are the grid spa

