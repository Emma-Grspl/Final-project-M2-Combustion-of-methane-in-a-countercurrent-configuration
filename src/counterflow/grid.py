import numpy as np

def create_grid(Lx: float, Nx: int):
    """
    Create a 2D uniform Cartesian grid.

    Parameters :
    Lx : float
        Domain length in x-direction (meters).
    Nx : int
        Number of grid points in x-direction.

    Returns :
    x, y, dx, dy : arrays and floats
        x and y are 2D arrays of grid coordinates.
        dx and dy are the grid spacings.
    """
    dx = Lx / (Nx - 1)
    dy = Lx / (Nx - 1)

    x_vec = np.linspace(0, Lx, Nx)
    y_vec = np.linspace(0, Lx, Nx)
    x, y = np.meshgrid(x_vec, y_vec, indexing='ij')

    return x, y, dx, dy


def allocate_field(Nx: int, value: float = 0.0):
    """
    Allocate a 2D scalar field (e.g., pressure, temperature, species).

    Parameters :
    Nx : int
        Dimensions of the field.
    value : float, optional
        Initial value of the field, by default 0.0.

    Returns :
    np.ndarray
        A (Nx, Nx) array filled with `value`.
    """
    return value * np.ones((Nx, Nx), dtype=float)

def allocate_velocity_fields(Nx: int):
    """
    Allocate the velocity fields u and v.

    Parameters :
    Nx: int
        Grid dimensions.

    Returns :
    u, v : np.ndarray
        Two (Nx, Nx) arrays initialized to zero.
    """
    u = np.zeros((Nx, Nx), dtype=float)
    v = np.zeros((Nx, Nx), dtype=float)
    return u, v


