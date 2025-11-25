import numpy as np

"""
config.py — Default simulation parameters for the counterflow combustion project.
"""

def get_config():
    """
    Return a dictionary with all default simulation parameters.

    Returns :
    dict
    Dictionary containing physical and numerical parameters.
    """
    cfg = {}

    # ---- Physical constants ----
    cfg["rho"] = 1.1614          # kg/m^3 (air)
    cfg["nu"] = 15e-6            # m^2/s kinematic viscosity
    cfg["D"] = 15e-6             # m^2/s species diffusivity
    cfg["cp"] = 1200             # J/kg/K

    # ---- Domain size ----
    cfg["Lx"] = 0.002            # m
    cfg["Ly"] = 0.002            # m

    # ---- Grid resolution ----
    cfg["Nx"] = 50
    cfg["Ny"] = 50

    # ---- Time parameters ----
    cfg["Tf"] = 0.01             # final time (s)
    cfg["Nt"] = 3000             # number of time steps

    cfg["dt"] = cfg["Tf"] / (cfg["Nt"] - 1)

    return cfg


def default_params():
    """
    Alias for get_config(), used for cleaner API.

    Returns:
    dict
    Default simulation parameters.
    """
    return get_config()
