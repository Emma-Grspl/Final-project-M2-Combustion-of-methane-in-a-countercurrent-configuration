import numpy as np
from counterflow.config import default_params
params = default_params()
Nx = params["Nx"]

def k(i, j):
    """
    Convert 2D (i,j) grid indices into a 1D flattened index.

    Parameters
    i : int
        Index along x.
    j : int
        Index along y.

    Returns
    int
    Flattened index corresponding to (i, j)
    """
    return i + j*Nx
    
# Steady-state detection 
def steady_time(Nx, Nt, p, dt, eps=1e-6):
    """
    Determine the iteration and physical time at which the pressure field reaches steady state. Steady state is evaluated by monitoring the temporal evolution of the sum of all pressure values. When the relative change between two  consecutive time steps becomes smaller than `eps`, the flow is considered steady.

    Parameters
    Nx : int
        Number of grid points in each spatial direction (square grid Nx × Nx).
    Nt : int
        Total number of time steps.
    p : ndarray of shape (Nx*Nx, Nt)
        Pressure field history computed by the fractional-step solver.
    dt : float
        Hydrodynamic time step size.
    eps : float, optional
        Relative tolerance for steady-state detection. Default is 1e-6.

    Returns
    n : int
        Time step index at which steady state is detected.
    tstat : float
        Physical time corresponding to that iteration (n * dt).

    Notes
    - If no steady state is detected before the final iteration, the function returns the last time step (Nt-1) and its associated time.
    - This criterion is global but robust: it only checks the total pressure content, not local variations.
    """

    Ysum = Nx*Nx

    for n in range(1, Nt):
        S = np.sum(p[:, n])
        diff = abs(S - Ysum) / S
        Ysum = S

        if diff < eps:
            return n, n*dt

    return Nt-1, (Nt-1)*dt

def compute_shear_evolution(v, Nx, Nt, dx, Nstat):
    """
    Compute dv/dy on the left boundary, its maximum at each time, and the global maxima (absolute + steady-state).

    Parameters
    v : ndarray (Nx*Nx, Nt)
        Vertical velocity field stored as flattened arrays over time.
    k : function
        Index mapping k(i,j,Nx).
    Nx : int
        Grid size.
    Nt : int
        Number of time steps.
    dx : float
        Spatial resolution.
    Nstat : int
        Time index representing steady state.

    Returns
    dv_dy : ndarray (Nx, Nt)
        Shear rate dv/dy on the left boundary as a function of y and t.

    dv_dy_max : ndarray (Nt,)
        Maximum shear rate over y for each time.

    abs_max_value : float
        Absolute maximum value of dv/dy in all space-time.

    abs_max_position : float
        y-position (meters) where the absolute maximum occurs.

    stat_max_value : float
        Maximum dv/dy at steady state.

    stat_max_position : float
        y-position (meters) of that steady-state maximum.
    """

    dv_dy = np.zeros((Nx, Nt))

    for t in range(Nt):
        for j in range(Nx):
            if j == 0:
                dv_dy[j, t] = abs((v[k(0, j+1), t] - v[k(0, j), t]) / dx)
            elif j == Nx - 1:
                dv_dy[j, t] = abs((v[k(0, j), t] - v[k(0, j-1), t]) / dx)
            else:
                dv_dy[j, t] = abs((v[k(0, j+1), t] - v[k(0, j-1), t]) / (2*dx))

    dv_dy_max = dv_dy.max(axis=0)
    abs_max_value = dv_dy.max()
    j_abs, t_abs = np.unravel_index(np.argmax(dv_dy), dv_dy.shape)
    abs_max_position = j_abs * dx
    dv_dy_stat_profile = dv_dy[:, Nstat]
    stat_max_value = dv_dy_stat_profile.max()
    j_stat = np.argmax(dv_dy_stat_profile)
    stat_max_position = j_stat * dx

    return (dv_dy,dv_dy_max,abs_max_value,abs_max_position,stat_max_value,stat_max_position)

import numpy as np
import matplotlib.pyplot as plt

def compute_diffusive_thickness(YN2_calc, Nx, dx,Ymin=0.1*0.79,Ymax=0.9*0.79,plot=True):
    """
    Compute the thickness of the diffusive zone on the left boundary based on Y_N2 at final time.

    Parameters
    YN2_calc : array of shape (Nx*Nx, 4)
        Stored Y_N2 values at 4 selected times.
    Nx : int
        Grid resolution.
    dx : float
        Spatial grid size.
    Ymin, Ymax : float
        Thresholds defining the diffusive zone.
    plot : bool
        If True, generate the YN2 profile plot.

    Returns
    thickness : float
        Physical thickness of the diffusive zone (meters).
    idx_start : int
        Index where diffusive zone begins on the wall.
    """

    YN2_mod = np.reshape(YN2_calc, (Nx, Nx, 4))
    Y_left = YN2_mod[:, 0, 3]
    mask = (Y_left > Ymin) & (Y_left < Ymax)
    count_cells = np.sum(mask)
    thickness = count_cells * dx
    if count_cells > 0:
        idx_start = Nx - count_cells
    else:
        idx_start = None

    if plot:
        h = np.linspace(0, Nx+1, Nx)

        plt.plot(h, Y_left, "r")
        plt.grid(True)
        plt.xlabel("y (grid index)")
        plt.ylabel("YN2 profile")
        plt.title("YN2 profile on the left wall")
        plt.ylim(0.4, 0.9)

        if idx_start is not None:
            plt.axvline(x=idx_start, linestyle="--", label="diffusive zone")
        plt.axvline(x=Nx, linestyle="--")

        plt.legend(loc="lower left")
        plt.show()

    return thickness, idx_start
