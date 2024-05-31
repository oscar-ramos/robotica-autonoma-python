import numpy as np


def attractive_potential(x, xd, d, k=1):
    """
    x  - initial point 
    xd - final (desired) point 
    d  - threshold to change the behavior (from conic to quadratic)
    k  - scale factor (positive)
    
    """
    dist = np.linalg.norm(x-xd)
    if (dist <= d):
        # Quadratic behavior
        U = 0.5*k*dist*dist
        f = -k*(x-xd)
    else:
        # Conic behavior
        U = d*k*dist - 0.5*k*d*d
        f = -d*k*(x-xd)/dist
    return U, f


def repulsive_potential(x, xobst, d, k=1, Umax=1):
    """
        x - initial point
    xobst - point obstacle
        d - threshold to determine the "influence" zone of the obstacle
        k - scale  factor
     Umax - max value used at the border of the obstacle
    
    """
    # This considers a point obstacle at qobst
    dist = np.linalg.norm(x-xobst)
    if (dist <= d):
        if (dist == 0):
            U = 100
            f = np.zeros(2,)
        else:
            U = 0.5*k*(1.0/dist - 1.0/d)**2
            f = k*(1.0/dist - 1.0/d)*(x-xobst)/(dist**3)
        if (U > Umax):
            U = Umax
    else:
        U = 0
        f = np.zeros(2,)
    return U, f


def repulsive_potentials(x, xobst, d, k=1, Umax=1):
    """
         x - initial point
    xobst - point obstacles (more than one)
        d - threshold to determine the "influence" zone of the obstacle
        k - scale  factor
     Umax - max value used at the border of the obstacle
       
    """
    Utot = 0
    ftot = np.zeros(2,)
    for i in range(len(xobst)):
        U, f = repulsive_potential(x, xobst[i,:], d, k, Umax)
        Utot += U
        ftot += f
    return Utot, ftot
