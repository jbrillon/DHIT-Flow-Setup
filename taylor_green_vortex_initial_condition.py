import numpy as np
'''
   Initial Condition Function: Taylor Green Vortex (uniform density)
*  Reference: (1) Gassner2016split, 
*             (2) de la Llave Plata et al. (2019). "On the performance of a high-order multiscale DG approach to LES at increasing Reynolds number."
*  These initial conditions are given in nondimensional form (free-stream as reference)
'''
def tgv_initial_condition_velocity_gradients(pos):
    x,y,z = pos
    du_dx = np.cos(x)*np.cos(y)*np.cos(z)
    dv_dy = -np.cos(x)*np.cos(y)*np.cos(z)
    dw_dz = 0.0
    return np.array([du_dx,dv_dy,dw_dz])

def tgv_initial_condition_primitive(pos):
    x,y,z = pos
    density = 1.0
    mach_inf = 0.1
    gamma = 1.4
    # primitive variables
    u = np.sin(x)*np.cos(y)*np.cos(z)
    v = -np.cos(x)*np.sin(y)*np.cos(z)
    w = 0.0
    P = (1.0/(gamma*mach_inf*mach_inf)) + (1.0/16.0)*(np.cos(2.0*x)+np.cos(2.0*y))*(np.cos(2.0*z)+2.0)
    return np.array([density,u,v,w,P])

def primitive_to_conservative(primitive_soln):
    density,u,v,w,P = primitive_soln
    gamma = 1.4
    # to conservative variables:
    xmomentum = density*u; ymomentum = density*v; zmomentum = density*w;
    energy = P/(gamma-1.0) + 0.5*density*(u*u + v*v + w*w)
    return np.array([density,xmomentum,ymomentum,zmomentum,energy])

def tgv_initial_condition(pos):
    primitive_soln = tgv_initial_condition_primitive(pos)
    conservative_soln = primitive_to_conservative(primitive_soln)
    return conservative_soln