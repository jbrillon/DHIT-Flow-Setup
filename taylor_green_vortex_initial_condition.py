import numpy as np
from scipy import integrate
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
    return np.array([du_dx,dv_dy,dw_dz],dtype=np.float64)

def tgv_initial_condition_primitive(pos,density_initialization_type="uniform"):
    x,y,z = pos

    # constants
    mach_inf = 0.1 # Freestream Mach number
    gamma = 1.4

    # density
    if(density_initialization_type=="uniform"):
        # results in a uniform density field
        density = 1.0
    elif(density_initialization_type=="isothermal"):
        # results in a uniform temperature field of T=1.0
        density = 1.0 + ((gamma*mach_inf*mach_inf)/16.0)*(np.cos(2.0*x)+np.cos(2.0*y))*(np.cos(2.0*z)+2.0)

    # primitive variables
    u = np.sin(x)*np.cos(y)*np.cos(z)
    v = -np.cos(x)*np.sin(y)*np.cos(z)
    w = 0.0
    P = (1.0/(gamma*mach_inf*mach_inf)) + (1.0/16.0)*(np.cos(2.0*x)+np.cos(2.0*y))*(np.cos(2.0*z)+2.0)
    return np.array([density,u,v,w,P],dtype=np.float64)

def primitive_to_conservative(primitive_soln):
    density,u,v,w,P = primitive_soln
    gamma = 1.4
    # to conservative variables:
    xmomentum = density*u; ymomentum = density*v; zmomentum = density*w;
    energy = P/(gamma-1.0) + 0.5*density*(u*u + v*v + w*w)
    return np.array([density,xmomentum,ymomentum,zmomentum,energy],dtype=np.float64)

def tgv_initial_condition(pos,density_initialization_type="uniform"):
    primitive_soln = tgv_initial_condition_primitive(pos,density_initialization_type)
    conservative_soln = primitive_to_conservative(primitive_soln)
    return conservative_soln

def get_kinetic_energy_from_primitive(primitive_soln):
    density,u,v,w,P = primitive_soln
    return 0.5*density*(u*u + v*v + w*w)

def get_kinetic_energy_from_zyx(z,y,x):
    pos = np.array([x,y,z])
    primitive_soln = tgv_initial_condition_primitive(pos,"isothermal") # hard-coded for isothermal
    return get_kinetic_energy_from_primitive(primitive_soln)

def get_integrated_initial_kinetic_energy():
    domain_volume = (2.0*np.pi)**3.0
    val,abserr = integrate.tplquad(get_kinetic_energy_from_zyx, -np.pi, np.pi, -np.pi, np.pi, -np.pi, np.pi, epsabs=1e-12, epsrel=1e-12)
    print("Numerical integration error: %.6e" % abserr)
    integrated_initial_kinetic_energy = val/domain_volume
    print("Integrated initial kinetic energy (non-dimensional): %1.16e" % integrated_initial_kinetic_energy)
    # print("Diff between DNS value: %1.16e" % (integrated_initial_kinetic_energy-1.2499999999999990e-01))
    # print("Diff between DNS value: %1.16e" % (integrated_initial_kinetic_energy-1.2500000000000425e-01))
    print("Diff between DNS value: %1.16e" % (integrated_initial_kinetic_energy-1.2499999999999613e-01))
    # integrate.tplquad(get_kinetic_energy_from_zyx, 0.0, (2.0*np.pi), 0.0, (2.0*np.pi), 0.0, (2.0*np.pi))
    return #integrated_initial_kinetic_energy

get_integrated_initial_kinetic_energy()

nDOF_per_dim = 21
z = np.linspace(-np.pi,np.pi,nDOF_per_dim)
y = np.linspace(-np.pi,np.pi,nDOF_per_dim)
x = np.linspace(-np.pi,np.pi,nDOF_per_dim)
for k in range(0,nDOF_per_dim):
    for j in range(0,nDOF_per_dim):
        for i in range(0,nDOF_per_dim):
            pos = np.array([x[i],y[j],z[k]])
            density,u,v,w,P = tgv_initial_condition_primitive(pos,"isothermal")
            if(density < 0.0):
                print("Density is negative")
            if(P < 0.0):
                print("Pressure is negative")