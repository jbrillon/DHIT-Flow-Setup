import numpy as np
from scipy.optimize import fsolve

# ================================================================
# Determine viscosity, temperature, Prandtl number, and pressure
# ================================================================
def sutherlands_law(T):
    # Sutherland's law for viscosity
    # Reference: https://www.cfd-online.com/Wiki/Sutherland%27s_law
    S = 110.4
    mu_ref = 1.716e-5
    T_ref = 273.15
    return mu_ref*((T/T_ref)**(1.5))*((T_ref+S)/(T+S))

def func(T):
    global mu
    return mu - sutherlands_law(T)

global mu
mean_velocity = 10.0 # m/s
Re_mesh = 34000.0
mesh_size = 5.08/100.0 # m
density = 1.225 # kg/m3
viscosity = density*mean_velocity*mesh_size/Re_mesh
mu = 1.0*viscosity

print(" ")
print("----------------------------------------------------------------------------")
print("Determined temperature, pressure, and Prandtl number from experiment:")
print("----------------------------------------------------------------------------")

# the dimensional reference (i.e. freestream temperature)
temperature = fsolve(func,298.15)
print(" - Freestream temperature: %3.13f [K]" % temperature)

# the dimensional pressure we initialize the flow with
R_gas = 287.057 # J/kg/K
pressure = density*R_gas*temperature
print(" - Uniform pressure: %3.11f [Pa]" % pressure)
print(" - Uniform pressure: %3.11f [bar]" % (pressure/100000.0))

# Prandtl number obtained by interpolating values from Engineering Toolbox at 1bar 
# https://www.engineeringtoolbox.com/air-prandtl-number-viscosity-heat-capacity-thermal-conductivity-d_2009.html
prandtl_number = ((0.707-0.709)/(300.0-289.0))*(temperature-289.0) + 0.709
print(" - Prandtl number: %0.5f" % prandtl_number)

# ================================================================
# Determine similarity parameters:
# ================================================================
# (1) Set reference values:
density_ref = 1.225 # kg
length_ref = 11.0*5.08/(2.0*np.pi) # cm
velocity_ref = np.sqrt(3.0/2.0)*22.2 # cm/s
time_ref = length_ref/velocity_ref # [s]
length_ref /= 100.0 # [m]
velocity_ref /= 100.0 # [m/s]
temperature_ref = 1.0*temperature #[K]
viscosity_ref = 1.0*viscosity
pressure_ref = density_ref*velocity_ref*velocity_ref
reynolds_number_ref = density_ref*velocity_ref*length_ref/viscosity_ref
gamma_gas = 1.4
mach_number_ref = velocity_ref/np.sqrt(gamma_gas*R_gas*temperature_ref)
print(" ")
print("----------------------------------------------------------------------------")
print("Similarity parameters for the nondimensionalized Navier-Stokes equations:")
print("----------------------------------------------------------------------------")
print(" - Reynolds number reference: %4.12f" % reynolds_number_ref)
print(" - Mach number reference: %0.16f" % mach_number_ref)
# check that mesh based Reynolds number using reference values is the same as the experiment
print(" * Mesh based Reynolds number using reference values (must be same as experiment): %f" % (density_ref*mean_velocity*mesh_size/viscosity_ref))

# (2) Determine the constant nondimensional values used for initializing the flow
nondim_mean_velocity = mean_velocity/velocity_ref
nondim_pressure = pressure/pressure_ref
nondim_density = density/density_ref
print(" ")
print("----------------------------------------------------------------------------")
print("Nondimensionalized quantities for initializing the flow:")
print("----------------------------------------------------------------------------")
print(" - Nondimensionalized density: %f" % nondim_density)
print(" - Nondimensionalized mean velocity: %f" % nondim_mean_velocity)
print(" - Nondimensionalized pressure: %f" % nondim_pressure)

# (3) write all the values to an output file
filename = "parameters_for_dhit_setup.txt"
print(" ")
print("----------------------------------------------------------------------------")
print("Writing parameters to file: %s" % filename)
print("----------------------------------------------------------------------------")
print("Order of the parameters is listed below: ")
print("temperature, prandtl_number, reynolds_number_ref, mach_number_ref, nondim_density, nondim_mean_velocity, nondim_pressure")
n_values = 7
values_for_setup = np.zeros((1,n_values),dtype=np.float64)
values_for_setup[0,0] = temperature
values_for_setup[0,1] = prandtl_number
values_for_setup[0,2] = reynolds_number_ref
values_for_setup[0,3] = mach_number_ref
values_for_setup[0,4] = nondim_density
values_for_setup[0,5] = nondim_mean_velocity
values_for_setup[0,6] = nondim_pressure
np.savetxt(filename,values_for_setup)
