import numpy as np

# From experiment
freestream_velocity_from_experiment = 10.0; # [m/s]
rms_turbulence_intensity_from_experiment = 0.222; # [m/s] 22.2cm/s
mesh_size_from_experiment = 0.0508; # [m] 5.08cm

# viscosity based on freestream mesh-based Reynolds number
mesh_based_reynolds_number = 34000.0;
kinematic_viscosity_experiment = (freestream_velocity_from_experiment*mesh_size_from_experiment)/mesh_based_reynolds_number;

# # viscosity based on the t=42 table values
# taylor_microscale_from_experiment = 0.00484; # [m]
# integral_lengthscale_from_experiment = 0.0127; # [m]
# reynolds_number_based_on_taylor_microscale_from_exp = 71.6; #27.3*(integral_lengthscale_from_experiment/taylor_microscale_from_experiment);
# kinematic_viscosity_experiment = (rms_turbulence_intensity_from_experiment*taylor_microscale_from_experiment)/reynolds_number_based_on_taylor_microscale_from_exp

print("Viscosity: %18.6e" % kinematic_viscosity_experiment)

reynolds_number_based_on_mesh_and_rms_velocity_from_experiment = (rms_turbulence_intensity_from_experiment*mesh_size_from_experiment)/kinematic_viscosity_experiment
print("Mesh based Reynolds number using r.m.s. velocity: %18.6f" % reynolds_number_based_on_mesh_and_rms_velocity_from_experiment)

# Reference values
# kinematic_viscosity_computation=6.206e-4;
kinematic_viscosity_computation=1.0*kinematic_viscosity_experiment;
size_of_computational_box = 11.0*mesh_size_from_experiment; # L=11M
reference_length = size_of_computational_box/(2.0*np.pi);
# reference_length = (2.0*np.pi)/size_of_computational_box;
reference_velocity = np.sqrt(3.0/2.0)*rms_turbulence_intensity_from_experiment;
reference_time = reference_length/reference_velocity;
reference_reynolds_number = (reference_velocity*reference_length)/kinematic_viscosity_computation;

print("Simulation Reynolds number: %18.6f" % reference_reynolds_number)


# # Nondimensional spectra output times for the computational model
# nondim_times_experiment = [42.0 98.0 171.0] - 42.0*ones(3); # normalized to t*=42
# dim_times_experiment = (mesh_size_from_experiment/freestream_velocity_from_experiment)*nondim_times_experiment;
# nondim_times_cfd = dim_times_experiment/reference_time
# # TO DO: Write these to the setup_more.dat
