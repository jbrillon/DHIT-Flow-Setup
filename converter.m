function converter()

fprintf("Importing data... ")
A   = importdata('velocity.fld');
fprintf("done.\n")

% Coordinates
B   = A(:,1); % x
B   = [B A(:,2)]; % y
B   = [B A(:,3)]; % z

% Fluctuating velocity components
C   = A(:,4); % u
C   = [C A(:,5)]; % v
C   = [C A(:,6)]; % w

% Compute the dot product of the velocity fluctation vector at each point
for i=1:max(size(C))
    C(i,4) = (C(i,1)^2+C(i,2)^2+C(i,3)^2);
end

turbulent_kinetic_energy = 0.5*mean(C(:,4)); % k
rms_turbulence_intensity = sqrt((2.0/3.0)*turbulent_kinetic_energy) % u' or u_{r.m.s.}

fprintf("Turbulent kinetic energy (TKE), k = %.6f\n",turbulent_kinetic_energy)
fprintf("r.m.s. turbulence intensity, u_{r.m.s.} = %.6f\n",rms_turbulence_intensity)


% From experiment
freestream_velocity_from_experiment = 10.0; % [m/s]
rms_turbulence_intensity_from_experiment = 0.222; % [m/s] 22.2cm/s
mesh_size_from_experiment = 0.0508; % [m] 5.08cm
mesh_based_reynolds_number = 34000.0;
kinematic_viscosity_experiment = (freestream_velocity_from_experiment*mesh_size_from_experiment)/mesh_based_reynolds_number;

taylor_microscale_from_experiment = 0.00484; % [m]
integral_lengthscale_from_experiment = 0.0127; % [m]
reynolds_number_based_on_taylor_microscale_from_exp = 71.6; %27.3*(integral_lengthscale_from_experiment/taylor_microscale_from_experiment);
% kinematic_viscosity_experiment = (rms_turbulence_intensity_from_experiment*taylor_microscale_from_experiment)/reynolds_number_based_on_taylor_microscale_from_exp

% Reference values
size_of_computational_box = 11.0*mesh_size_from_experiment; % L=11M
reference_length = size_of_computational_box/(2.0*pi);
reference_velocity = sqrt(3.0/2.0)*rms_turbulence_intensity_from_experiment;
reference_time = reference_length/reference_velocity;
reference_reynolds_number = (reference_velocity*reference_length)/kinematic_viscosity_experiment;

% Nondimensional spectra output times for the computational model
nondim_times_experiment = [42.0 98.0 171.0] - 42.0*ones(3); % normalized to t*=42
dim_times_experiment = (mesh_size_from_experiment/freestream_velocity_from_experiment)*nondim_times_experiment;
nondim_times_cfd = dim_times_experiment/reference_time
% TO DO: Write these to the setup_more.dat


M = 0.5*mean(C(:,4));

disp('Mach Number')
V0 = sqrt(2*M)
M  = V0/sqrt(40*1.4)

disp('Turnover Time')
tao = 1/(V0*4)

fid = fopen('setup.dat','w');
N = max(size(A))
fprintf(fid,'%d\n',max(size(A)));
for k=1:N
    fprintf(fid,'%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n',B(k,1),B(k,2),B(k,3),C(k,1),C(k,2),C(k,3));
end
fclose(fid)

% Write additional setup details
fid2 = fopen('setup_more.dat','w');
fprintf(fid2,"Mean velocity\n");
fprintf(fid2,'%18.16e\n',V0);
fprintf(fid2,"Mach number\n");
fprintf(fid2,'%18.16e\n',M);
fprintf(fid2,"Eddy turnover time\n");
fprintf(fid2,'%18.16e\n',tao);
fclose(fid2)

end

