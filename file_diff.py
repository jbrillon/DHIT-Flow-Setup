import numpy as np

nRows=13824

base_path = "/home/julien/Codes/2022-09-01/PHiLiP/build_release/tests/integration_tests_control_files/decaying_homogeneous_isotropic_turbulence/"
file_new = np.loadtxt(base_path+"velocity_vorticity-0-00000.dat",skiprows=1,usecols=(0,1,2,3,4,5),max_rows=nRows,dtype=np.float64)

# base_path = "/home/julien/Codes/PHiLiP/build_release/tests/integration_tests_control_files/decaying_homogeneous_isotropic_turbulence/"
# file_original = np.loadtxt(base_path+"velocity_vorticity-0-00000.dat",skiprows=1,max_rows=nRows,dtype=np.float64)
# base_path = "/home/julien/Codes/PHiLiP/build_release/tests/integration_tests_control_files/decaying_homogeneous_isotropic_turbulence/"
# file_original = np.loadtxt(base_path+"velocity_vorticity-0-00000.dat",skiprows=1,max_rows=nRows,dtype=np.float64)

base_path = "/home/julien/Codes/2022-09-01/PHiLiP/build_release/tests/integration_tests_control_files/decaying_homogeneous_isotropic_turbulence/"
file_original = np.loadtxt(base_path+"velocity-0-00000-equid-new.dat",skiprows=1,max_rows=nRows,dtype=np.float64)

for i in range(0,nRows):
    err = np.abs(np.linalg.norm(file_original[i,:]-file_new[i,:]))
    if(err>3.0e-15):
        print("Error: %1.6e" % err)
