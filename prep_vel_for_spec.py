import numpy as np
from var import *

# Filename for resulting averaged velocity field (output of this code)
averaged_velocity_field_filename = "velocity.fld"

print("Input nDOF: %i" % nDOF)
print("Resulting reduced nDOF: %i" % reduced_nDOF)

print(" ")
# Specify velocity field file
velocity_field_filename = input('Input velocity field filename: ')
print("Averaging velocity field given by file: %s" % velocity_field_filename)

all_coordinates = np.loadtxt(velocity_field_filename,skiprows=0,usecols=(0,1,2),dtype=np.float64)
all_velocities = np.loadtxt(velocity_field_filename,skiprows=0,usecols=(3,4,5),dtype=np.float64)
unique_coordinates = -1.0*np.ones((reduced_nDOF,3),dtype=np.float64)
averaged_velocities = np.zeros((reduced_nDOF,3),dtype=np.float64)
non_averaged_velocities = np.zeros((reduced_nDOF,3),dtype=np.float64) # for verification purposes
number_of_points_to_average_with = np.ones(reduced_nDOF,dtype=np.float64)

print("... averaging ...")

j=0
for i in range(0,nDOF):
    check = np.equal(unique_coordinates,all_coordinates[i,:]).all(1)
    if(any(check)):
        index_of_repeated_point = np.where(check)[0][0]
        averaged_velocities[index_of_repeated_point,:] += all_velocities[i,:]
        number_of_points_to_average_with[index_of_repeated_point] += 1.0
        continue
    else:
        unique_coordinates[j,:] = all_coordinates[i,:]
        averaged_velocities[j,:] = 1.0*all_velocities[i,:]
        non_averaged_velocities[j,:] = 1.0*all_velocities[i,:]
        j += 1
# average the values
for j in range(0,reduced_nDOF):
    averaged_velocities[j,:] = averaged_velocities[j,:]/number_of_points_to_average_with[j]

print("done.")

# UNIT TEST:
print(" ")
print("Performing unit test for checking coordinates and averaged values...")
# Check that the coordinates match and averaged velocities match the expected values
expected_unique_coordinates = np.loadtxt("expected_input.fld",skiprows=0,usecols=(0,1,2),dtype=np.float64)
for i in range(0,reduced_nDOF):
    err_coords = np.linalg.norm(unique_coordinates[j,:]- expected_unique_coordinates[j,:])
    err_averaging = np.linalg.norm(non_averaged_velocities[j,:]-averaged_velocities[j,:])
    if((err_coords > 1.0e-16) or (err_averaging > 1.0e-16)):
        if(err_coords > 1.0e-16):
            print("Error: Coordinates do not match within tolerance")
            print("Index: %i | Error val: %18.16e"  % (i,err_averaging))
            print("Aborting...")
            exit()
        if(err_averaging > 1.0e-16):
            print("Error: Averaging do not match within tolerance")
            print("Index: %i | Error val: %18.16e"  % (i,err_averaging))
            print("Aborting...")
            exit()
print("passed.")
# END OF UNIT TEST

# UNIT TEST:
# file1 = open("unique_coords_only.dat","w")
# file2 = open("expected_unique_coords_only.dat","w")
# for i in range(0,reduced_nDOF):
#     wstr1 = " %21.18f %21.18f %21.18f\n" % \
#         (unique_coordinates[i,0],unique_coordinates[i,1],unique_coordinates[i,2])
#     wstr2 = " %21.18f %21.18f %21.18f\n" % \
#         (expected_unique_coordinates[i,0],expected_unique_coordinates[i,1],expected_unique_coordinates[i,2])
#     file1.write(wstr1)
#     file2.write(wstr2)
# file1.close()
# file2.close()
# Note: To test --> diff expected_unique_coords_only.dat unique_coords_only.dat
# END OF UNIT TEST

# write the input file for the fortran code
print(" ")
print("Writing averaged velocity field to file: %s" % averaged_velocity_field_filename)
file = open(averaged_velocity_field_filename,"w")
print("... writing ...")
for i in range(0,reduced_nDOF):
    wstr = " %21.18f %21.18f %21.18f %21.18f %21.18f %21.18f\n" % \
        (unique_coordinates[i,0],unique_coordinates[i,1],unique_coordinates[i,2],\
            averaged_velocities[i,0],averaged_velocities[i,1],averaged_velocities[i,2])
    file.write(wstr)
file.close()
print("done.")

print(" ")
print("%s is ready to be passed to box.for code to extract spectra." % averaged_velocity_field_filename)