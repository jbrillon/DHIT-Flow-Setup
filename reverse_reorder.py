import numpy as np
from var import *
#===========================================================
#            REORDER DATA FOR FORTRAN CODE
#===========================================================
#-------------------------------------------------------------
# Assemble mpi files for velocity field outputted from PHiLiP
#-------------------------------------------------------------
# get padded mpi rank string
def get_padded_mpi_rank_string(mpi_rank):
    padding_length = 5
    mpi_rank_string = '%i' % mpi_rank
    mpi_rank_string = mpi_rank_string.zfill(padding_length)
    return mpi_rank_string

subdir = "/Users/Julien/DHIT-Flow-Setup/philip_outputs/test"
prefix = "velocity_vorticity-0"

filename=subdir+"/"+prefix+".dat"
fout = open(filename, "w")
for i in range(0,num_procs):
    mpi_rank_string = get_padded_mpi_rank_string(i)
    tempfile = subdir+"/"+prefix+"-"+mpi_rank_string+".dat"
    fin = open(tempfile,"r")
    for line in fin:
        fout.write(line)
    fin.close()
fout.close()
#-------------------------------------------------------------
# Store the velocity field + coordinates
#-------------------------------------------------------------
stored_data = np.zeros((nElements_per_direction,nElements_per_direction,nElements_per_direction,nQuadPoints_per_element,nQuadPoints_per_element,nQuadPoints_per_element,1,nValues_per_row),dtype=np.float64)

velocity_file_from_philip = subdir+"/"+prefix+".dat"
fin = open(velocity_file_from_philip,"r")

# First line: Number of DOFs
nDOF_expected = int(fin.readline())
if(nDOF!=nDOF_expected):
    print("Error: nDOF does not match expected nDOF from file %s, check var.py",velocity_file_from_philip)
    print("Aborting...")
    exit()

''' must add more nested for loops for higher
    number of elements per direction
    currently can handle up to 32 (i.e. 2,4,8,16,32)
'''

iproc = 0
iDOF_per_proc = 0
start_new_file=True

ez_L_base_base_base = 0
for z_base_base_base in range(0,loop_bounds[3]):
    ey_L_base_base_base = 0
    for y_base_base_base in range(0,loop_bounds[3]):
        ex_L_base_base_base = 0
        for x_base_base_base in range(0,loop_bounds[3]):
            ez_L_base_base = ez_L_base_base_base
            for z_base_base in range(0,loop_bounds[2]):
                ey_L_base_base = ey_L_base_base_base
                for y_base_base in range(0,loop_bounds[2]):
                    ex_L_base_base = ex_L_base_base_base
                    for x_base_base in range(0,loop_bounds[2]):
                        ez_L_base = ez_L_base_base
                        for z_base in range(0,loop_bounds[1]):
                            ey_L_base = ey_L_base_base
                            for y_base in range(0,loop_bounds[1]):
                                ex_L_base = ex_L_base_base
                                for x_base in range(0,loop_bounds[1]):
                                    # algorithm for a cube with 64 (4^3) elements:
                                    ez_L = ez_L_base
                                    for cz in range(0,loop_bounds[0]):
                                        ez_R = ez_L + 1
                                        ey_L = ey_L_base
                                        for cy in range(0,loop_bounds[0]):
                                            ey_R = ey_L + 1
                                            ex_L = ex_L_base
                                            for cx in range(0,loop_bounds[0]):
                                                ex_R = ex_L + 1
                                                for ez in range(ez_L,ez_R+1):
                                                    for ey in range(ey_L,ey_R+1):
                                                        for ex in range(ex_L,ex_R+1):
                                                            for qz in range(0,nQuadPoints_per_element):
                                                                for qy in range(0,nQuadPoints_per_element):
                                                                    for qx in range(0,nQuadPoints_per_element):
                                                                        row_string = fin.readline()
                                                                        row_data = np.fromstring(row_string, dtype=np.float64, sep=' ')
                                                                        for iValue in range(0,nValues_per_row):
                                                                            stored_data[ez,ey,ex,qz,qy,qx,0,iValue] = row_data[iValue] # modify to read in vorticity
                                                ex_L += 2
                                            ey_L += 2
                                        ez_L += 2
                                    ex_L_base = ex_L
                                ey_L_base = ey_L
                            ez_L_base = ez_L
                        ex_L_base_base = ex_L_base
                    ey_L_base_base = ey_L_base
                ez_L_base_base = ez_L_base
            ex_L_base_base_base = ex_L_base_base
        ey_L_base_base_base = ey_L_base_base
    ez_L_base_base_base = ez_L_base_base
#-------------------------------------------------------------
# Store the velocity field + coordinates
#-------------------------------------------------------------
# expected coords are now cartesian
expected_coords = np.loadtxt("velocity_equidistant_nodes.fld",skiprows=0,usecols=(0,1,2),dtype=np.float64)
np.savetxt("setup_coords_only-2.dat",expected_coords,fmt="%18.16e")

reordered_velocity_file = subdir+"/"+prefix+"_reordered_for_spectra"+".dat"
# reordered_velocity_file = "reverse_order_test-2.dat"
file = open(reordered_velocity_file,"w")

# file = open("velocity.fld","w")

# file.write('Number of degrees of freedom:\n')
# wstr = "%i\n" % nDOF
# file.write(wstr)

for ez in range(0,nElements_per_direction):
    for qz in range(0,nQuadPoints_per_element):
        for ey in range(0,nElements_per_direction):
            for qy in range(0,nQuadPoints_per_element):
                for ex in range(0,nElements_per_direction):
                    for qx in range(0,nQuadPoints_per_element):
                        wstr = "%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % \
                                (stored_data[ez,ey,ex,qz,qy,qx,0,0],stored_data[ez,ey,ex,qz,qy,qx,0,1],\
                                    stored_data[ez,ey,ex,qz,qy,qx,0,2],stored_data[ez,ey,ex,qz,qy,qx,0,3],\
                                    stored_data[ez,ey,ex,qz,qy,qx,0,4],stored_data[ez,ey,ex,qz,qy,qx,0,5])
                        # wstr = "%18.16e %18.16e %18.16e \n" % \
                        #         (stored_data[ez,ey,ex,qz,qy,qx,0,0],stored_data[ez,ey,ex,qz,qy,qx,0,1],\
                        #             stored_data[ez,ey,ex,qz,qy,qx,0,2]) # coords only
                        file.write(wstr)
file.close()
# NOTE: check that it works by doing 'diff reverse_order_test.dat setup_coords_only.dat'
exit()
# ================================================================
# check that it works
# ================================================================
setup_data = np.loadtxt("setup_coords_only-2.dat",dtype=np.float64)
readin_data_from_philip = np.loadtxt(reordered_velocity_file,dtype=np.float64)
file = open("check_reverse_reordering_of_philip_files_vs_setup-2.dat","w")
for i in range(0,nDOF):
    check = readin_data_from_philip[i,:]
    ref = setup_data[i,:]
    err = np.linalg.norm(check-ref)
    if(err > 2.0e-15):
        err_msg = "%i %18.16e \n" % (i,err)
        file.write(err_msg)
file.close()
# ================================================================