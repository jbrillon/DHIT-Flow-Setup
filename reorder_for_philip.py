import numpy as np
from var import *

def get_conservative_from_primitive(primitive_soln):
    global gamma_gas_minus_one
    density = 1.0*primitive_soln[0]
    velocity = 1.0*primitive_soln[1:4]
    pressure = 1.0*primitive_soln[4]

    conservative_soln = np.zeros(5,dtype=np.float64)
    conservative_soln[0] = 1.0*density
    conservative_soln[1:4] = density*velocity
    conservative_soln[4] = pressure/(gamma_gas_minus_one) + 0.5*density*np.dot(velocity,velocity)
    return conservative_soln

# Constants
characteristic_length = 1.0 # L

# Global constants
global gamma_gas, gamma_gas_minus_one
gamma_gas = 1.4
gamma_gas_minus_one = 0.4

# TO DO: Read in setup_more.dat and assign the following variables
freestream_velocity = 1.0
freestream_mach_number = 1.0
eddy_turnover_time = 1.0

# TO DO: discuss with Brian the proper nondimensionalization
#        since the Navier-Stokes equations themselves have to
#        be initialized with nondimensionalized values based
#        on freestream quantities

# Determine final time for PHiLiP
nondimensional_eddy_turnover_time = eddy_turnover_time/(characteristic_length/freestream_velocity)
nondimensionalized_density = 1.0
nondimensionalized_pressure = 1.0

# nDOF_expected = np.loadtxt("setup.dat",max_rows=1,dtype='int')
nDOF_expected = 13824 # for python 2.7

if(nDOF!=nDOF_expected):
    print("Error: nDOF does not match expected nDOF from file, check var.py")
    print("Aborting...")
    exit()

raw_data = np.loadtxt("setup.dat",skiprows=1,dtype=np.float64)

stored_data = np.zeros((nElements_per_direction,nElements_per_direction,nElements_per_direction,nQuadPoints_per_element,nQuadPoints_per_element,nQuadPoints_per_element,1,nValues_per_row),dtype=np.float64)
nondimensionalized_conservative_solution = np.zeros((nElements_per_direction,nElements_per_direction,nElements_per_direction,nQuadPoints_per_element,nQuadPoints_per_element,nQuadPoints_per_element,1,5),dtype=np.float64)

file = open("read_test.dat","w")
# file.write('Number of degrees of freedom:\n')
wstr = "%i\n" % nDOF
file.write(wstr)
i = 0
for ez in range(0,nElements_per_direction):
    for qz in range(0,nQuadPoints_per_element):
        for ey in range(0,nElements_per_direction):
            for qy in range(0,nQuadPoints_per_element):
                for ex in range(0,nElements_per_direction):
                    for qx in range(0,nQuadPoints_per_element):
                        row_data = raw_data[i,:]
                        for iValue in range(0,nValues_per_row):
                            stored_data[ez,ey,ex,qz,qy,qx,0,iValue] = row_data[iValue]

                        wstr = "%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % \
                                (stored_data[ez,ey,ex,qz,qy,qx,0,0],stored_data[ez,ey,ex,qz,qy,qx,0,1],\
                                    stored_data[ez,ey,ex,qz,qy,qx,0,2],stored_data[ez,ey,ex,qz,qy,qx,0,3],\
                                    stored_data[ez,ey,ex,qz,qy,qx,0,4],stored_data[ez,ey,ex,qz,qy,qx,0,5])
                        file.write(wstr)
                        i += 1

                        # ------------------------------ START -------------------------------
                        # ON THE FLY PRE-PROCESSING
                        # --------------------------------------------------------------------
                        # Primitive solution at current point
                        nondimensionalized_primitive_sol_at_q_point = np.zeros(5,dtype=np.float64)
                        nondimensionalized_primitive_sol_at_q_point[0] = nondimensionalized_density
                        nondimensionalized_primitive_sol_at_q_point[4] = nondimensionalized_pressure
                        # -- Apply freestream non-dimensionalization to velocity components
                        nondimensionalized_primitive_sol_at_q_point[1] = stored_data[ez,ey,ex,qz,qy,qx,0,3]/freestream_velocity
                        nondimensionalized_primitive_sol_at_q_point[2] = stored_data[ez,ey,ex,qz,qy,qx,0,4]/freestream_velocity
                        nondimensionalized_primitive_sol_at_q_point[3] = stored_data[ez,ey,ex,qz,qy,qx,0,5]/freestream_velocity

                        # Conservative solution at current point
                        nondimensionalized_conservative_sol_at_q_point = 1.0*get_conservative_from_primitive(nondimensionalized_primitive_sol_at_q_point)
                        # Store conservative solution
                        for state in range(0,5):
                            nondimensionalized_conservative_solution[ez,ey,ex,qz,qy,qx,0,state] = 1.0*nondimensionalized_conservative_sol_at_q_point[state]
                        # ------------------------------- END --------------------------------
file.close()
# NOTE: Do 'diff read_test.dat setup.dat' to make sure we're reading this properly

#===========================================================
#                 REORDER DATA FOR PHiLiP
#===========================================================

num_procs = 4
nDOF_per_proc = nDOF/num_procs
philip_prefix="setup_philip"

# TO DO:
# if(nDOF % nDOFs_per_proc != 0):
#     print("ERROR: Must use a number of processors that evenly divides the ")

file = open("reordered_data.dat","w")
wstr = "%i\n" % nDOF
file.write(wstr)

''' must add more nested for loops for higher
    number of elements per direction
    currently can handle up to 16 (i.e. 2,4,8,16)
'''

iproc = 0
iDOF_per_proc = 0
start_new_file=True

ez_L_base_base = 0
for z_base_base in range(0,loop_bounds[2]):
    ey_L_base_base = 0
    for y_base_base in range(0,loop_bounds[2]):
        ex_L_base_base = 0
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
                                                            # wstr = "%i %i %i \n" % (ex,ey,ez)
                                                            # wstr = "%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % \
                                                            #         (stored_data[ez,ey,ex,qz,qy,qx,0,0],stored_data[ez,ey,ex,qz,qy,qx,0,1],\
                                                            #             stored_data[ez,ey,ex,qz,qy,qx,0,2],stored_data[ez,ey,ex,qz,qy,qx,0,3],\
                                                            #             stored_data[ez,ey,ex,qz,qy,qx,0,4],stored_data[ez,ey,ex,qz,qy,qx,0,5])
                                                            wstr = "%18.16e %18.16e %18.16e \n" % \
                                                                    (stored_data[ez,ey,ex,qz,qy,qx,0,0],stored_data[ez,ey,ex,qz,qy,qx,0,1],\
                                                                        stored_data[ez,ey,ex,qz,qy,qx,0,2])
                                                            file.write(wstr)

                                                for state in range(0,5):
                                                    for qz in range(0,nQuadPoints_per_element):
                                                        for qy in range(0,nQuadPoints_per_element):
                                                            for qx in range(0,nQuadPoints_per_element):
                                                                if(state==0 and start_new_file):
                                                                    filename_for_philip="%s-0000%i.dat" % (philip_prefix,iproc)
                                                                    file_for_philip = open(filename_for_philip,"w")
                                                                    wstr = "%i\n" % nDOF
                                                                    file_for_philip.write(wstr)
                                                                    start_new_file=False

                                                                wstr2 = "%18.16e %18.16e %18.16e %i %18.16e\n" % \
                                                                        (stored_data[ez,ey,ex,qz,qy,qx,0,0],\
                                                                            stored_data[ez,ey,ex,qz,qy,qx,0,1],\
                                                                            stored_data[ez,ey,ex,qz,qy,qx,0,2],\
                                                                            state,\
                                                                            nondimensionalized_conservative_solution[ez,ey,ex,qz,qy,qx,0,state])
                                                                file_for_philip.write(wstr2)

                                                                if(state==4):
                                                                    iDOF_per_proc += 1

                                                                if(state==4 and (iDOF_per_proc == nDOF_per_proc)):
                                                                    file_for_philip.close()
                                                                    iDOF_per_proc = 0
                                                                    iproc += 1
                                                                    start_new_file=True

                                    ex_L += 2
                                ey_L += 2
                            ez_L += 2
                        ex_L_base = ex_L
                    ey_L_base = ey_L
                ez_L_base = ez_L
            ex_L_base_base = ex_L_base
        ey_L_base_base = ey_L_base
    ez_L_base_base = ez_L_base

file.close()

# ================================================================
# check that it works
# ================================================================
data_dir = "philip_outputs/"
# filename="1procs/coord_check_%i_elements_p%i-proc_0.txt" % (nElements_per_direction,poly_degree)
filename="8procs/assembled_coords.txt"
philip_data = np.loadtxt(data_dir + filename,skiprows=1,dtype=np.float64)
reordered_data = np.loadtxt("reordered_data.dat",skiprows=1,dtype=np.float64)
file = open("check_reordering_vs_philip_output.dat","w")
for i in range(0,nDOF):
    check = reordered_data[i,:]
    ref = philip_data[i,:]
    err = np.linalg.norm(check-ref)
    if(err > 2.0e-15):
        err_msg = "%i %18.16e \n" % (i,err)
        file.write(err_msg)
file.close()
# ================================================================
