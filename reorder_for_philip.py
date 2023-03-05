import numpy as np
from var import *

# get padded mpi rank string -- TO DO: Put this function somewhere since its being copied to a few files
def get_padded_mpi_rank_string(mpi_rank):
    padding_length = 5
    mpi_rank_string = '%i' % mpi_rank
    mpi_rank_string = mpi_rank_string.zfill(padding_length)
    return mpi_rank_string

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

def compute_mach_number(primitive_soln):
    global gamma_gas_minus_one
    density = 1.0*primitive_soln[0]
    velocity = 1.0*primitive_soln[1:4]
    pressure = 1.0*primitive_soln[4]

    speed_of_sound = np.sqrt(gamma_gas*pressure/density)
    velocity_magnitude = np.sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2])
    mach_number = velocity_magnitude/speed_of_sound
    return mach_number

# Global constants
global gamma_gas, gamma_gas_minus_one
gamma_gas = 1.4
gamma_gas_minus_one = gamma_gas - 1.0

maximum_desired_turbulent_mach_number = 0.3

print(" ")
print("----------------------------------------------------------------------------")
print("Nondimensionalized quantities for initializing the flow:")
print("----------------------------------------------------------------------------")
temperature, prandtl_number, reynolds_number_ref, mach_number_ref, nondim_density, nondim_mean_velocity, nondim_pressure = \
    np.loadtxt("parameters_for_dhit_setup.txt",dtype=np.float64,unpack=True)
nondimensionalized_pressure = 1.0*nondim_pressure
nondimensionalized_density = 1.0*nondim_density
nondimensionalized_mean_velocity = 1.0*nondim_mean_velocity

print(" - Nondimensionalized density: %f" % nondimensionalized_density)
print(" - Nondimensionalized mean velocity: %f" % nondimensionalized_mean_velocity)
print(" - Nondimensionalized pressure: %f" % nondimensionalized_pressure)

# (1) Pre-process file to generate the setup.dat file -- same way the deprecated MATLAB script did
subdir = "dofs128_p3_velocity"
input_vel_field_file = subdir+"/velocity_gl_nodes.fld"
input_data = np.loadtxt(input_vel_field_file,dtype=np.float64)
nDOF_input = input_data.shape[0]
file = open("setup.dat","w")
# file.write('Number of degrees of freedom:\n')
wstr = "%i\n" % nDOF_input
file.write(wstr)
i = 0
for i in range(0,nDOF_input):
    wstr = "%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % \
            (input_data[i,0],input_data[i,1],\
                input_data[i,2],input_data[i,3],\
                input_data[i,4],input_data[i,5])
    file.write(wstr)
file.close()

# (2) Read in the data from setup.dat
nDOF_expected = np.loadtxt("setup.dat",max_rows=1,dtype='int')
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
                        # - Nondimensionalized velocity components; add nondimensionalized mean velocity to x-component
                        nondimensionalized_primitive_sol_at_q_point[1] = nondimensionalized_mean_velocity + stored_data[ez,ey,ex,qz,qy,qx,0,3]
                        nondimensionalized_primitive_sol_at_q_point[2] = stored_data[ez,ey,ex,qz,qy,qx,0,4]
                        nondimensionalized_primitive_sol_at_q_point[3] = stored_data[ez,ey,ex,qz,qy,qx,0,5]

                        # Check turbulent mach number -- may want to remove this
                        mach_number_at_q_point = compute_mach_number(nondimensionalized_primitive_sol_at_q_point)
                        if(mach_number_at_q_point > maximum_desired_turbulent_mach_number):
                            print("Error: a local mach_number exceeds the maximum_desired_turbulent_mach_number")
                            print("Value: %.3e" % mach_number_at_q_point)
                            print("Aborting...")
                            exit()


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

nDOF_per_proc = nDOF/num_procs
philip_prefix=subdir+"/setup_files/setup"

# TO DO:
# if(nDOF % nDOFs_per_proc != 0):
#     print("ERROR: Must use a number of processors that evenly divides the ")

file = open("reordered_data.dat","w")
wstr = "%i\n" % nDOF
file.write(wstr)

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
                                                                                padded_mpi_rank_string = get_padded_mpi_rank_string(iproc)
                                                                                filename_for_philip="%s-%s.dat" % (philip_prefix,padded_mpi_rank_string)
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
            ex_L_base_base_base = ex_L_base_base
        ey_L_base_base_base = ey_L_base_base
    ez_L_base_base_base = ez_L_base_base

file.close()
exit()
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
