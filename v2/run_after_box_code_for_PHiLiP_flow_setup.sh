# TO DO: 
# Generate the box.in script automatically with the variables below (optional)

# After running the box.for code to generate the velocity field by running:
# # echo -e "box.in\n 1" | ./main.exe # generate the expected coordinates input
# # echo -e "box.in\n 0" | ./main.exe # run code as normal
# two files are generated:
# - (1) velocity_equidistant_nodes.fld
# - (2) velocity_gl_nodes.fld
# - (3) expected_input.fld

#====================================
# INPUTS
#====================================
# - number of elements per direction (must match value in box.in)
ELEMENTS_PER_DIRECTION=4
# - polynomial degree (must match value in box.in)
POLY_DEGREE=5
# - number of processors to run PHiLiP with
NPROCS=8

# STEP 1: Get the velocity field at unique points only by
#         averaging the values at the intersecting nodes
# # echo -e "velocity_equidistant_nodes.fld" | python3 prep_vel_for_spec.py 
# This will generate the following file:
# - (1) velocity_equidistant_nodes_unique.fld

# STEP 2: Determine the turbulent quantities
# Mean rms velocity, taylor microscale length


#====================================
# TESTING THE DHIT CASE
