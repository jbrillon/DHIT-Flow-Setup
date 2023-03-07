import numpy as np
import sys; sys.path.append("../submodules/TurboGenPY"); import TurboGenPY_HighOrderFEM as turboFEM
import convert_equidistant_to_gauss_lobatto_nodes as eq2gll
# import sys; sys.path.append("../submodules/quickplotlib/lib"); import quickplotlib as qp
#-----------------------------------------------------
def get_DOF_vars(nElements_per_direction,poly_degree):
    nElements = nElements_per_direction*nElements_per_direction*nElements_per_direction
    nQuadPoints_per_element = poly_degree + 1
    nQuadPoints = nQuadPoints_per_element*nElements_per_direction
    nDOF = nQuadPoints*nQuadPoints*nQuadPoints
    reduced_nQuadPoints = nElements_per_direction*nQuadPoints_per_element - (nElements_per_direction-1)
    reduced_nDOF = reduced_nQuadPoints*reduced_nQuadPoints*reduced_nQuadPoints
    return [nElements,nQuadPoints_per_element,nQuadPoints,nDOF,reduced_nQuadPoints,reduced_nDOF]
#-----------------------------------------------------
def get_reduced_nDOF_and_nQuadPoints(nElements_per_direction,poly_degree):
    reduced_nQuadPoints = get_DOF_vars(nElements_per_direction,poly_degree)[-2]
    reduced_nDOF = get_DOF_vars(nElements_per_direction,poly_degree)[-1]
    return [reduced_nQuadPoints,reduced_nDOF]
#-----------------------------------------------------
#*****************************************************
#               GENERAL VARIABLES
#*****************************************************
#=====================================================
# INPUTS
#=====================================================
nValues_per_row = 6 # 3 coords, 3 vel components
nElements_per_direction = 4 # number of elements in each direction/dimension
poly_degree = 5 # solution polynomial degree
num_procs = 64 # 1024 # number of mpi processors to run PHiLiP with
output_directory = "./"
# Get DOF variables
nElements,nQuadPoints_per_element,nQuadPoints,nDOF,reduced_nQuadPoints,reduced_nDOF = get_DOF_vars(nElements_per_direction,poly_degree)
#-----------------------------------------------------
# Fixed variables
#-----------------------------------------------------
nLoops = 8
loop_bounds = np.ones(nLoops,dtype=np.int64)
if(nElements_per_direction>=4):
    loop_bounds[0] = 2
if(nElements_per_direction>=8):
    loop_bounds[1] = 2
if(nElements_per_direction>=16):
    loop_bounds[2] = 2
if(nElements_per_direction>=32):
    loop_bounds[3] = 2
if(nElements_per_direction>=64):
    loop_bounds[4] = 2
if(nElements_per_direction>=128):
    loop_bounds[5] = 2
if(nElements_per_direction>=256):
    loop_bounds[6] = 2
if(nElements_per_direction>=512):
    loop_bounds[7] = 2
#-----------------------------------------------------
#=====================================================
# Generate iso turb files for PHiLiP code
#=====================================================

import os; print(os.getcwd())
# generate the velocity flow field
turboFEM.generate_isotropic_turbulence_high_order_fem(
    nElements_per_direction,
    poly_degree,
    output_filename=(output_directory+"velocity_equidistant_nodes.fld"),
    number_of_modes=5000, # suggested in TurboGenPY paper
    spectra_name='ml'
    )

# convert to GLL nodes
eq2gll.convert_equidistant_to_gauss_lobatto_nodes(
    (output_directory+"velocity_equidistant_nodes.fld"),
    nElements_per_direction,
    nQuadPoints_per_element,
    nValues_per_row,
    poly_degree,
    nDOF,
    output_filename=(output_directory+"velocity_gl_nodes.fld"),
    test_reading=False # set a true if testing that this function works
    )

# run the ML scaling calc


