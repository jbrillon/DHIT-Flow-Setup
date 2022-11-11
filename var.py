import numpy as np
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
nValues_per_row = 6 # 3 coords, 3 vel components
nElements_per_direction = 4
poly_degree = 5
# Get DOF variables
nElements,nQuadPoints_per_element,nQuadPoints,nDOF,reduced_nQuadPoints,reduced_nDOF = get_DOF_vars(nElements_per_direction,poly_degree)
#-----------------------------------------------------
nLoops = 4
loop_bounds = np.ones(nLoops,dtype=np.int64)
if(nElements_per_direction>=4):
    loop_bounds[0] = 2
if(nElements_per_direction>=8):
    loop_bounds[1] = 2
if(nElements_per_direction>=16):
    loop_bounds[2] = 2
if(nElements_per_direction>=32):
    loop_bounds[3] = 2
# if(nElements_per_direction>=64):
#     loop_bounds[4] = 2
# if(nElements_per_direction>=128):
#     loop_bounds[5] = 2
# if(nElements_per_direction>=256):
#     loop_bounds[6] = 2
# if(nElements_per_direction>=512):
#     loop_bounds[7] = 2
#-----------------------------------------------------