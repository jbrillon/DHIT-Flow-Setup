import numpy as np
#*****************************************************
#               GENERAL VARIABLES
#*****************************************************
nValues_per_row = 6 # 3 coords, 3 vel components
nElements_per_direction = 4
poly_degree = 5
nElements = nElements_per_direction*nElements_per_direction*nElements_per_direction
nQuadPoints_per_element = poly_degree + 1
nQuadPoints = nQuadPoints_per_element*nElements_per_direction
nDOF = (nElements_per_direction*nQuadPoints_per_element)*\
        (nElements_per_direction*nQuadPoints_per_element)*\
        (nElements_per_direction*nQuadPoints_per_element)
reduced_nDOF = (nElements_per_direction*nQuadPoints_per_element - (nElements_per_direction-1))*\
                (nElements_per_direction*nQuadPoints_per_element - (nElements_per_direction-1))*\
                (nElements_per_direction*nQuadPoints_per_element - (nElements_per_direction-1))

nLoops = 3
loop_bounds = np.ones(nLoops,dtype=np.int64)
if(nElements_per_direction>=4):
    loop_bounds[0] = 2
if(nElements_per_direction>=8):
    loop_bounds[1] = 2
if(nElements_per_direction>=16):
    loop_bounds[2] = 2
# if(nElements_per_direction>=32):
#     loop_bounds[3] = 2
# if(nElements_per_direction>=64):
#     loop_bounds[4] = 2
# if(nElements_per_direction>=128):
#     loop_bounds[5] = 2
# if(nElements_per_direction>=256):
#     loop_bounds[6] = 2
# if(nElements_per_direction>=512):
#     loop_bounds[7] = 2