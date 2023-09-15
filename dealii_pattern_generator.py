import numpy as np

data_dir = "philip_outputs/1procs/"

global nElements_per_direction
nElements_per_direction = 4
poly_degree = 5
nQuadPoints_per_element = poly_degree+1
nLines_per_pattern = nQuadPoints_per_element*nQuadPoints_per_element*nQuadPoints_per_element
nElements_total = nElements_per_direction*nElements_per_direction*nElements_per_direction
nDOF = (nElements_per_direction*nQuadPoints_per_element)*(nElements_per_direction*nQuadPoints_per_element)*(nElements_per_direction*nQuadPoints_per_element)
nPatternLines = int(nDOF/nLines_per_pattern)
# print(nDOF)
# print(nLines_per_pattern)
# print(nPatternLines)
# exit()

filename="coord_check_%i_elements_p%i-proc_0.txt" % (nElements_per_direction,poly_degree)
data = np.loadtxt(data_dir + filename,skiprows=1,dtype=np.float64)

global bounds
bounds = np.linspace(0.0,2.0*np.pi,nElements_per_direction+1)

def get_element_index(p1,p2):
    global bounds,nElements_per_direction
    tol = 1.0e-6
    for el in range(0,nElements_per_direction):
        # if(np.abs(p1-bounds[e])<=tol and np.abs(p2-bounds[e+1])<=tol):
        # if(p1<=bounds[e] and p2<=bounds[e+1]):
        if(p2<bounds[el+1]):
            return el
    return -1

file = open(data_dir + "pattern_gen.dat","w")
ibase = 0
for p in range(0,nPatternLines):
    point1 = np.zeros(3,dtype=np.float64)
    point2 = np.zeros(3,dtype=np.float64)
    
    # point1 for determining pattern
    for d in range(0,3):
        point1[d] = data[ibase,d]
    
    # point2 for determining pattern
    point2[0] = data[ibase+1,0]
    point2[1] = data[ibase+nQuadPoints_per_element,1]
    point2[2] = data[ibase+nQuadPoints_per_element*nQuadPoints_per_element,2]

    for d in range(0,3):
        wstr = "%i " % get_element_index(point1[d],point2[d])
        file.write(wstr)
    wstr = "\n"
    file.write(wstr)
    ibase += nLines_per_pattern
file.close()

# Check pattern

file = open(data_dir + "pattern_check.dat","w")

nLoops = 3
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
# if(nElements_per_direction>=128):
#     loop_bounds[5] = 2
# if(nElements_per_direction>=256):
#     loop_bounds[6] = 2
# if(nElements_per_direction>=512):
#     loop_bounds[7] = 2

''' must add more nested for loops for higher
    number of elements per direction
    currently can handle up to 16 (i.e. 2,4,8,16)
'''

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
                        # algorithm for a cube with 64 (4^3) elements
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
                                                wstr = "%i %i %i \n" % (ex,ey,ez)
                                                file.write(wstr)
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