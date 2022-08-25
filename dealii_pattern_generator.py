import numpy as np

data_dir = "philip_outputs/1procs/"
nElements_per_direction = 4
poly_degree = 5
nQuadPoints_per_element = poly_degree+1
nLines_per_pattern = nQuadPoints_per_element*nQuadPoints_per_element*nQuadPoints_per_element
nDOF = (nElements_per_direction*nQuadPoints_per_element)*(nElements_per_direction*nQuadPoints_per_element)*(nElements_per_direction*nQuadPoints_per_element)
nPatternLines = int(nDOF/nLines_per_pattern)
# print(nDOF)
# print(nLines_per_pattern)
# print(nPatternLines)
# exit()

data = np.loadtxt(data_dir + "coord_check-proc_0.txt",skiprows=1,dtype=np.float64)
# bounds = np.linspace(0.0,2.0*np.pi,nElements_per_direction)
def get_element_index(p1,p2):
    bounds = np.array([0.0000000000000000e+00,\
            1.5707963267948966e+00,\
            3.1415926535897931e+00,\
            4.7123889803846897e+00,\
            6.2831853071795862e+00])
    tol = 1.0e-10
    for el in range(0,nElements_per_direction):
        # if(np.abs(p1-bounds[e])<=tol and np.abs(p2-bounds[e+1])<=tol):
        # if(p1<=bounds[e] and p2<=bounds[e+1]):
        if(p2<=bounds[el+1]):
            return el
    return -1

file = open(data_dir + "pattern_test.dat","w")
# # file.write('Number of degrees of freedom:\n')
# wstr = "%i\n" % nDOF
# file.write(wstr)
ibase = 0
for p in range(0,nPatternLines):
    point1 = np.zeros(3)
    point2 = np.zeros(3)
    
    # point1 for determining pattern
    for d in range(0,3):
        point1[d] = data[ibase,d]
    
    # point2 for determining pattern
    point2[0] = data[ibase+1,0]
    point2[1] = data[ibase+nQuadPoints_per_element,1]
    point2[2] = data[ibase+nQuadPoints_per_element*nQuadPoints_per_element,2]
    
    if(p==22):
        print(point1)
        print(point2)

    for d in range(0,3):
        wstr = "%i " % get_element_index(point1[d],point2[d])
        file.write(wstr)
    wstr = "\n"
    file.write(wstr)
    ibase += nLines_per_pattern

# i += 1
file.close()