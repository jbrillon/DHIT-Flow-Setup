import numpy as np

# x = np.loadtxt("check.txt",usecols=0)
# npoints = np.size(x)

# for i in range(1,npoints):
#   print(x[i]-x[i-1])

# nDOF = np.loadtxt("setup.dat",max_rows=1,dtype='int')
nDOF = 13824
raw_data = np.loadtxt("setup.dat",skiprows=1,dtype=np.float64)
np.savetxt("reference_data.dat",raw_data)
nValues_per_row = 6
reordered_data = np.zeros((nDOF,nValues_per_row))

nElements = 4
poly_degree = 5
nQuadPoints_per_element = poly_degree + 1
nQuadPoints = nQuadPoints_per_element*nElements

stored_data = np.zeros((nQuadPoints,nQuadPoints,nQuadPoints,1,nValues_per_row))

file = open("read_test.dat","w")
# file.write('Number of degrees of freedom:\n')
wstr = "%i\n" % nDOF
file.write(wstr)
i = 0
for qz in range(0,nQuadPoints):
    for qy in range(0,nQuadPoints):
        for qx in range(0,nQuadPoints):
            row_data = raw_data[i,:]
            for iValue in range(0,nValues_per_row):
                stored_data[qz,qy,qx,0,iValue] = row_data[iValue]

            wstr = "%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % \
                    (stored_data[qz,qy,qx,0,0],stored_data[qz,qy,qx,0,1],\
                        stored_data[qz,qy,qx,0,2],stored_data[qz,qy,qx,0,3],\
                        stored_data[qz,qy,qx,0,4],stored_data[qz,qy,qx,0,5])
            # wstr = "%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % (row_data[0],row_data[1],row_data[2],row_data[3],row_data[4],row_data[5])
            file.write(wstr)
            i += 1
file.close()


file = open("reordered_data.dat","w")
# file.write('Number of degrees of freedom:\n')
wstr = "%i\n" % nDOF
file.write(wstr)
i = 0
for el in range(0,nElements):
    index_start = el*nQuadPoints_per_element
    index_end = (el+1)*nQuadPoints_per_element
    for qz in range(index_start,nQuadPoints):
    

    for qy in range(0,nQuadPoints):
        for qx in range(0,nQuadPoints):

            wstr = "%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % \
                    (stored_data[qz,qy,qx,0,0],stored_data[qz,qy,qx,0,1],\
                        stored_data[qz,qy,qx,0,2],stored_data[qz,qy,qx,0,3],\
                        stored_data[qz,qy,qx,0,4],stored_data[qz,qy,qx,0,5])
            row_data = raw_data[i,:]
            wstr = "%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % (row_data[0],row_data[1],row_data[2],row_data[3],row_data[4],row_data[5])
            file.write(wstr)
            i += 1
file.close()



# nQuadPoints = poly_degree + 1
# for q2 in range(0,nQuadPoints*nElements):
#     print("==================================================")
    
#     for q1 in range(0,nQuadPoints):
#         print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
#         for el in range(0,nElements):
#             # index of raw data
#             row_index_start = q2*(nQuadPoints*(nElements*nQuadPoints)) + q1*(nElements*nQuadPoints) + el*nQuadPoints
#             row_index_end = row_index_start + nQuadPoints
#             raw_values = raw_data[row_index_start:(row_index_end),:]

#             print(raw_values[:,0:3])
#             print("--------------------------------------------")

#             # index of reordered data
#             # row_index_start_reordered_data = el*(3*nQuadPoints) + d*nQuadPoints
#             # row_index_end_reordered_data = row_index_start_reordered_data + nQuadPoints

#             # reordered_data[row_index_start_reordered_data:(row_index_end_reordered_data),:] = raw_values

# # np.savetxt("reordered_data.dat", reordered_data, fmt='%.18e', delimiter=' ')
        
