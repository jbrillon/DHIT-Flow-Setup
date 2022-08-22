import numpy as np

# x = np.loadtxt("check.txt",usecols=0)
# npoints = np.size(x)

# for i in range(1,npoints):
#   print(x[i]-x[i-1])

nDOF = np.loadtxt("setup.dat",max_rows=1,dtype='int')
raw_data = np.loadtxt("setup.dat",skiprows=1,dtype=np.float64)

reordered_data = np.zeros((nDOF,6))

nElements = 4
poly_degree = 5
nQuadPoints = poly_degree + 1

for q2 in range(0,nQuadPoints):
    print("==================================================")
    
    for q1 in range(0,nQuadPoints):
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
        for el in range(0,nElements):
            # index of raw data
            row_index_start = q2*(nQuadPoints*(nElements*nQuadPoints)) + q1*(nElements*nQuadPoints) + el*nQuadPoints
            row_index_end = row_index_start + nQuadPoints
            raw_values = raw_data[row_index_start:(row_index_end),:]

            print(raw_values[:,0:3])
            print("--------------------------------------------")

            # index of reordered data
            # row_index_start_reordered_data = el*(3*nQuadPoints) + d*nQuadPoints
            # row_index_end_reordered_data = row_index_start_reordered_data + nQuadPoints

            # reordered_data[row_index_start_reordered_data:(row_index_end_reordered_data),:] = raw_values

# np.savetxt("reordered_data.dat", reordered_data, fmt='%.18e', delimiter=' ')
        
