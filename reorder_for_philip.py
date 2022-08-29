import numpy as np

# x = np.loadtxt("check.txt",usecols=0)
# npoints = np.size(x)

# for i in range(1,npoints):
#   print(x[i]-x[i-1])

nDOF = np.loadtxt("setup.dat",max_rows=1,dtype='int')
# nDOF = 13824 # for python 2.7
raw_data = np.loadtxt("setup.dat",skiprows=1,dtype=np.float64)
np.savetxt("reference_data.dat",raw_data)
nValues_per_row = 6
reordered_data = np.zeros((nDOF,nValues_per_row))

nElements_per_direction = 4
nElements = nElements_per_direction*nElements_per_direction*nElements_per_direction
poly_degree = 5
nQuadPoints_per_element = poly_degree + 1
nQuadPoints = nQuadPoints_per_element*nElements_per_direction

stored_data = np.zeros((nElements_per_direction,nElements_per_direction,nElements_per_direction,nQuadPoints_per_element,nQuadPoints_per_element,nQuadPoints_per_element,1,nValues_per_row))

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
                        # wstr = "%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % (row_data[0],row_data[1],row_data[2],row_data[3],row_data[4],row_data[5])
                        file.write(wstr)
                        i += 1
file.close()

#===========================================================
#                 REORDER DATA FOR PHiLiP
#===========================================================

file = open("reordered_data.dat","w")

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


data_dir = "philip_outputs/1procs/"
filename="coord_check_%i_elements_p%i-proc_0.txt" % (nElements_per_direction,poly_degree)
philip_data = np.loadtxt(data_dir + filename,skiprows=1,dtype=np.float64)

reordered_data = np.loadtxt("reordered_data.dat",dtype=np.float64)

file = open("check_reordering_vs_philip_output.dat","w")
for i in range(0,nDOF):
    check = reordered_data[i,:]
    ref = philip_data[i,:]
    err = np.linalg.norm(check-ref)
    if(err > 1.0e-12):
        err_msg = "%i %18.16e \n" % (i,err)
        file.write(err_msg)

file.close()

# file = open("reordered_data.dat","w")
# # file.write('Number of degrees of freedom:\n')
# wstr = "%i\n" % nDOF
# file.write(wstr)
# i = 0
# for el in range(0,nElements):
#     index_start = el*nQuadPoints_per_element
#     index_end = (el+1)*nQuadPoints_per_element
#     for qz in range(index_start,nQuadPoints):
    

#     for qy in range(0,nQuadPoints):
#         for qx in range(0,nQuadPoints):

#             wstr = "%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % \
#                     (stored_data[qz,qy,qx,0,0],stored_data[qz,qy,qx,0,1],\
#                         stored_data[qz,qy,qx,0,2],stored_data[qz,qy,qx,0,3],\
#                         stored_data[qz,qy,qx,0,4],stored_data[qz,qy,qx,0,5])
#             row_data = raw_data[i,:]
#             wstr = "%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % (row_data[0],row_data[1],row_data[2],row_data[3],row_data[4],row_data[5])
#             file.write(wstr)
#             i += 1
# file.close()



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
        
