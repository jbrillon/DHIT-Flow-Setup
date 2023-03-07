import numpy as np

def convert_equidistant_to_gauss_lobatto_nodes(
    equidistant_data_filename,
    nElements_per_direction,
    nQuadPoints_per_element,
    nValues_per_row,
    poly_degree,
    nDOF,
    output_filename="velocity_gl_nodes.fld",
    test_reading=False
    ):

    # future work: print progress in terms of iDOF/nDOF ever x steps or print the percentage complete
    # this will be slow when doing high DOFs

    nDOF_expected = np.loadtxt(equidistant_data_filename,skiprows=0,max_rows=1,dtype=np.float64)#raw_data.shape[0]
    if(nDOF!=nDOF_expected):
        print("Error: nDOF does not match expected nDOF from file %s, check var.py",equidistant_data_filename)
        print("Aborting...")
        exit()

    # load data at equidistant nodes
    raw_data = np.loadtxt(equidistant_data_filename,skiprows=1,dtype=np.float64)

    # Store the data
    stored_data = np.zeros((nElements_per_direction,nElements_per_direction,nElements_per_direction,nQuadPoints_per_element,nQuadPoints_per_element,nQuadPoints_per_element,1,nValues_per_row),dtype=np.float64)
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
                            i += 1

    # FROM WRITE FLD
    ng = 1*nElements_per_direction*nQuadPoints_per_element
    ne = 1*nElements_per_direction
    pp = nQuadPoints_per_element #poly_degree + 1
    Ni1 = 1
    Nj1 = 1
    Nk1 = 1
    Ni2 = 1*ng
    Nj2 = 1*ng
    Nk2 = 1*ng
    Ni3 = 1*ne
    Nj3 = 1*ne
    Nk3 = 1*ne
    # read it in the way it was written
    i_global = 0
    if(test_reading):
        Xg_in = np.zeros((ng,ng,ng),dtype=np.float64) # for testing
        Yg_in = np.zeros((ng,ng,ng),dtype=np.float64) # for testing
        Zg_in = np.zeros((ng,ng,ng),dtype=np.float64) # for testing
    Vel = np.zeros((ng,ng,ng,3),dtype=np.float64)
    Vel2 = np.zeros((ng,ng,ng,3),dtype=np.float64)
    for k in range(Nk1-1, Nk2):
        for j in range(Nj1-1, Nj2):
            for i in range(Ni1-1, Ni2):
                if(test_reading):
                    Xg_in[i,j,k] = raw_data[i_global,0]
                    Yg_in[i,j,k] = raw_data[i_global,1]
                    Zg_in[i,j,k] = raw_data[i_global,2]
                Vel[i,j,k,0] = raw_data[i_global,3]
                Vel[i,j,k,1] = raw_data[i_global,4]
                Vel[i,j,k,2] = raw_data[i_global,5]
                i_global += 1
                # do same for velocities if it works
    # NOTE: Do 'diff read_test.dat setup.dat' to make sure we're reading this properly

    #======================================================
    # Generate Grid
    #======================================================

    # CARTESIAN GRIDPOINTS
    coeff1 = np.zeros((5,6),dtype=np.float64)
    coeff1[1-1,1-1] = -1.0000000000000000e0
    coeff1[1-1,2-1] = +1.0000000000000000e0
    coeff1[2-1,1-1] = -1.0000000000000000e0
    coeff1[2-1,2-1] = +0.0000000000000000e0
    coeff1[2-1,3-1] = +1.0000000000000000e0
    coeff1[3-1,1-1] = -1.0000000000000000e0
    coeff1[3-1,2-1] = -4.47213595499957939e-1
    coeff1[3-1,3-1] = +4.47213595499957939e-1
    coeff1[3-1,4-1] = +1.0000000000000000e0
    coeff1[4-1,1-1] = -1.0000000000000000e0
    coeff1[4-1,2-1] = -6.54653670707977144e-1
    coeff1[4-1,3-1] = +0.0000000000000000e0
    coeff1[4-1,4-1] = +6.54653670707977144e-1
    coeff1[4-1,5-1] = +1.0000000000000000e0
    coeff1[5-1,1-1] = -1.0000000000000000e0
    coeff1[5-1,2-1] = -7.6505532392946469e-1
    coeff1[5-1,3-1] = -2.85231516480645096e-1
    coeff1[5-1,4-1] = +2.85231516480645096e-1
    coeff1[5-1,5-1] = +7.6505532392946469e-1
    coeff1[5-1,6-1] = +1.0000000000000000e0

    # OUTPUT GRIDPOINTS
    coeff2 = np.zeros((5,6),dtype=np.float64)
    coeff2[1-1,1-1] = -1.0000000000000000e0
    coeff2[1-1,2-1] = +1.0000000000000000e0
    coeff2[2-1,1-1] = -1.0000000000000000e0
    coeff2[2-1,2-1] = +0.0000000000000000e0
    coeff2[2-1,3-1] = +1.0000000000000000e0
    coeff2[3-1,1-1] = -1.0000000000000000e0
    coeff2[3-1,2-1] = -3.3333333333333333e-1
    coeff2[3-1,3-1] = +3.3333333333333333e-1
    coeff2[3-1,4-1] = +1.0000000000000000e0
    coeff2[4-1,1-1] = -1.0000000000000000e0
    coeff2[4-1,2-1] = -5.0000000000000000e-1
    coeff2[4-1,3-1] = +0.0000000000000000e0
    coeff2[4-1,4-1] = +5.0000000000000000e-1
    coeff2[4-1,5-1] = +1.0000000000000000e0
    coeff2[5-1,1-1] = -1.0000000000000000e0
    coeff2[5-1,2-1] = -6.0000000000000000e-1
    coeff2[5-1,3-1] = -2.0000000000000000e-1
    coeff2[5-1,4-1] = +2.0000000000000000e-1
    coeff2[5-1,5-1] = +6.0000000000000000e-1
    coeff2[5-1,6-1] = +1.0000000000000000e0

    # safeguard
    if(poly_degree>5):
        print("ERROR: Convert equidistant to GLL nodes not implemented for poly_degree>5.")
        print("Please add quadrature weights/nodes for higher poly_degree and update this file.")
        print("Aborting...")
        exit()

    ng = 1*nElements_per_direction*nQuadPoints_per_element
    ne = 1*nElements_per_direction
    pp = nQuadPoints_per_element
    Ni1 = 1
    Nj1 = 1
    Nk1 = 1
    Ni2 = 1*ng
    Nj2 = 1*ng
    Nk2 = 1*ng
    Ni3 = 1*ne
    Nj3 = 1*ne
    Nk3 = 1*ne

    Xg = np.zeros((ng,ng,ng),dtype=np.float64)
    Yg = np.zeros((ng,ng,ng),dtype=np.float64)
    Zg = np.zeros((ng,ng,ng),dtype=np.float64)
    Xs = np.zeros((ng,ng,ng),dtype=np.float64)
    Ys = np.zeros((ng,ng,ng),dtype=np.float64)
    Zs = np.zeros((ng,ng,ng),dtype=np.float64)

    dh = 2.0*np.pi
    dh = dh/np.float64(ne)

    if(test_reading):
        file = open("test_reading_equidistant_to_gl_nodes.txt","w") # for testing
        i_check = 0 # for testing

    p = 1*poly_degree-1 # minus one because of indexing
    for k in range(Nk1-1, Nk3):
        for j in range(Nj1-1, Nj3):
            for i in range(Ni1-1, Ni3):
                for l in range(1-1, pp):
                    for m in range(1-1, pp):
                        for n in range(1-1, pp):

                            ri = np.float64(i)
                            rj = np.float64(j)
                            rk = np.float64(k)
                            rl = np.float64(l)
                            rm = np.float64(m)
                            rn = np.float64(n)
                            rp = np.float64(p)

                            indi = (i)*pp+l
                            indj = (j)*pp+m
                            indk = (k)*pp+n

                            #At Cartesian Points
                            Xg[indi,indj,indk] = (ri)*dh + (coeff2[p,l]+1.0)/2.0*dh;
                            Yg[indi,indj,indk] = (rj)*dh + (coeff2[p,m]+1.0)/2.0*dh;
                            Zg[indi,indj,indk] = (rk)*dh + (coeff2[p,n]+1.0)/2.0*dh;

                            #At Solution Points
                            Xs[indi,indj,indk] = (ri)*dh + (coeff1[p,l]+1.0)/2.0*dh;
                            Ys[indi,indj,indk] = (rj)*dh + (coeff1[p,m]+1.0)/2.0*dh;
                            Zs[indi,indj,indk] = (rk)*dh + (coeff1[p,n]+1.0)/2.0*dh;

                            if(test_reading):
                                # for testing
                                check = np.array([Xg[indi,indj,indk],Yg[indi,indj,indk],Zg[indi,indj,indk]])
                                ref = np.array([Xg_in[indi,indj,indk],Yg_in[indi,indj,indk],Zg_in[indi,indj,indk]])
                                err = np.linalg.norm(check-ref)
                                # err_msg = "%7.5e %7.5e %7.5e | %7.5e %7.5e %7.5e \n" % (Xg[indi,indj,indk],Yg[indi,indj,indk],Zg[indi,indj,indk],raw_data[i_check,0],raw_data[i_check,1],raw_data[i_check,2])
                                # err_msg = "%7.5e %7.5e %7.5e | %7.5e %7.5e %7.5e \n" % (Xg[indi,indj,indk],Yg[indi,indj,indk],Zg[indi,indj,indk],Xg_in[indi,indj,indk],Yg_in[indi,indj,indk],Zg_in[indi,indj,indk])
                                # file.write(err_msg)
                                if(err > -1.0):
                                    err_msg = "%i %18.16e \n" % (i_check,err)
                                    file.write(err_msg)
                                i_check += 1
    if(test_reading):
        file.close()

    #======================================================
    # Generate velocity at GL nodes
    #======================================================
    # CONVERT FROM CARTESIAN GRID TO THE SOLUTION GRID
    ng = 1*nElements_per_direction*nQuadPoints_per_element
    ne = 1*nElements_per_direction
    pp = nQuadPoints_per_element #poly_degree + 1
    Ni1 = 1
    Nj1 = 1
    Nk1 = 1
    Ni2 = 1*ng
    Nj2 = 1*ng
    Nk2 = 1*ng
    Ni3 = 1*ne
    Nj3 = 1*ne
    Nk3 = 1*ne

    dh = 2.0*np.pi/np.float64(ne)

    for k in range(Nk1-1, Nk3):
        for j in range(Nj1-1, Nj3):
            for i in range(Ni1-1, Ni3):
                for l in range(1-1, pp):
                    for m in range(1-1, pp):
                        for n in range(1-1, pp):

                            indi = (i)*pp+l
                            indj = (j)*pp+m
                            indk = (k)*pp+n

                            #RESET INTERP VECTOR
                            Vel2[indi,indj,indk,1-1] = 0.0
                            Vel2[indi,indj,indk,2-1] = 0.0
                            Vel2[indi,indj,indk,3-1] = 0.0

                            # POSITION IN THE SOLUTION SPACE
                            X1 = Xs[indi,indj,indk]
                            Y1 = Ys[indi,indj,indk]
                            Z1 = Zs[indi,indj,indk]

                            # BUILD THE LAGRANGE INTERPOLATION
                            for ll in range(1-1, pp):
                                for mm in range(1-1, pp):
                                    for nn in range(1-1, pp):
                                        indi = (i)*pp+ll
                                        indj = (j)*pp+mm
                                        indk = (k)*pp+nn

                                        # POSITION IN THE CARTESIAN SPACE
                                        X2 = Xg[indi,indj,indk]
                                        Y2 = Yg[indi,indj,indk]
                                        Z2 = Zg[indi,indj,indk]

                                        U2 = Vel[indi,indj,indk,0]
                                        V2 = Vel[indi,indj,indk,1]
                                        W2 = Vel[indi,indj,indk,2]

                                        # BUILD THE WEIGHTS
                                        wx = 1.0
                                        wy = 1.0
                                        wz = 1.0

                                        for o in range(1-1, pp):
                                            indi = (i)*pp+o
                                            indj = (j)*pp+mm
                                            indk = (k)*pp+nn

                                            Xo = Xg[indi,indj,indk]
                                            if(o != ll):
                                                wx = wx*(X1-Xo)/(X2-Xo)

                                        for o in range(1-1, pp):
                                            indi = (i)*pp+ll
                                            indj = (j)*pp+o
                                            indk = (k)*pp+nn

                                            Yo = Yg[indi,indj,indk]
                                            if(o != mm):
                                                wy = wy*(Y1-Yo)/(Y2-Yo)

                                        for o in range(1-1, pp):
                                            indi = (i)*pp+ll
                                            indj = (j)*pp+mm
                                            indk = (k)*pp+o

                                            Zo = Zg[indi,indj,indk]
                                            if(o != nn):
                                                wz = wz*(Z1-Zo)/(Z2-Zo)

                                        # ADD THEM UP
                                        indi = (i)*pp+l
                                        indj = (j)*pp+m
                                        indk = (k)*pp+n

                                        Vel2[indi,indj,indk,0] = Vel2[indi,indj,indk,0] + wx*wy*wz*U2
                                        Vel2[indi,indj,indk,1] = Vel2[indi,indj,indk,1] + wx*wy*wz*V2
                                        Vel2[indi,indj,indk,2] = Vel2[indi,indj,indk,2] + wx*wy*wz*W2

    # WRITE TO FILE
    file = open(output_filename,"w")
    wstr = "%i\n" % nDOF
    file.write(wstr)
    for k in range(Nk1-1, Nk2):
        for j in range(Nj1-1, Nj2):
            for i in range(Ni1-1, Ni2):
                wstr = " %18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % (Xs[i,j,k],Ys[i,j,k],Zs[i,j,k],Vel2[i,j,k,0],Vel2[i,j,k,1],Vel2[i,j,k,2])
                file.write(wstr)
    file.close()
    return

# # ================================================================
# # check that it works
# # ================================================================
# # data_dir = "philip_outputs/"
# # filename="1procs/coord_check_%i_elements_p%i-proc_0.txt" % (nElements_per_direction,poly_degree)
# # filename="8procs/assembled_coords.txt"
# data1 = np.loadtxt("dummy_delete_me.txt",skiprows=0,dtype=np.float64)
# data2 = np.loadtxt("velocity_gl_nodes.fld",skiprows=0,dtype=np.float64)
# file = open("diff_gl_nodes.dat","w")
# for i in range(0,nDOF):
#     check = data1[i,:]
#     ref = data2[i,:]
#     err = np.linalg.norm(check-ref)
#     if(err > 2.0e-15):
#         err_msg = "%i %18.16e \n" % (i,err)
#         file.write(err_msg)
# file.close()
# # ================================================================