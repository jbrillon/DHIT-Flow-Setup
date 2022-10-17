from includes import *
from var import get_reduced_nDOF_and_nQuadPoints
from taylor_green_vortex_initial_condition import *
from turbulent_flow_3D_data_processor import *

nElements_per_direction = 8
poly_degree = 5
reduced_nQuadPoints, reduced_nDOF = get_reduced_nDOF_and_nQuadPoints(nElements_per_direction,poly_degree)
#----------------------------------------------------------------------
# Generate coordinates
#----------------------------------------------------------------------
coordinates1D = np.linspace(0., 2.0*np.pi, num=reduced_nQuadPoints)
coordinates = np.zeros((reduced_nDOF,3),dtype=np.float64)
p = 0
for k in range(0,reduced_nQuadPoints):
    for j in range(0,reduced_nQuadPoints):
        for i in range(0,reduced_nQuadPoints):
            coordinates[p,0] = coordinates1D[i] # x
            coordinates[p,1] = coordinates1D[j] # y
            coordinates[p,2] = coordinates1D[k] # z
            p += 1

#----------------------------------------------------------------------
# Set velocity field to TGV as a test
#----------------------------------------------------------------------
velocities = np.zeros((reduced_nDOF,3),dtype=np.float64)
for i in range(0,reduced_nDOF):
    velocities[i,:] = tgv_initial_condition_primitive(coordinates[i])[1:4]
#----------------------------------------------------------------------

# Get turbulent quantities
rms_velocity,taylor_micro_scale = \
    compute_turbulent_quantities(
        coordinates,velocities,
        return_rms_velocity=True,
        return_tke=False,
        return_taylor_micro_scale=True)

# Compute exact velocity gradients
du_dx_exact = np.zeros(reduced_nDOF,dtype=np.float64)
dv_dy_exact = np.zeros(reduced_nDOF,dtype=np.float64)
dw_dz_exact = np.zeros(reduced_nDOF,dtype=np.float64)
for i in range(0,reduced_nDOF):
    du_dx_exact[i],dv_dy_exact[i],dw_dz_exact[i] = tgv_initial_condition_velocity_gradients(coordinates[i])

# Compute TMS from exact gradients
taylor_micro_scale_exact = get_taylor_micro_scale(rms_velocity,du_dx_exact,dv_dy_exact,dw_dz_exact)

err = np.abs(taylor_micro_scale_exact - taylor_micro_scale)/taylor_micro_scale_exact
print("Error in computing Taylor micro-scale is: %1.3e" % err)

# print(du_dx)
# print(du_dx_exact)
# sum_gradients = du_dx+dv_dy+dw_dz
# print(sum_gradients)
# sum_gradients_exact = du_dx_exact+dv_dy_exact+dw_dz_exact
# print(sum_gradients_exact)
# print(np.linalg.norm(sum_gradients))
# print(np.linalg.norm(sum_gradients_exact))

# get_taylor_microscale_reynolds_number(rms_velocity,du_dx,dv_dy,dw_dz)

# print(get_taylor_microscale(rms_velocity,du_dx,dv_dy,dw_dz))
# print(get_taylor_microscale(rms_velocity,du_dx_exact,dv_dy_exact,dw_dz_exact))
# error_du_dx = np.linalg.norm(du_dx-du_dx_exact)
# error_dv_dy = np.linalg.norm(dv_dy-dv_dy_exact)
# error_dw_dz = np.linalg.norm(dw_dz-dw_dz_exact)
# tol = 1.0e-13
# if((error_du_dx>tol) or (error_dv_dy>tol) or (error_dw_dz>tol)):
#     print("Error in computed velocity gradients is larger than specified tolerance of %1.1e" % tol)
#     print("Error norms (du_dx,dv_dy,dw_dz): %1.6e %1.6e %1.6e" % (error_du_dx,error_dv_dy,error_dw_dz))
#     print("Aborting...")
#     exit()