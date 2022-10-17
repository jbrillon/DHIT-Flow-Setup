import numpy as np
import pyvista as pv
import meshio
from write_vtk_file import write_vtk_file_uniform_cube

'''
 NOTE: this is called after having ran prep_vel_for_spec.py so
 that were working with the "averaged" (at element faces) velocity field
'''

###############################################################################
# ``mesh_g["gradient"]`` is an ``N`` by 9 NumPy array of the gradients, so we
# could make a dictionary of NumPy arrays of the gradients like:
def gradients_to_dict(arr):
    """A helper method to label the gradients into a dictionary."""
    keys = np.array(
        ["du/dx", "du/dy", "du/dz", "dv/dx", "dv/dy", "dv/dz", "dw/dx", "dw/dy", "dw/dz"]
    )
    keys = keys.reshape((3, 3))[:, : arr.shape[1]].ravel()
    return dict(zip(keys, arr.T))
###############################################################################

###############################################################################
def get_taylor_micro_scale(rms_velocity,du_dx,dv_dy,dw_dz):
    den = np.mean(du_dx**2.0 + dv_dy**2.0 + dw_dz**2.0)
    taylor_micro_scale = np.sqrt(rms_velocity/den)
    return taylor_micro_scale
###############################################################################

def generate_pyvista_mesh_object(coordinates,velocities):
    # write the vtk file
    write_vtk_file_uniform_cube(coordinates,velocities,filename="solution.vtk")
    # read in vtk file to get the mesh object
    mesh = pv.read("solution.vtk",file_format="vtk")
    # mesh.plot(show_edges=True)
    return mesh

def gradients_of_pyvista_mesh_object(mesh):
    #----------------------------------------------------------------------
    # ---- Gradients of velocity field using PyVista
    #----------------------------------------------------------------------
    # Reference: https://docs.pyvista.org/examples/01-filter/gradients.html
    '''
    Compute the gradients of the ``velocity`` vector 
    field in the point data of that mesh. 
    This is as simple as calling
    :func:`pyvista.DataSetFilters.compute_derivative`.

    TO DO: Generate Q-criterion of TGV using this
    # .. note:: You can also use :func:`pyvista.DataSetFilters.compute_derivative` for
    #   computing other derivative based quantities, such as divergence, vorticity,
    #   and Q-criterion. See function documentation for options.
    '''
    mesh_g = mesh.compute_derivative(scalars="velocity",qcriterion=True,vorticity=True)
    mesh_g["gradient"]
    gradients = gradients_to_dict(mesh_g["gradient"])

    ###############################################################################
    # add all of those components as individual arrays back to the mesh object:
    mesh_g.point_data.update(gradients)
    ###############################################################################

    return mesh_g,gradients

def visualize_turbulent_flow_from_mesh_and_gradient_objects(mesh,mesh_g,gradients):
    # TO DO: Clean this up

    # keys = np.array(list(gradients.keys())).reshape(3, 3)

    # p = pv.Plotter(shape=keys.shape)
    # for i in range(keys.shape[0]):
    #     for j in range(keys.shape[1]):
    #         name = keys[i, j]
    #         p.subplot(i, j)
    #         p.add_mesh(mesh_g.contour(scalars=name), scalars=name, opacity=0.75)
    #         p.add_mesh(mesh_g.outline(), color="k")
    # p.link_views()
    # p.view_isometric()
    # p.show()

    p = pv.Plotter()
    # p.add_mesh(mesh_g.contour(scalars="qcriterion"), scalars="qcriterion", opacity=1.0)
    p.add_mesh(mesh_g.contour(scalars="vorticity"), scalars="vorticity", opacity=1.0)
    p.add_mesh(mesh_g.outline(), color="k")
    p.view_isometric()
    p.show()

    # to plot both quantities
    quantities=["vorticity","qcriterion"]
    n_quantities = len(quantities)
    p = pv.Plotter(shape=(1,n_quantities))
    for i in range(0,n_quantities):
        name = quantities[i]
        p.subplot(0,i)
        p.add_mesh(mesh_g.contour(scalars=name), scalars=name, opacity=1.0)
        p.add_mesh(mesh_g.outline(), color="k")

    p.link_views()
    p.view_isometric()
    p.show()

    # mesh_g.plot(scalars="vorticity")

    # ###############################################################################
    # # And there you have it, the gradients for a vector field! We could also do
    # # this for a scalar  field like for the ``scalars`` field in the given dataset.
    # mesh_g = mesh.compute_derivative(scalars="scalars")

    # gradients = gradients_to_dict(mesh_g["gradient"])
    # gradients

    # ###############################################################################

    # mesh_g.point_data.update(gradients)

    # keys = np.array(list(gradients.keys())).reshape(1, 3)

    # p = pv.Plotter(shape=keys.shape)

    # for i in range(keys.shape[0]):
    #     for j in range(keys.shape[1]):
    #         name = keys[i, j]
    #         p.subplot(i, j)
    #         p.add_mesh(mesh_g.contour(scalars=name), scalars=name, opacity=0.75)
    #         p.add_mesh(mesh_g.outline(), color="k")
    # p.link_views()
    # p.view_isometric()
    # p.show()
    return

def visualize_turbulent_flow_from_coordinates_and_velocity(coordinates,velocities):
    mesh = generate_pyvista_mesh_object(coordinates,velocities)
    mesh_g,gradients = gradients_of_pyvista_mesh_object(mesh)
    visualize_turbulent_flow_from_mesh_and_gradient_objects(mesh,mesh_g,gradients)
    return

def compute_turbulent_quantities(coordinates,velocities,
    return_rms_velocity=True,
    return_tke=True,
    return_taylor_micro_scale=True):
    #======================================================================
    # ---- Compute flow quantities
    #======================================================================
    reduced_nDOF = int(coordinates.shape[0])
    # velocity magnitude squared at each point
    velocity_magnitude_squared = np.zeros(reduced_nDOF,dtype=np.float64)
    for i in range(0,reduced_nDOF):
        for d in range(0,3):
            velocity_magnitude_squared[i] += velocities[i,d]*velocities[i,d]

    # mean velocity magnitude squared
    mean_vel_mag_squared = np.mean(velocity_magnitude_squared)

    # turbulent kinetic energy (TKE)
    turbulent_kinetic_energy = 0.5*mean_vel_mag_squared

    # r.m.s. velocity
    rms_velocity = np.sqrt(mean_vel_mag_squared/3.0)

    #----------------------------------------------------------------------
    # ---- Gradients of velocity field using PyVista
    #----------------------------------------------------------------------
    # Reference: https://docs.pyvista.org/examples/01-filter/gradients.html
    mesh = generate_pyvista_mesh_object(coordinates,velocities)
    mesh_g,gradients = gradients_of_pyvista_mesh_object(mesh)

    # extract gradients
    du_dx = np.array(gradients["du/dx"],dtype=np.float64)
    dv_dy = np.array(gradients["dv/dy"],dtype=np.float64)
    dw_dz = np.array(gradients["dw/dz"],dtype=np.float64)

    # Compute Taylor micro-scale
    taylor_micro_scale = get_taylor_micro_scale(rms_velocity,du_dx,dv_dy,dw_dz)

    # Decide what to return
    turbulent_quantities = []
    if(return_rms_velocity==True):
        turbulent_quantities.append(rms_velocity)
    if(return_tke==True):
        turbulent_quantities.append(turbulent_kinetic_energy)
    if(return_taylor_micro_scale==True):
        turbulent_quantities.append(taylor_micro_scale)

    return turbulent_quantities

# # Filename of averaged velocity field to read in
# averaged_velocity_field_filename = "velocity.fld"
# coordinates = np.loadtxt(averaged_velocity_field_filename,skiprows=0,usecols=(0,1,2),dtype=np.float64)
# velocities = np.loadtxt(averaged_velocity_field_filename,skiprows=0,usecols=(3,4,5),dtype=np.float64)

# # Get turbulent quantities
# rms_velocity,turbulent_kinetic_energy,taylor_micro_scale = compute_turbulent_quantities(coordinates,velocities)


