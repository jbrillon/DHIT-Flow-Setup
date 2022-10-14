import numpy as np
import pyvista as pv
import meshio
from var import *
from write_vtk_file import write_vtk_file_uniform_cube


def tgv_initial_condition_primitive(pos):
    x,y,z = pos
    density = 1.0
    mach_inf = 0.1
    gamma = 1.4
    # primitive variables
    u = np.sin(x)*np.cos(y)*np.cos(z)
    v = -np.cos(x)*np.sin(y)*np.cos(z)
    w = 0.0
    P = (1.0/(gamma*mach_inf*mach_inf)) + (1.0/16.0)*(np.cos(2.0*x)+np.cos(2.0*y))*(np.cos(2.0*z)+2.0)
    return np.array([density,u,v,w,P])

def primitive_to_conservative(primitive_soln):
    density,u,v,w,P = primitive_soln
    gamma = 1.4
    # to conservative variables:
    xmomentum = density*u; ymomentum = density*v; zmomentum = density*w;
    energy = P/(gamma-1.0) + 0.5*density*(u*u + v*v + w*w)
    return np.array([density,xmomentum,ymomentum,zmomentum,energy])

def tgv_initial_condition(pos):
    primitive_soln = tgv_initial_condition_primitive(pos)
    conservative_soln = primitive_to_conservative(primitive_soln)
    return conservative_soln

'''
 this is called after having ran prep_vel_for_spec.py so
 that were working with the "averaged" (at element faces) velocity field
'''

# Filename of averaged velocity field to read in
averaged_velocity_field_filename = "velocity.fld"
coordinates = np.loadtxt(averaged_velocity_field_filename,skiprows=0,usecols=(0,1,2),dtype=np.float64)
velocities = np.loadtxt(averaged_velocity_field_filename,skiprows=0,usecols=(3,4,5),dtype=np.float64)

#======================================================================
# ---- Compute flow quantities
#======================================================================

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
# Set velocity field to TGV as a test
# for i in range(0,reduced_nDOF):
#     velocities[i,:] = tgv_initial_condition(coordinates[i])[1:4]
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# ---- Gradients of velocity field using PyVista
#----------------------------------------------------------------------
# Reference: https://docs.pyvista.org/examples/01-filter/gradients.html

# write the vtk file
write_vtk_file_uniform_cube(coordinates,velocities,filename="solution.vtk")
# read in vtk file to get the mesh object
mesh = pv.read("solution.vtk",file_format="vtk")
# mesh.plot(show_edges=True)

'''
Now compute the gradients of the ``velocity`` vector 
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


gradients = gradients_to_dict(mesh_g["gradient"])

###############################################################################
# And we can add all of those components as individual arrays back to the mesh
# by:
mesh_g.point_data.update(gradients)

###############################################################################

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

mesh_g.plot(scalars="vorticity")


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

