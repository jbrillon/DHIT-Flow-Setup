import numpy as np
import sys; sys.path.append("submodules/quickplotlib/lib"); import quickplotlib as qp
# import sys; sys.path.append("/Users/Julien/Python/quickplotlib/lib"); import quickplotlib as qp

z_cut=3.141592653589793116

velocity_field = np.loadtxt("dofs048_p5_velocity/velocity_equidistant_nodes.fld",skiprows=0,usecols=(0,1,2,3,4,5),dtype=np.float64)
# velocity_field = np.loadtxt("dofs128_p3_velocity/velocity_equidistant_nodes.fld",skiprows=0,usecols=(0,1,2,3,4,5),dtype=np.float64)

row_indices_of_plane = np.where(np.equal(velocity_field[:,2],z_cut))[0]
npoints = row_indices_of_plane.size
x = np.zeros(npoints,dtype=np.float64)
y = np.zeros(npoints,dtype=np.float64)
u = np.zeros(npoints,dtype=np.float64)
v = np.zeros(npoints,dtype=np.float64)
for i in range(npoints):
    index = row_indices_of_plane[i]
    x[i] = velocity_field[index,0]
    y[i] = velocity_field[index,1]
    u[i] = velocity_field[index,3]
    v[i] = velocity_field[index,4]

qp.plotfield(xdata=x,ydata=y,udata=u,vdata=v,ylabel="$y$",xlabel="$x$",
    # title_label="DHIT Initial Velocity Field at $z=\\pi$\n $24^{3}$ DOF (P5, $N_{el}=4^{3}$)",
    fig_directory=".",figure_filename="velocity_field_48_init",
    xlimits=[0.0,2.0*np.pi],ylimits=[0.0,2.0*np.pi])
