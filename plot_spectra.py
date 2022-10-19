import numpy as np
import sys; sys.path.append("submodules/quickplotlib/lib"); import quickplotlib as qp
# import sys; sys.path.append("/Users/Julien/Python/quickplotlib/lib"); import quickplotlib as qp

x=[];y=[];labels=[];

# Original spectra file
spectra = np.loadtxt("energy.prf",skiprows=1,dtype=np.float64)
x.append(spectra[:,0])
y.append(spectra[:,1])
labels.append("$E(k)_{1}\\longrightarrow (u,v,w)$")

spectra = np.loadtxt("cbc_experiment_table3.txt",skiprows=3,usecols=(0,1),dtype=np.float64)
# Misra and Lund non-dimensionalization:
M = 5.08 # [cm] (mesh size from experiment)
u_rms = 22.2 # [cm/s] (rms velocity from experiment)
L_ref = 11.0*M/(2.0*np.pi) # cm
U_ref = np.sqrt(3.0/2.0)*u_rms # cm/s
energy_ref = U_ref*U_ref*L_ref # cm3/s2
# -- nondimensionalize experiment values
spectra[:,0] *= L_ref # non-dimensionalize wavenumber
spectra[:,1] /= energy_ref # non-dimensionalize energy
# add to plot
x.append(spectra[:,0])
y.append(spectra[:,1])
labels.append("CBC Experiment")

# Generate spec.tec before running the commented code below:
# spectra = np.loadtxt("spec.tec",dtype=np.float64)
# x.append(spectra[:,0])
# y.append(spectra[:,1])
# labels.append("$(u,v,w)\\longrightarrow E(k)_{2}$")

qp.plotfxn(xdata=x,ydata=y,xlabel="$k$",ylabel="$E(k)$",
    title_label="DHIT Initialization Check\n P5, $N_{el}=4^{3}$ ($24^{3}$ DOF)",
    fig_directory=".",figure_filename="spectra",log_axes="both",figure_filetype="pdf",
    #xlimits=[1e0,1e2],ylimits=[1e-4,1e-1],
    markers=True,legend_on=True,legend_labels_tex=labels)



# fix quickplotlib so that if no legend_labels are provided, it either prints a message or automatically sets legend_on==False