import numpy as np
import sys; sys.path.append("submodules/quickplotlib/lib"); import quickplotlib as qp
# import sys; sys.path.append("/Users/Julien/Python/quickplotlib/lib"); import quickplotlib as qp

x=[];y=[];labels=[];

spectra = np.loadtxt("energy.prf",skiprows=1,dtype=np.float64)
x.append(spectra[:,0])
y.append(spectra[:,1])
labels.append("$E(k)_{1}\\longrightarrow (u,v,w)$")

spectra = np.loadtxt("spec.tec",dtype=np.float64)
x.append(spectra[:,0])
y.append(spectra[:,1])
labels.append("$(u,v,w)\\longrightarrow E(k)_{2}$")

qp.plotfxn(xdata=x,ydata=y,xlabel="$k$",ylabel="$E(k)$",
    title_label="DHIT Initialization Check\n P5, $N_{el}=4^{3}$ ($24^{3}$ DOF)",
    fig_directory=".",figure_filename="spectra",log_axes="both",
    xlimits=[],ylimits=[],
    markers=True,legend_on=True,legend_labels_tex=labels)



# fix quickplotlib so that if no legend_labels are provided, it either prints a message or automatically sets legend_on==False