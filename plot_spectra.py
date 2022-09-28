import numpy as np
import sys; sys.path.append("submodules/quickplotlib/lib"); import quickplotlib as qp
# import sys; sys.path.append("/Users/Julien/Python/quickplotlib/lib"); import quickplotlib as qp

x=[];y=[];labels=[];

spectra = np.loadtxt("energy.prf",skiprows=1,dtype=np.float64)
x.append(spectra[:,0])
y.append(spectra[:,1])
labels.append("input")

spectra = np.loadtxt("spec.tec",dtype=np.float64)
x.append(spectra[:,0])
y.append(spectra[:,1])
labels.append("output")

qp.plotfxn(xdata=x,ydata=y,
    fig_directory=".",figure_filename="spectra",log_axes="both",
    xlimits=[],ylimits=[],
    markers=True,legend_on=True,legend_labels_tex=labels)

# fix quickplotlib so that if no legend_labels are provided, it either prints a message or automatically sets legend_on==False