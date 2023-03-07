import numpy as np
import sys; sys.path.append("../submodules/TurboGenPY"); import TurboGenPY_HighOrderFEM as turboFEM
import sys; sys.path.append("../submodules/PHiLiP-Post-Processing/src/tools"); import generate_spectra_files as gsf
import sys; sys.path.append("../submodules/quickplotlib/lib"); import quickplotlib as qp
import convert_equidistant_to_gauss_lobatto_nodes as eq2gll
import flow_parameter_calc as fpc
#-----------------------------------------------------
def get_DOF_vars(nElements_per_direction,poly_degree):
    nElements = nElements_per_direction*nElements_per_direction*nElements_per_direction
    nQuadPoints_per_element = poly_degree + 1
    nQuadPoints = nQuadPoints_per_element*nElements_per_direction
    nDOF = nQuadPoints*nQuadPoints*nQuadPoints
    reduced_nQuadPoints = nElements_per_direction*nQuadPoints_per_element - (nElements_per_direction-1)
    reduced_nDOF = reduced_nQuadPoints*reduced_nQuadPoints*reduced_nQuadPoints
    return [nElements,nQuadPoints_per_element,nQuadPoints,nDOF,reduced_nQuadPoints,reduced_nDOF]
#-----------------------------------------------------
def get_reduced_nDOF_and_nQuadPoints(nElements_per_direction,poly_degree):
    reduced_nQuadPoints = get_DOF_vars(nElements_per_direction,poly_degree)[-2]
    reduced_nDOF = get_DOF_vars(nElements_per_direction,poly_degree)[-1]
    return [reduced_nQuadPoints,reduced_nDOF]
#-----------------------------------------------------
#*****************************************************
#               GENERAL VARIABLES
#*****************************************************
#=====================================================
# INPUTS
#=====================================================
nValues_per_row = 6 # 3 coords, 3 vel components
nElements_per_direction = 32 # number of elements in each direction/dimension
poly_degree = 3 # solution polynomial degree
num_procs = 64 # 1024 # number of mpi processors to run PHiLiP with
output_directory = "./"
spectra_name='ml'
# Get DOF variables
nElements,nQuadPoints_per_element,nQuadPoints,nDOF,reduced_nQuadPoints,reduced_nDOF = get_DOF_vars(nElements_per_direction,poly_degree)
#-----------------------------------------------------
#=====================================================
# Generate iso turb files for PHiLiP code
#=====================================================

# generate the velocity flow field
turboFEM.generate_isotropic_turbulence_high_order_fem(
    nElements_per_direction,
    poly_degree,
    output_filename=(output_directory+"velocity_equidistant_nodes.fld"),
    number_of_modes=5000, # suggested in TurboGenPY paper
    spectra_name=spectra_name
    )

# verify that resulting spectra from generated flow field is what we expect
file_without_extension = (output_directory+"velocity_equidistant_nodes")
file_extension = "fld"
# equidistant_velocity_file = file_without_extension+"."+file_extension
gsf.generate_spectra_file_from_flow_field_file(
    file_without_extension,
    file_extension,
    n_skiprows=1,
    use_TurboGenPy=True)

# plot the spectra
x = [];y=[];labels=[];title_label="Spectra Comparison";

# load data
cbcspec = np.loadtxt("../submodules/TurboGenPY/cbc_spectrum.txt")
if(spectra_name=="ml"):
    # Misra and Lund non-dimensionalization:
    M = 5.08 # [cm] (mesh size from experiment)
    u_rms = 22.2 # [cm/s] (rms velocity from experiment)
    L_ref = 11.0*M/(2.0*np.pi) # cm
    U_ref = np.sqrt(3.0/2.0)*u_rms # cm/s
    energy_ref = U_ref*U_ref*L_ref # cm3/s2
    # -- nondimensionalize experiment values
    kcbc=cbcspec[:,0]*L_ref
    ecbc=cbcspec[:,1]/energy_ref
    x.append(kcbc);y.append(ecbc);labels.append("Misra and Lund")
elif(spectra_name=="cbc"):
    kcbc=cbcspec[:,0]*100 # [1/m]
    ecbc=cbcspec[:,1]*1e-6 # [m3/s/s]
    x.append(kcbc);y.append(ecbc);labels.append("CBC")
else:
    print("ERROR: Invalid spectra name. Aborting...")
    exit()

x_file,y_file = np.loadtxt(output_directory+"velocity_equidistant_nodes_spectra.fld",skiprows=1,unpack=True)
x.append(x_file);y.append(y_file);labels.append("Generated")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
            title_label=title_label,
            fig_directory=".",figure_filename="spectra",log_axes="both",figure_filetype="pdf",
            nlegendcols=1,
            # xlimits=[2e0,1e2],ylimits=[1e-7,3e-2],
            markers=False,legend_on=True,legend_labels_tex=labels,
            which_lines_black=[0],
            transparent_legend=True,legend_border_on=False,grid_lines_on=False
            # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
            )
exit()

# convert to GLL nodes
eq2gll.convert_equidistant_to_gauss_lobatto_nodes(
    (output_directory+"velocity_equidistant_nodes.fld"),
    nElements_per_direction,
    nQuadPoints_per_element,
    nValues_per_row,
    poly_degree,
    nDOF,
    output_filename=(output_directory+"velocity_gl_nodes.fld"),
    test_reading=False # set a true if testing that this function works
    )

# run the ML scaling calc
fpc.determine_flow_parameters(spectra_name=spectra_name,filename=(output_directory+"parameters_for_dhit_setup.txt"))
