# get padded mpi rank string
def get_padded_mpi_rank_string(mpi_rank):
    padding_length = 5
    mpi_rank_string = '%i' % mpi_rank
    mpi_rank_string = mpi_rank_string.zfill(padding_length)
    return mpi_rank_string

num_procs = 4
subdir = "/home/julien/Codes/2022-09-01/PHiLiP/build_release/tests/integration_tests_control_files/decaying_homogeneous_isotropic_turbulence"
prefix = "velocity-0"

filename=subdir+"/"+prefix+".dat"
fout = open(filename, "w")
for i in range(0,num_procs):
    mpi_rank_string = get_padded_mpi_rank_string(i)
    tempfile = subdir+"/"+prefix+"-"+mpi_rank_string+".dat"
    fin = open(tempfile,"r")
    for line in fin:
        fout.write(line)
    fin.close()
fout.close()