# get padded mpi rank string
def get_padded_mpi_rank_string(mpi_rank):
    padding_length = 5
    mpi_rank_string = '%i' % i_proc
    mpi_rank_string.rjust(padding_length, "0")
    return mpi_rank_string

num_procs = 4
subdir = "./philip_outputs"
prefix = "velocity"

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