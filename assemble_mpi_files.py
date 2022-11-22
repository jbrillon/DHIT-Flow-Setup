# get padded mpi rank string
def get_padded_mpi_rank_string(mpi_rank):
    padding_length = 5
    mpi_rank_string = '%i' % mpi_rank
    mpi_rank_string = mpi_rank_string.zfill(padding_length)
    return mpi_rank_string

# assemble mpi files
def assemble_mpi_files(path_and_prefix,num_procs):
    filename=path_and_prefix+".dat"
    fout = open(filename, "w")
    for i in range(0,num_procs):
        mpi_rank_string = get_padded_mpi_rank_string(i)
        tempfile = path_and_prefix+"-"+mpi_rank_string+".dat"
        fin = open(tempfile,"r")
        for line in fin:
            fout.write(line)
        fin.close()
    fout.close()