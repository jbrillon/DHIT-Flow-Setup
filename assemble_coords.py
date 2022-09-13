# test code for assembling coordinates outputted by PHiLiP

tempfiles=[]

num_procs = 4
subdir = "./philip_outputs/%iprocs" % num_procs
for i in range(0,num_procs):
    filename = subdir+"/"+"coord_check-proc_%i.txt" % i
    tempfiles.append(filename)

# tempfiles is an ordered list of temp files (open for reading)
fout = open(subdir+"/"+"bigfile.txt", "w")
for file_index,tempfile in enumerate(tempfiles):
    fin = open(tempfile,"r")
    i=0
    for line in fin:
        if(i==0 and file_index==0):
            fout.write(line) # write first line
        if(i>0): # skip first line
            fout.write(line)
        i += 1
    fin.close()

fout.close()