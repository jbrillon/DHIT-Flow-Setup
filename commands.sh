# (1) running the DHIT code
# for Linux:
# f77 box.for -o main.exe
# for MAC OS X
gfortran -malign-double box.for -o main.exe
./main.exe

# # (2) running the converter.m
# PATH_TO_MATLAB_SCRIPTS="/home/julien/Codes/DHIT/example"
# matlab -nodisplay -nosplash -nodesktop -r "run('${PATH_TO_MATLAB_SCRIPTS}/converter.m');exit;" | tail -n +11