# Script for (1) running the DHIT code
# -- Linux:
# f77 box.for -o main.exe
# -- MAC OS X
gfortran -malign-double -ffixed-line-length-144 box.for -o main.exe
./main.exe