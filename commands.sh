# Script for (1) running the DHIT code
# -- Linux:
# f77 box.for -o main.exe
# -- MAC OS X
gfortran -malign-double -ffixed-line-length-144 box.for -o main.exe
# echo -e "box.in" | ./main.exe
echo -e "box.in\n 1" | ./main.exe # generate the expected coordinates input
echo -e "box.in\n 0" | ./main.exe # run code as normal