#! /bin/bash
#
gfortran -c -Wall r8lib.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8lib.o ~/proj/particledynamics1D/source/tools/r8lib.o
#
echo "Normal end of execution."
