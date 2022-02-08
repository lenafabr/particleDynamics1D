#! /bin/bash
#
gfortran -c -Wall test_mat.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv test_mat.o ~/proj/particledynamics1D/source/tools/
#
echo "Normal end of execution."
