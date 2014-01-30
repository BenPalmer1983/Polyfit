#!/bin/bash
gfortran -o polyfit.x src/kinds.f90 src/initialise.f90 src/maths.f90 src/input.f90 src/polyfit.f90 -L/usr/lib64 -llapack -lblas
sleep 1 
rm *.mod