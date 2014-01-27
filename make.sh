#!/bin/bash
gfortran -o polyfit.x src/kinds.f90 src/maths.f90 src/polyfit.f90 
sleep 1 
rm *.mod