#!/bin/bash
gfortran -fno-range-check ./MT/mt19937.f90 mc.f90
rm ./output/coordinate.xyz
./a.out
rm a.out *mod
