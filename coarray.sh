#!/bin/bash
clear
echo "Compiling Co-Array Program"
ifort -coarray=share -c jacobiCAF.f90 coarrayjacobi.f90
ifort -coarray=share jacobiCAF.o coarrayjacobi.o -o coarrayjacobi
echo "Compilation Complete"
echo "Creating 16 images..."
export FOR_COARRAY_NUM_IMAGES=16 
echo "Running Program"
./coarrayjacobi
