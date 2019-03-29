#!/bin/csh
MSTL_LOC = /home/srichers/software/mstl
NULIB_LOC = /home/srichers/software/NuLib
HDF5_LOC = /usr/local/hdf5-1.8.20_gnu7.3.0

INCLUDE  = -I${MSTL_LOC}/include/ -I${MSTL_LOC}/include/physics/ -I${MSTL_LOC}/include/math2/ -I${MSTL_LOC}/include/math2/algebra -I${MSTL_LOC}/include/math2/analysis -I${MSTL_LOC}/include/math2/data -I${MSTL_LOC}/include/math2/geometry -I${MSTL_LOC}/include/math2/group -I${MSTL_LOC}/include/math2/spline -I"${MSTL_LOC}/include/math2/probability and statistics" -I${HDF5_LOC}/include

LIBRARY = ${NULIB_LOC}/src/nulib.a -L${HDF5_LOC}/lib -lhdf5 -lhdf5_fortran -lhdf5_cpp -lgfortran -L${MSTL_LOC}/lib/ -lmstl.g++

WHICHCODE = sqa2.psma

COMP = g++ -fopenmp -g -std=gnu++11 -O2

${WHICHCODE}.o: ${WHICHCODE}.cpp
	rm -f sqa2.psma.x sqa.tar.gz
	tar --exclude=*.txt --exclude=*.lum --exclude=*.cyl -cvzf sqa.tar.gz ./*
	${COMP} ${WHICHCODE}.cpp -o ${WHICHCODE}.x ${INCLUDE} ${LIBRARY}
