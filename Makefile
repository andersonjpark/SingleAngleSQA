NULIB_LOC = ./NuLib
HDF5_LOC = /usr/local/hdf5-1.8.20_gnu7.3.0

INCLUDE = -I${HDF5_LOC}/include

LIBRARY = ${NULIB_LOC}/src/nulib.a -L${HDF5_LOC}/lib -lhdf5 -lhdf5_fortran -lhdf5_cpp -lgfortran

COMP = g++ -fopenmp -g -std=gnu++11 -O2 -Wall

sqa.x: src/sqa.cpp src/*.h src/*.cpp
	cp NuLib_make.inc NuLib/make.inc
	make -C NuLib
	rm -f sqa.x sqa.tar.gz
	tar --exclude=*.txt --exclude=*.lum --exclude=*.cyl -cvzf sqa.tar.gz ./*
	${COMP} src/sqa.cpp -o sqa.x ${INCLUDE} ${LIBRARY}

clean:
	rm -f sqa.x
	$(MAKE) -C NuLib clean
