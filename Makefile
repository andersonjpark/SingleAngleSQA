NULIB_LOC = /Users/jspark971/external/NuLib
HDF5_LOC  = /Users/jspark971/external/hdf5-1.10.5/hdf5

INCLUDE = -I${HDF5_LOC}/include

LIBRARY = ${NULIB_LOC}/src/nulib.a -L${HDF5_LOC}/lib -lhdf5 -lhdf5_fortran -lhdf5_cpp -L/usr/local/Cellar/gcc/9.2.0_3/lib/gcc/9 -lgfortran -lomp -lstdc++

COMP = g++ -Xpreprocessor -fopenmp -g -std=gnu++11 -O0 -Wall

sqa.x: src/sqa.cpp src/*.h src/*.cpp
	cp NuLib_make.inc NuLib/make.inc
	cp requested_interactions.inc NuLib/src/requested_interactions.inc
	make -C NuLib
	rm -f sqa.x sqa.tar.gz
	tar --exclude=*.txt --exclude=*.lum --exclude=*.cyl -cvzf sqa.tar.gz ./*
	${COMP} src/sqa.cpp -o sqa.x ${INCLUDE} ${LIBRARY}

clean:
	rm -f sqa.x
	$(MAKE) -C NuLib clean
