# system choice
# CONFIG = "intel"
 CONFIG = "gfortran"

ifeq ($(CONFIG),"gfortran")
FC      = gfortran    # compiler
LD      = gfortran -lgfortran     # linker
LDFLAGS =           # linker options
SFLAGS  = -c -mcmodel=large -fbacktrace -fopenmp # standard compiler flags
OFLAGS  = -O2  # optimization compiler flags
DFLAGS  = -O0 -g -ggdb -fbounds-check -Wuninitialized -Wall -Wextra -C # debug compiler flags
# mpi library
MPIINCL =
MPILIB  =
MPI     =
# Openmp library
OMPINCL =
OMPLIB  = -fopenmp
OMP     =
# hdf5 library
HDF5INCL =
HDF5LIB  =
HDF5     =
endif

ifeq ($(CONFIG),"intel")
FC      = ifort     # compiler
LD      = ifort     # linker
LDFLAGS =           # linker options
SFLAGS  = -c -r8 -mcmodel=large -i-dynamic -traceback -fpe0 -openmp  # standard compiler flags
OFLAGS  = -O2 -axW -xW -ip -ipo -unroll=3 # optimization compiler flags
DFLAGS  = -O0 -g -ggdb -check bounds -WB -check uninit -debug all -warn all -C # debug compiler flags
# mpi library
MPIINCL =
MPILIB  =
MPI     =
# Openmp library
OMPINCL =
OMPLIB  = -openmp
OMP     =
# hdf5 library
HDF5INCL =
HDF5LIB  =
HDF5     =
endif

# executable name
EXEC = ../vnu

# put compiler flags together
 CFLAGS = $(SFLAGS) $(OFLAGS) $(OPT) $(MPIINCL) $(HDF5INCL) # optimize
# CFLAGS = $(SFLAGS) $(DFLAGS) $(OPT) $(MPIINCL) $(HDF5INCL) # debug
# CFLAGS = $(SFLAGS) $(OPT) $(MPIINCL) $(HDF5INCL) # standard

# put libraries together
LDLIBS = $(MPILIB) $(MPI) $(HDF5LIB) $(HDF5) $(OMPLIB) $(OMP)

# source directory
outer    = 

OBJS =\
  main.o\
  nud_module.o\
  parameter_module.o\
  read_luminosity.o\
  read_2Davg.o\
  read_2Dtau.o\
  tool_module.o\
  units_module.o\

$(EXEC):\
  $(OBJS); $(LD) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $(EXEC)

clean: 
	rm -f *.o *.mod $(EXEC)

nud_module.o:\
  parameter_module.o\
  tool_module.o\
  units_module.o\
  $(outer)$(@:.o=.f90); $(FC) $(CFLAGS) $(outer)$(@:.o=.f90)

main.o:\
  nud_module.o\
  parameter_module.o\
  units_module.o\
  $(outer)$(@:.o=.f90); $(FC) $(CFLAGS) $(outer)$(@:.o=.f90)

parameter_module.o:\
  units_module.o\
  $(outer)$(@:.o=.f90); $(FC) $(CFLAGS) $(outer)$(@:.o=.f90)

read_luminosity.o:\
  parameter_module.o\
  $(outer)$(@:.o=.f90); $(FC) $(CFLAGS) $(outer)$(@:.o=.f90)

read_2Davg.o:\
  parameter_module.o\
  units_module.o\
  $(outer)$(@:.o=.f90); $(FC) $(CFLAGS) $(outer)$(@:.o=.f90)

read_2Dtau.o:\
  parameter_module.o\
  units_module.o\
  $(outer)$(@:.o=.f90); $(FC) $(CFLAGS) $(outer)$(@:.o=.f90)

tool_module.o:\
  parameter_module.o\
  units_module.o\
  $(outer)$(@:.o=.f90); $(FC) $(CFLAGS) $(outer)$(@:.o=.f90)

units_module.o:\
  $(outer)$(@:.o=.f); $(FC) $(CFLAGS) $(outer)$(@:.o=.f)

