#! /bin/bash

# compile
rm main.x
rm main.o
gfortran -c -fdefault-real-8 main.f90
gfortran -o main.x main.o

# run
./main.x
