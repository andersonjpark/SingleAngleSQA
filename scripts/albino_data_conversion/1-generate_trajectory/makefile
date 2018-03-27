all:
	mkdir -p output
	gfortran -c -fdefault-real-8 main.f90
	gfortran -o main.x main.o
	./main.x

clean:
	rm main.x
	rm main.o
