#compiling VlasovAF

FC = gfortran
PY = python3 

FCFLAGS = -ffree-line-length-512

OBJ = parameters.o advection.o poisson.o output.o main.o

all: main
	@echo "Compilation completed"

main: $(OBJ)
	$(FC) -o $@ $(OBJ)

advection.o: advection.F90
	$(FC) $(FCFLAGS) -c advection.F90

poisson.o: poisson.F90
	$(FC) $(FCFLAGS) -c poisson.F90

main.o: main.F90
	$(FC) $(FCFLAGS) -c main.F90

output.o: output.F90
	$(FC) $(FCFLAGS) -c output.F90

parameters.o: parameters.F90
	$(FC) $(FCFLAGS) -c parameters.F90
	
clean:
	rm *.o *.mod *.dat main
	
plot:
	$(PY) plotResults.py 
	$(PY) plotDistribution.py
	$(PY) plotSlices.py
