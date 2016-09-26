include ../../../example.mk

CC=mpic++

LDIR =

OBJ_EIGEN = main_eigen.o
OBJ_PETSC = main_petsc.o

%.o: %.cpp
	$(CC) -fext-numeric-literals  -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

all: stokes_3d_eigen stokes_3d_petsc

stokes_3d_eigen: $(OBJ_EIGEN)
	$(CC) -fopenmp -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

stokes_3d_petsc: $(OBJ_PETSC)
	$(CC) -fopenmp -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all
	source $$HOME/openfpm_vars; mpirun -np 3 ./stokes_3d_eigen; mpirun -np 3 ./stokes_3d_petsc
	
.PHONY: clean all run

clean:
	rm -f *.o *~ core stokes_3d_eigen stokes_3d_petsc

