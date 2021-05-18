include ../../../example.mk

CC=mpic++

LDIR =

OBJ_EIGEN = main_eigen.o
OBJ_PETSC = main_petsc.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

all: stokes_3d_eigen stokes_3d_petsc

stokes_3d_eigen: $(OBJ_EIGEN)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

stokes_3d_petsc: $(OBJ_PETSC)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all
	mpirun --oversubscribe -np 3 ./stokes_3d_eigen && mpirun --oversubscribe -np 3 ./stokes_3d_petsc
	
.PHONY: clean all run

clean:
	rm -f *.o *~ core stokes_3d_eigen stokes_3d_petsc

