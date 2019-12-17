include ../../../example.mk

CC=mpic++

LDIR =

OBJ_MODIF = main_modif.o
OBJ_EIGEN = main_eigen.o
OBJ_PETSC = main_petsc.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

all: stokes_2d_eigen stokes_2d_petsc main_modif

stokes_2d_eigen: $(OBJ_EIGEN)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

stokes_2d_petsc: $(OBJ_PETSC)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

main_modif: $(OBJ_MODIF)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all
	mpirun -np 3 ./stokes_2d_eigen && mpirun -np 3 ./stokes_2d_petsc

.PHONY: clean all

clean:
	rm -f *.o *~ core stokes_2d_eigen stokes_2d_petsc main_modif

