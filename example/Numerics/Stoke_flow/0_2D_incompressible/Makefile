include ../../../common.mk
include ../../../example.mk

OBJ_EIGEN = main_eigen.o
OBJ_PETSC = main_petsc.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

all: stokes_2d_eigen stokes_2d_petsc

stokes_2d_eigen: $(OBJ_EIGEN)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

stokes_2d_petsc: $(OBJ_PETSC)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all
	mpirun --oversubscribe -np 3 ./stokes_2d_eigen && mpirun --oversubscribe -np 3 ./stokes_2d_petsc
	
.PHONY: clean all

clean:
	rm -f *.o *~ core stokes_2d_eigen stokes_2d_petsc

