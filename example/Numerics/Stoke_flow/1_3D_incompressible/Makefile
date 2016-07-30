include ../../../example.mk

CC=mpic++

LDIR =

OBJ_EIGEN = main_eigen.o
OBJ_PETSC = main_petsc.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

all: stokes_3d_eigen stokes_3d_petsc

stokes_3d_eigen: $(OBJ_EIGEN)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

stokes_3d_petsc: $(OBJ_PETSC)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)


.PHONY: clean all

clean:
	rm -f *.o *~ core stokes_3d_eigen stokes_3d_petsc

