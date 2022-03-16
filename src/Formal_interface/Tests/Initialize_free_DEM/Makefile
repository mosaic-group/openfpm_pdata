include paths.mk

CC=mpic++

all: DEM_Test

%.o: %.cpp %.hpp
	$(CC) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

DEM_Test: DEM_Test.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: DEM_Test
	mpirun --oversubscribe -np 2 ./DEM_Test;

debug: DEM_Test
	mpirun -n 2 xterm -hold -e gdb -ex run  ./DEM_Test;

.PHONY: clean all run

clean:
	rm -f *.o *~ core DEM_Test *.pvtp *.vtp *.vtk

