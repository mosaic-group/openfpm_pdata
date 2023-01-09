include paths.mk

CC=mpic++

all: SPH

%.o: %.cpp %_Algorithm.hpp %_Instance.hpp
	$(CC) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

SPH: SPH.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: SPH
	mpirun --oversubscribe -np 4 ./SPH;

debug: SPH
	mpirun -n 2 xterm -hold -e gdb -ex run  ./SPH;

.PHONY: clean all run

clean:
	rm -f *.o *~ core SPH *.pvtp *.vtp *.vtk

