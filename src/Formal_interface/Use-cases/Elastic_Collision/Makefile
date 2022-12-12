include paths.mk

CC=mpic++

all: EC

%.o: %.cpp %_Algorithm.hpp %_Instance.hpp
	$(CC) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

EC: EC.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: EC
	mpirun --oversubscribe -np 4 ./EC;

debug: EC
	mpirun -n 2 xterm -hold -e gdb -ex run  ./EC;

.PHONY: clean all run

clean:
	rm -f *.o *~ core SPH *.pvtp *.vtp *.vtk

