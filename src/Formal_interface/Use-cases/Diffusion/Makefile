include paths.mk

CC=mpic++

all: Dif

%.o: %.cpp %_Algorithm.hpp %_Instance.hpp
	$(CC) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

Dif: Dif.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: Dif
	mpirun --oversubscribe -np 2 ./Dif;

debug: Dif
	mpirun -n 2 xterm -hold -e gdb -ex run  ./Dif;

.PHONY: clean all run

clean:
	rm -f *.o *~ Dif *.pvtp *.vtp *.vtk

