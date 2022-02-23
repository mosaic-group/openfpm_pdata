include paths.mk

CC=mpic++

all: Test1 DEM_Test PSE_Test testing

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

Test1: Test1.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

DEM_Test: DEM_Test.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

PSE_Test: PSE_Test.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

testing: testing.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run_%: %
	mpirun --oversubscribe -np 2 ./$^;

debug_%: %
	mpirun -n 2 xterm -hold -e gdb -ex run  ./$^;

.PHONY: clean all run

clean:
	rm -f *.o *~ core Test1 *.pvtp *.vtp *.vtk

