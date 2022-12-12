include paths.mk

CC=mpic++

all: OpenFPM_version FI_parallel_version

%.o: %.cpp %.hpp
	$(CC) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

OpenFPM_version: OpenFPM_version.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

FI_parallel_version: FI_parallel_version.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: OpenFPM_version FI_parallel_version
	mpirun --oversubscribe -np 4 ./OpenFPM_version;
	mpirun --oversubscribe -np 4 ./FI_parallel_version;
	mpirun --oversubscribe -np 2 ./OpenFPM_version;
	mpirun --oversubscribe -np 2 ./FI_parallel_version;
	mpirun --oversubscribe -np 1 ./OpenFPM_version;
	mpirun --oversubscribe -np 1 ./FI_parallel_version;

.PHONY: clean all run

clean:
	rm -f *.o *~ core OpenFPM_version FI_parallel_version  *.pvtp *.vtp *.vtk

