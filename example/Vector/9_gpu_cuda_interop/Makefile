include ../../example.mk

CC=mpic++

OBJ = main.o

gpu_interop:

%.o: %.cu
	nvcc -O3 -g -c -isystem=/home/i-bird/MPI/include --std=c++11 -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

gpu_interop: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: gpu_interop

run: gpu_interop
	mpirun -np 2 ./gpu_interop

.PHONY: clean all run

clean:
	rm -f *.o *~ core gpu_interop

