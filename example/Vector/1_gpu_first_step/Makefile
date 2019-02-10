include ../../example.mk

CUDA_CC=
ifeq (, $(shell which nvcc))
        CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
        INCLUDE_PATH_NVCC=
else
        CUDA_CC=nvcc
endif


CC=mpic++

OBJ = main.o

gpu_fstep:

%.o: %.cu
	$(CUDA_CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

gpu_fstep: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: gpu_fstep

run: gpu_fstep
	mpirun -np 2 ./gpu_fstep

.PHONY: clean all run

clean:
	rm -f *.o *~ core gpu_fstep

