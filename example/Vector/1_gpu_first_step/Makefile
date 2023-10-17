include ../../common.mk
include ../../example.mk

OBJ = main.o

gpu_fstep:

%.o: %.cu
	$(CUDA_CC) -O3 -g $(CUDA_OPTIONS) -c --std=c++14 -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

gpu_fstep: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: gpu_fstep

run: gpu_fstep
	mpirun --oversubscribe -np 2 ./gpu_fstep

.PHONY: clean all run

clean:
	rm -f *.o *~ core gpu_fstep

