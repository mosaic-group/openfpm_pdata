include ../../common.mk
include ../../example.mk

OBJ = main.o

gpu_interop:

%.o: %.cu
	$(CUDA_CC) -O3 $(OPT) $(CUDA_OPTIONS) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

gpu_interop: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: gpu_interop

run: gpu_interop
	mpirun --oversubscribe -np 2 ./gpu_interop

.PHONY: clean all run

clean:
	rm -f *.o *~ core gpu_interop

