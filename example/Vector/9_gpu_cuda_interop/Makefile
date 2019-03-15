include ../../example.mk

ifeq (, $(shell which nvcc))
        CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
        INCLUDE_PATH_NVCC=
        CUDA_CC_LINK=mpic++
        CUDA_OPTIONS=
else
        CUDA_CC=nvcc -ccbin=mpic++
        CUDA_CC_LINK=nvcc -ccbin=mpic++
        CUDA_OPTIONS=-use_fast_math  -arch=sm_61 -lineinfo
endif


OBJ = main.o

gpu_interop:

%.o: %.cu
	$(CUDA_CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CUDA_CC) -O3 $(OPT) -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

gpu_interop: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: gpu_interop

run: gpu_interop
	mpirun -np 2 ./gpu_interop

.PHONY: clean all run

clean:
	rm -f *.o *~ core gpu_interop

