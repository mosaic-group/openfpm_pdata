include ../../example.mk

### internally the example disable with the preprocessor its code if not compiled with nvcc 
CUDA_CC=
CUDA_CC_LINK=
ifdef CUDA_ON_CPU
        CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
        INCLUDE_PATH_NVCC=
        CUDA_CC_LINK=mpic++
else
        ifeq (, $(shell which nvcc))
                CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
                INCLUDE_PATH_NVCC=
                CUDA_CC_LINK=mpic++
        else
                CUDA_CC=nvcc -ccbin=mpic++
                CUDA_CC_LINK=nvcc -ccbin=mpic++
        endif
endif

gray_scott_sparse_gpu_test: OPT += -DTEST_RUN
gray_scott_sparse_gpu_test: gray_scott_sparse_gpu

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cu
	$(CUDA_CC) $(OPT) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH_NVCC)

gray_scott_sparse_gpu: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(LIBS_PATH) $(LIBS)

all: gray_scott_sparse_gpu

run: gray_scott_sparse_gpu_test
	mpirun -np 4 ./gray_scott_sparse_gpu

.PHONY: clean all run

clean:
	rm -f *.o *~ core gray_scott_sparse_gpu

