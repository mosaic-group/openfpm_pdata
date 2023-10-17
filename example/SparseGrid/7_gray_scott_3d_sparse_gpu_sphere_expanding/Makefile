include ../../common.mk
include ../../example.mk

OBJ = main.o

%.o: %.cu
	$(CUDA_CC) $(OPT) $(CUDA_OPTIONS) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH_NVCC)

gray_scott_sparse_gpu: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(LIBS_PATH) $(LIBS_SELECT)

all: gray_scott_sparse_gpu

run: gray_scott_sparse_gpu
	mpirun --oversubscribe -np 4 ./gray_scott_sparse_gpu

.PHONY: clean all run

clean:
	rm -f *.o *~ core gray_scott_sparse_gpu

