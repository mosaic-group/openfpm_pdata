include ../../common.mk
include ../../example.mk

OBJ = main.o

gray_scott_sparse_gpu_test: OPT += -DTEST_RUN
gray_scott_sparse_gpu_test: gray_scott_sparse_gpu

%.o: %.cu
	$(CUDA_CC) $(OPT) $(CUDA_OPTIONS) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH_NVCC)

gray_scott_sparse_gpu: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(LIBS_PATH) $(LIBS_SELECT)

all: gray_scott_sparse_gpu

run: gray_scott_sparse_gpu_test
	mpirun --oversubscribe -np 4 ./gray_scott_sparse_gpu

.PHONY: clean all run

clean:
	rm -f *.o *~ core gray_scott_sparse_gpu

