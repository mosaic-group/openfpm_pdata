include ../../example.mk
include ../../common.mk

gray_scott_sparse_gpu_test: OPT += -DTEST_RUN
gray_scott_sparse_gpu_test: gray_scott_sparse_gpu

OBJ = main.o

%.o: %.cu
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)

gray_scott_sparse_gpu: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(LIBS_PATH) $(LIBS_SELECT)

all: gray_scott_sparse_gpu

run: gray_scott_sparse_gpu_test
	mpirun --oversubscribe -np 4 ./gray_scott_sparse_gpu

.PHONY: clean all run

clean:
	rm -f *.o *~ core gray_scott_sparse_gpu

