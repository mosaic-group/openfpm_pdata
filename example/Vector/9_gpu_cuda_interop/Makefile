include ../../example.mk
include ../../common.mk

OBJ = main.o

gpu_interop:

%.o: %.cu
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

gpu_interop: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: gpu_interop

run: gpu_interop
	mpirun --oversubscribe -np 2 ./gpu_interop

.PHONY: clean all run

clean:
	rm -f *.o *~ core gpu_interop

