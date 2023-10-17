include ../../common.mk
include ../../example.mk

OBJ = main.o

miniBUDE:

%.o: %.cu
	$(CUDA_CC) -O3 $(CUDA_OPTIONS) $(OPT) -c --std=c++14 -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

miniBUDE: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: miniBUDE

run: miniBUDE
	mpirun --oversubscribe -np 1 ./miniBUDE

.PHONY: clean all run

clean:
	rm -f *.o *~ core miniBUDE

