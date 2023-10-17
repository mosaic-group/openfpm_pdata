include ../../common.mk
include ../../example.mk

OBJ = main.o

nbody:

%.o: %.cu
	$(CUDA_CC) -O3 $(CUDA_OPTIONS) $(OPT) -c --std=c++14  -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) -g -O3 $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

nbody: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: nbody

run: nbody
	mpirun --oversubscribe -np 1 ./nbody

.PHONY: clean all run

clean:
	rm -f *.o *~ core nbody

