include ../../example.mk
include ../../common.mk

OBJ = main.o

nbody:

%.o: %.cu
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

nbody: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: nbody

run: nbody
	mpirun --oversubscribe -np 1 ./nbody

.PHONY: clean all run

clean:
	rm -f *.o *~ core nbody

